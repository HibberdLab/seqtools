#!/usr/bin/ruby

#
# Draw  plots to show coverage over introns for different splice variants
#   of a gene
#
# Requires:
#   a set of fastq files OR a sam file
#   a genome
#   a gff file
#   a gene of interest
#
# Chris Boursnell (cmb211@cam.ac.uk)
# v1 created: 12/07/2013
# v3 updated: 17/07/2013
#

# How?
# run bowtie
# load genome
# load gff for gene of interest
# create an array representation for each splice variant
#    use genomic coordinates
#    start will all values as 0
# for each line in the sam file check if it falls completely with in any of the introns
#   for each splice variant
# draw a diagram for each splice variant with blocks for introns and a graph for coverage
#
# TODO: sort the sam files
#       can use to speed up allocation in exons in the gff file (which is also sorted)
#
#

require 'rubygems'
require 'trollop'
require 'rasem'

opts = Trollop::options do
  version "v0.0.1a"
  opt :sam, "Sam file", :required => true, :type => String
  opt :sam2, "Second sam file", :required => true, :type => String
  opt :gff, "GFF file", :required => true, :type => String
  #opt :gene, "Gene name (don't include the .N)", :required => true, :type => String
end
Trollop::die :sam, "must exist" if !File.exist?(opts[:sam]) if opts[:sam]
Trollop::die :sam2, "must exist" if !File.exist?(opts[:sam2]) if opts[:sam2]
Trollop::die :gff, "must exist" if !File.exist?(opts[:gff]) if opts[:gff]

yscale = 2
xscale = 0.2

class Protein
  attr_accessor :chromosome, :start, :stop, :strand, :exons, :name, :coverage, :coverage2

  def initialize(chromosome, start, stop, strand, name)
    chromosome.tr_s!("Chr", "")
    @chromosome = chromosome.to_i
    @start = start.to_i
    @stop = stop.to_i
    @strand = strand
    @name = name
    @exons = Array.new
    @coverage = Array.new(@stop-@start+1, 0)
    @coverage2 = Array.new(@stop-@start+1, 0)
  end

  def addExon(exon)
    @exons << exon
  end

  def addCoverage(position, length)
    #puts "#{name} - adding coverage at position #{position-@start} for #{length} bases"
    s = position-@start
    (s..(s+length)).each do |i|
      if i >= 0 and i < @coverage.length
        @coverage[i]+=1
      end
    end
  end

  def addCoverage2(position, length)
    #puts "#{name} - adding coverage2 at position #{position-@start} for #{length} bases"
    s = position-@start
    (s..(s+length)).each do |i|
      if i >= 0 and i < @coverage2.length
        @coverage2[i]+=1
      end
    end
  end

  def to_s
    "#{name} #{@chromosome} #{@start} #{@stop} #{@strand} #{@exons.size}exons"
  end

  def contains?(chrom, position)
    if chrom == @chromosome
      if position >= @start and position < @stop
        found=false
        exons.each do |ex|
          if ex.contains?(position)
            found=true
          end
        end
        return found
      else
        return false
      end
    else
      return false
    end
  end

end

class Exon
  attr_accessor :start, :stop, :type

  def initialize(start, stop, type)
    @start = start.to_i
    @stop = stop.to_i
    @type = type
  end

  def contains?(position)
    if position >= @start and position < @stop
      return true
    end
  end

  def to_s
    "#{start} #{stop}"
  end
end

sam = opts.sam
saml = sam.split(/\//)
samname = saml[saml.length-1]
saml = samname.split(/\./)
samname = saml[0]

sam2 = opts.sam2
saml2 = sam2.split(/\//)
samname2 = saml2[saml2.length-1]
saml2 = samname2.split(/\./)
samname2 = saml2[0]

proteins=Hash.new

deinfo=Hash.new
Splice = Struct.new(:name, :gdc, :s35, :cell, :log2)
print "Reading At_DE.txt file..."
File.open("At_DE.txt", "r").each_line do |line| # TODO put this back to the full file
  line.chomp!
  cols = line.split(/\t/)
  if cols[0] =~ /(\S+)\.\d+/
    gene = $1
  else
    gene = cols[0]
  end
  if deinfo.has_key?(gene)
    deinfo[gene] << Splice.new(cols[0], cols[9], cols[10], cols[12], cols[13])
  else
    deinfo[gene] = []
    deinfo[gene] << Splice.new(cols[0], cols[9], cols[10], cols[12], cols[13])
  end
end
puts "Done"

count=1
out = `wc -l #{opts.gff}`
lines = out.split(/\s+/)[0].to_i

print "Reading gff file..."
File.open("#{opts.gff}", "r").each_line do |line|
  if count%10000 == 0
    print "#{(100*count/lines).round(0)}%.."
  end
  line.chomp!
  cols = line.split(/\t/)
  deinfo.each_key do |genename|
    if cols[2]=~ /mRNA/
      if cols[8] =~ /ID=(#{genename}\S+?);/
        # key is in form AT1G15530.1
        proteins[$1] = Protein.new(cols[0], cols[3], cols[4], cols[6], $1)
      end
    elsif cols[2]=~ /UTR/
      if cols[8] =~ /=(#{genename}\S+)/
        exon = Exon.new(cols[3], cols[4], "utr")
        if proteins.has_key?($1)
          proteins[$1].addExon(exon)
          #puts "adding an exon to protein #{$1}" 
        else
          puts "Couldn't find #{$1} in protein hash"
        end
      end
    elsif cols[2] =~/exon/
      if cols[8] =~ /=(#{genename}\S+)/
        exon = Exon.new(cols[3], cols[4], "exon")
        if proteins.has_key?($1)
          proteins[$1].addExon(exon)
          #puts "adding an exon to protein #{$1}" 
        else
          puts "Couldn't find #{$1} in protein hash"
        end
      end
    end
  end
  count+=1
end
puts "Done"

count=1
j=0
print "Reading sam file.."
File.open("#{opts.sam}", "r").each_line do |line|
  if count%1000000 == 0
    print "#{j}.."
    j+=1
  end
  line.chomp!
  cols = line.split(/\t/)
  chrom = cols[2].to_i
  pos = cols[3].to_i
  proteins.each_key do |k|
    if proteins[k].chromosome == chrom
      if proteins[k].contains?(chrom, pos)
        proteins[k].addCoverage(pos, cols[9].length)
      end
    end
  end
  count+=1
end
puts "Done"


count=1
j=0
print "Reading 2nd sam file.."
File.open("#{opts.sam2}", "r").each_line do |line|
  if count%1000000 == 0
    print "#{j}.."
    j+=1
  end
  line.chomp!
  cols = line.split(/\t/)
  chrom = cols[2].to_i
  pos = cols[3].to_i
  proteins.each_key do |k|
    if proteins[k].chromosome == chrom
      if proteins[k].contains?(chrom, pos)
        proteins[k].addCoverage2(pos, cols[9].length)
      end
    end
  end
  count+=1
end
puts "Done"

puts "Drawing svg files..."

########################################################################################################

deinfo.each_key do |genename| # genename is in the form AT4G70410

  yscale = 2
  print "  Drawing #{genename}.."

  #find maxcoverage and scale height
  maxlength=0

  deinfo[genename].each do |struct|    
    if proteins.has_key?(struct.name)
      p = proteins[struct.name]
      p.exons.each do |ex|
        if ex.stop-p.start > maxlength
          maxlength = ex.stop-p.start
        end
      end
    end
  end

  #puts "maxlength for #{genename} is #{maxlength}"

  xoffset = 110 # offset from left edge
  width = maxlength*xscale + xoffset*2 + 50
  yoffset = 36
  exonheight = 50
  height=0
  maxx=0
  ystep=10 # 
  xoffset = 110 # offset from left edge
  gap = 60

  c = deinfo[genename].length
  #puts "Number of splice variants #{c}"
  totalheight = c*(100 + exonheight + gap + 30) + xoffset
  img = Rasem::SVGImage.new(width,totalheight) do

    text 10, 24, "Mapping #{samname} and #{samname2}reads to #{genename}", "font-size"=>18

    deinfo[genename].each_with_index do |struct, index|
      splice = struct.name

      #puts splice

      prot = proteins[splice]

      maxcov = 0
      prot.coverage.each do |i|
        if i > maxcov
          maxcov = i
        end
      end
      prot.coverage2.each do |i|
        if i > maxcov
          maxcov = i
        end
      end
      if maxcov * yscale > 100
        yscale = 100.0 / maxcov
        #puts "Setting yscale to #{yscale}"
      end
      while yscale * ystep < 20
        ystep+=5
      end
      
      height = 30 + yscale * maxcov

      list = prot.coverage
      list2 = prot.coverage2
      # coverage graph
      list.each_index do |i|
        if i>0
          x1 = xoffset + (i-1)*xscale
          x2 = xoffset + i*xscale
          y1 = yoffset + height - list[i-1] * yscale
          y2 = yoffset + height - list[i] * yscale
          if list[i]>0 and list[i-1]>0
            line x1, y1, x2, y2, :stroke=>"red"
          else
            line x1, y1, x2, y2, :stroke=>"gray"
          end
        end
      end
      list2.each_index do |i|
        if i>0
          x1 = xoffset + (i-1)*xscale
          x2 = xoffset + i*xscale
          y1 = yoffset + height - list2[i-1] * yscale
          y2 = yoffset + height - list2[i] * yscale
          if list2[i]>0 and list2[i-1]>0
            line x1, y1, x2, y2, :stroke=>"blue"
          else
            #line x1, y1, x2, y2, :stroke=>"gray"
          end
        end
      end

      #exon blocks
      x1 = xoffset
      x2 = xoffset + xscale * (prot.stop - prot.start)
      y = yoffset + gap + height - exonheight
      y2 = y1 = yoffset + gap + height - exonheight + (exonheight/2)
      line x1, y1, x2, y2, :stroke=>"black"
      #puts "number of exons = #{prot.exons.length}"
      prot.exons.each do |ex|
        x = xoffset + xscale * (ex.start - prot.start)
        l = xscale * (ex.stop - ex.start)
        if ex.type=="exon"
          rectangle x, y, l, exonheight, :fill=>"black"  # , :stroke=>"black"
        end
      end
      prot.exons.each do |ex|
        x = xoffset + xscale * (ex.start - prot.start)
        l = xscale * (ex.stop - ex.start)
        if ex.type=="utr"
          rectangle x, y, l, exonheight, :fill=>"gray" # , :stroke=>"black"
        end
      end

      #vertical scale
      x2 = x1 = xoffset-10
      y1 = yoffset + height
      y2 = yoffset + height - (maxcov+ystep) * yscale
      line x1,y1,x2,y2, :stroke=>"black"
      (0..maxcov+ystep).step(ystep) do |i|
        y = yoffset + height - i*yscale
        textoffset = 25
        if i>999
          textoffset = 40
        elsif i>99
          textoffset = 35
        elsif i>9
          textoffset = 30
        end
        text xoffset-textoffset, y+4, "#{i}", "font-size"=>10
        line x1, y, x1-8, y, :stroke=>"black"
      end

      #horizontal scale
      x1 = xoffset
      x2 = xoffset + (prot.stop - prot.start)*xscale
      y1 = y2 = yoffset + gap + height + 5
      line x1,y1,x2,y2, :stroke=>"black"
      (0..(prot.stop-prot.start)).step(200) do |i|
        textoffset=3
        if i > 99
          textoffset=10
        end
        text x1+i*xscale-textoffset, y1+20, "#{i}", "font-size"=>10
        line x1+i*xscale, y1, x1+i*xscale, y1+8, :stroke=>"black"
      end

      #splice label
      text 10, yoffset+height+gap-60, "#{prot.strand}", "font-size"=>24
      text 10, yoffset+height+gap-30, "#{splice}", "font-size"=>12

      #cell type
      text 10, yoffset+height+gap-10, "#{deinfo[genename][index].cell}", "font-size"=>18
      text xoffset+maxlength*xscale+10, yoffset+height+gap+4, "GDC.mean: #{deinfo[genename][index].gdc}", "font-size"=>14
      text xoffset+maxlength*xscale+10, yoffset+height+gap+18, "35S.mean: #{deinfo[genename][index].s35}", "font-size"=>14

      # legend
      line xoffset+maxlength*xscale+10, yoffset+height, xoffset+maxlength*xscale+30, yoffset+height, :stroke=>"red"
      line xoffset+maxlength*xscale+10, yoffset+height+14, xoffset+maxlength*xscale+30, yoffset+height+14, :stroke=>"blue"

      text xoffset+maxlength*xscale+33, yoffset+height+4, "#{samname}", "font-size"=>12
      text xoffset+maxlength*xscale+33, yoffset+height+18, "#{samname2}", "font-size"=>12


      yoffset += maxcov*yscale + exonheight + gap + 30
    end
  end

  File.open("#{genename}-#{samname}-#{samname2}.svg", 'w') { |file| file.write(img.output) }

  puts "Done"
end
