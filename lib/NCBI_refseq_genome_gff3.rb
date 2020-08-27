require 'bio'
require 'pp'

class GffFeature

  def initialize(gffline)
    @seqname
    @source
    @feature
    @start
    @end
    @score
    @strand
    @frame
    @attributes
    @data = Hash.new
    @dbxref = Hash.new
    @synonyms = []
    @children = []
    @gr = nil
    parse_gffline(gffline)
  end
  
  attr_reader :seqname, :source, :feature, :start, :end, :score, :strand, :frame, :attributes
  attr_reader :data, :dbxref, :synonyms
  attr_reader :gr
  attr_accessor :children
  alias_method :gff_line, :gr

  def length
    len = @end - @start + 1
    raise unless len > 0
    return len
  end

  def parse_gffline(gffline)
    @gr = Bio::GFF::GFF3::Record.new(gffline)
    @seqname = @gr.seqname
    @source = @gr.source
    @feature = @gr.feature
    @start = @gr.start
    @end = @gr.end
    @score = @gr.score
    @strand = @gr.strand
    @frame = @gr.frame
    @attributes = @gr.attributes
    parse_attributes
  end

  def parse_attributes
    @gr.attributes.each do |a|
      k, v = a
      if k == "Dbxref"
        m = /^(.+?):(.+)/.match(v)
        k2 = m[1]
        v2 = m[2]
        @dbxref[k2] = v2
        next
      elsif k == "gene_synonym"
        @synonyms << v
        next
      elsif k == "end_range" || k == "start_range" || k == "transl_except"
        ## do nothing
        next
      else

      end

      if @data.has_key?(k)
        raise "duplicated key: #{a.inspect}"
      else
        @data[k] = v
      end
    end

  end

  def id
    data['ID']
  end

  def parent
    data['Parent']
  end
end


class Gene < GffFeature

#  def initialize(gffline)
#  end
  

#  def parse_attributes
#  end

end

class Mrna < GffFeature
  def initialize(gffline)
    super
    @exons = []
    @cds_parts = []
  end
  attr_accessor :exons, :cds_parts

  def length
    exons.map{|e| e.length}.inject{|i, j| i+=j}
  end

  def cds_length
    cds_parts.map{|e| e.length}.inject{|i, j| i+=j}
  end

  def protein
    h = Hash.new
    h['num_cds_exons'] = cds_parts.size
    a = cds_parts.map{|c| c.id}.sort.uniq
    unless a.size == 1
      p cds_parts
      p a
      raise "Not unique protein ID" 
    end
#    cds_parts[0].dbxref
    h['protein_id'] = cds_parts[0].data["protein_id"]
    h['product'] = cds_parts[0].data["product"]
    h
  end

end

class Exon < GffFeature

end

class Cds_part < GffFeature
end

class Rrna < GffFeature
  def initialize(gffline)
    super
    @exons = []
  end
  attr_accessor :exons
end

class Lncrna < GffFeature
  def initialize(gffline)
    super
    @exons = []
  end
  attr_accessor :exons
end

#====

class GffGeneDb

  def self.load(file)
    gdb = Marshal.load(File.open(file).read)
    gdb
  end


  def initialize(gff_file)
    gff = gff_file

    genes = Hash.new
    mrnas = Hash.new
    rrnas = Hash.new
    lncrnas = Hash.new
    transcripts = Hash.new

    File.open(gff).each do |l|
      next if /^#/.match(l)
      puts l
      gr = Bio::GFF::GFF3::Record.new(l)

      if gr.feature == "gene"
        gn = Gene.new(l)
        genes[gn.id] = gn

      elsif gr.feature == "mRNA"
        mrna = Mrna.new(l)
        mrnas[mrna.id] = mrna
        transcripts[mrna.id] = mrna
        parental_gene = genes[mrna.parent]
        parental_gene.children << mrna
      elsif gr.feature == "exon"
        exon = Exon.new(l)
        parental_transcript = transcripts[exon.parent]
        next unless parental_transcript
        parental_transcript.exons << exon
      elsif gr.feature == "CDS"
        cds_part = Cds_part.new(l)
        parental_transcript = transcripts[cds_part.parent]
        parental_transcript.cds_parts << cds_part

      elsif gr.feature == "rRNA"
        rna = Rrna.new(l)
        rrnas[rna.id] = rna
        transcripts[rna.id] = rna
        parental_gene = genes[rna.parent]
        parental_gene.children << rna

      elsif gr.feature == "lnc_RNA"
        rna = Lncrna.new(l)
        lncrnas[rna.id] = rna
        transcripts[rna.id] = rna
        parental_gene = genes[rna.parent]
        parental_gene.children << rna
      end
    end
    @genes = genes
  end

  attr_reader :genes
  
  def dump(file)
    File.open(file, "w"){|o|
      Marshal.dump(self, o)
    }
  end


end

if __FILE__ == $0

  gdb = GffGeneDb.new(ARGV[0])

  gdb.genes.keys.each do |k|
    g = gdb.genes[k]
    if g.data['gene_biotype'] == 'protein_coding'
      #    p [k, g.children.map{|c| c.id}]
      g.children.each do |tr|
        puts [g.data['Name'], tr.data['Name'], tr.length, tr.cds_length, tr.protein['protein_id'], tr.protein['product'], g.dbxref["APHIDBASE"]].join("\t")
      end
    end
  end

  dumpfile = File.basename(ARGV[0], ".gff3") + ".GffGeneDb.dump"
  gdb.dump(dumpfile)
end


