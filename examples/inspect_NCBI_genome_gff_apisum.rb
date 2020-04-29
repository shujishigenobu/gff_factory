require './gff_factory/lib/NCBI_refseq_genome_gff3'

dumpfile = ARGV[0]

gdb = GffGeneDb.load(dumpfile)

puts "#" + %w{gene mRNA mRNA_length CDS_length protein product AphidBase}.join("\t")

gdb.genes.keys.each do |k|
  g = gdb.genes[k]
  if g.data['gene_biotype'] == 'protein_coding'
    #    p [k, g.children.map{|c| c.id}]
    g.children.each do |tr|
      puts [g.data['Name'], 
            tr.data['Name'], 
            tr.length, 
            tr.cds_length, 
            tr.protein['protein_id'], 
            tr.protein['product'], 
            g.dbxref["APHIDBASE"]
           ].join("\t")
    end
  end
end
