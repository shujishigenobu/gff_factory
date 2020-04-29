require './gff_factory/lib/NCBI_refseq_genome_gff3'

gff = ARGV[0]

gdb = GffGeneDb.new(gff)
dumpfile = File.basename(ARGV[0], ".gff3") + ".GffGeneDb.dump"
gdb.dump(dumpfile)

