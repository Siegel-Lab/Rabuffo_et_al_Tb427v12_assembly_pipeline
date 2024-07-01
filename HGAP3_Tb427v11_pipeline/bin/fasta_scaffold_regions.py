#!/usr/bin/python
from Bio import SeqIO
import sys
cmdargs = str(sys.argv)



with open(sys.argv[1]) as in_f, open(sys.argv[2],'w') as out_f:
	out_f.write('##gff-version 3\n')
	for record in SeqIO.parse(in_f, 'fasta'):
		out_f.write('{}\tkDNA\tgene\t1\t{}\t.\t+\t.\tID=g{};\n'.format(record.id,len(record),record.id))
		out_f.write('{}\tkDNA\tmRNA\t1\t{}\t.\t+\t.\tID=g{}.1;Parent=g{};\n'.format(record.id,len(record),record.id,record.id))
		out_f.write('{}\tkDNA\texon\t1\t{}\t.\t+\t.\tID=exon_g{};Parent=g{}.1;\n'.format(record.id,len(record),record.id,record.id))
		

