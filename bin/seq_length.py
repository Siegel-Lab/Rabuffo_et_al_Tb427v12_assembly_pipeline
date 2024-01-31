#!/usr/bin/python3
from Bio import SeqIO
import sys
cmdargs = str(sys.argv)
for seq_record in SeqIO.parse(str(sys.argv[1]), "fasta"):
 output_line = '%s\t%i' % \
(seq_record.id, len(seq_record))
 print(output_line)
 
 # To run: chmod +x seq_length.py
 # python seq_length.py input_file.fasta
 
 # Raul's note: Prints in a table format from the input_file.fasta the name and the length
 # of each fasta sequence in the file.