import fileinput
import sys
import re

def iterate_contigs(fasta_in):
    contig = ""
    contig_name = None
    with fileinput.input(fasta_in) as in_file:
        for line in in_file:
            if line[0] == ">":
                if not contig_name is None:
                    yield contig_name, contig
                contig = ""
                contig_name = line[1:-1]
                continue
            contig += line[:-1]
        if not contig_name is None:
            yield contig_name, contig




def join_contigs_with_n(fasta_in, num_n):
    str_n = ""
    for _ in range(int(num_n)):
        str_n += "N"
    print(">whole_genome")
    genome = str_n.join([contig for _, contig in iterate_contigs(fasta_in)])
    print(re.sub("(.{60})", "\\1\n", genome, 0, re.DOTALL))



if __name__ == "__main__":
    join_contigs_with_n(*sys.argv[1:])

