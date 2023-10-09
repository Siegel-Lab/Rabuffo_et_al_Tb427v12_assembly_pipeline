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




def fixup_number_of_n(fasta_in, num_in, num_out):
    str_in = ""
    for _ in range(int(num_in)):
        str_in += "N"
    str_out = ""
    for _ in range(int(num_out)):
        str_out += "N"
    for contig_name, contig in iterate_contigs(fasta_in):
        contig = contig.upper().replace(str_in, str_out)
        print(">" + contig_name)
        print(re.sub("(.{60})", "\\1\n", contig, 0, re.DOTALL))



if __name__ == "__main__":
    fixup_number_of_n(*sys.argv[1:])

