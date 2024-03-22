import fileinput
import sys

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

def identify_telomeric_repeats(genome_in):
        print("contig", "start_TTAGGG", "start_CCCTAA", "end_TTAGGG", "end_CCCTAA", sep="\t")
        for contig_name, contig in iterate_contigs(genome_in):
            print(contig_name, 
                  "yes" if "TTAGGG" in contig[:100] else "no", 
                  "yes" if "CCCTAA" in contig[:100] else "no", 
                  "yes" if "TTAGGG" in contig[-100:] else "no", 
                  "yes" if "CCCTAA" in contig[-100:] else "no", 
                  sep="\t")


if __name__ == "__main__":
    identify_telomeric_repeats(*sys.argv[1:])