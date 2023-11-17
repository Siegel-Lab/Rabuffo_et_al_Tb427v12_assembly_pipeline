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

def load_contigs(fasta_in):
    contigs = {}
    for contig_name, contig in iterate_contigs(fasta_in):
        contigs[contig_name] = contig
    return contigs

def extract_region_from_annotation(genome, annotation):
    contigs = load_contigs(genome)
    with fileinput.input(annotation) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            ctg, a, b, start, end, *extra = line.strip().split()
            seq = contigs[ctg][int(start):int(end)]
            print(ctg, start, end, seq, sep="\t")

if __name__ == "__main__":
    extract_region_from_annotation(*sys.argv[1:])
