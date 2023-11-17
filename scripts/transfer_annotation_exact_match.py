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

def transfer_annotation_exact_match(old_genome, new_genome, annotation, annos_not_found):
    old_contigs = load_contigs(old_genome)
    new_contigs = load_contigs(new_genome)
    with open(annos_not_found, "w") as out_file:
        with fileinput.input(annotation) as in_file:
            for line in in_file:
                if line[0] == "#":
                    print(line, end="")
                    continue
                ctg, a, b, start, end, *extra = line.strip().split()
                seq = old_contigs[ctg][int(start):int(end)]
                if ctg in new_contigs and seq in new_contigs[ctg]:
                    new_start = new_contigs[ctg].index(seq)
                    print(ctg, a, b, new_start, new_start + len(seq), *extra, sep="\t")
                else:
                    out_file.write(line)



if __name__ == '__main__':
    transfer_annotation_exact_match(*sys.argv[1:])