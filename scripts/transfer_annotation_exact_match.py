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
    print("##gff-version 3")
    for name, seq in new_contigs.items():
        print("##sequence-region", name, "1", len(seq), sep="\t")
    with open(annos_not_found, "w") as out_file:
        with fileinput.input(annotation) as in_file:
            for line in in_file:
                if line[0] == "#":
                    continue
                ctg, a, b, start, end, *extra = line.strip().split("\t")
                if ctg in old_contigs:
                    seq = old_contigs[ctg][int(start):int(end)]
                    if seq in new_contigs[ctg]:
                        new_start = new_contigs[ctg].index(seq)
                        if seq in new_contigs[ctg][new_start + len(seq):]:
                            out_file.write("multiple_matches:\t" + line)
                        else:
                            print(ctg, a, b, new_start + 1, new_start + len(seq) + 1, *extra, sep="\t")
                    else:
                        out_file.write("no_match:\t" + line)
                else:
                    out_file.write("contig_not_found:\t" + line)



if __name__ == '__main__':
    transfer_annotation_exact_match(*sys.argv[1:])