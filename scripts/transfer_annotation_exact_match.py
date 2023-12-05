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

def nth_index(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+1)
        n -= 1
    return start

def get_n(haystack, needle, start):
    n = 0
    last_start = 0
    while haystack.find(needle, start) != start:
        n += 1
        start = haystack.find(needle, start) + len(needle)

    return n

def transfer_annotation_exact_match(old_genome, new_genome, annotation, annos_not_found):
    old_contigs = load_contigs(old_genome)
    new_contigs = load_contigs(new_genome)
    print("##gff-version 3")
    for name, seq in new_contigs.items():
        print("##sequence-region", name, "1", len(seq), sep="\t")
    with open(annos_not_found, "w") as out_file:
        out_file.write("##gff-version 3\n")
        for name, seq in old_contigs.items():
            out_file.write(" ".join(["##sequence-region", name, "1", str(len(seq))]) + "\n")
        with fileinput.input(annotation) as in_file:
            for line in in_file:
                if line[0] == "#":
                    continue
                ctg, a, b, start, end, *extra = line.strip().split("\t")
                if ctg in old_contigs and ctg in new_contigs:
                    seq = old_contigs[ctg][int(start):int(end)]
                    if seq in new_contigs[ctg]:
                        if new_contigs[ctg].count(seq) != old_contigs[ctg].count(seq):
                            while len(extra) <= 3:
                                extra.append("")
                            if extra[3] != "":
                                extra[3] += ";"
                            extra[3] += "cannot_transfer_reason=different_number_of_occurrences"
                            out_file.write("\t".join([ctg, a, b, start, end, *extra]) + "\n")
                        else:
                            new_start = nth_index(new_contigs[ctg], seq, get_n(old_contigs[ctg], seq, int(start)))
                            print(ctg, a, b, new_start + 1, new_start + len(seq) + 1, *extra, sep="\t")
                    else:
                        while len(extra) <= 3:
                            extra.append("")
                        if extra[3] != "":
                            extra[3] += ";"
                        extra[3] += "cannot_transfer_reason=no_match"
                        out_file.write("\t".join([ctg, a, b, start, end, *extra]) + "\n")
                else:
                    while len(extra) <= 3:
                        extra.append("")
                    if extra[3] != "":
                        extra[3] += ";"
                    extra[3] += "cannot_transfer_reason=contig_not_found"
                    out_file.write("\t".join([ctg, a, b, start, end, *extra]) + "\n")



if __name__ == '__main__':
    transfer_annotation_exact_match(*sys.argv[1:])