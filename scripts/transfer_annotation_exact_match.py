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
    has_hap = False
    for contig_name, contig in iterate_contigs(fasta_in):
        ctg_suffix = contig_name.split("_")[-1]
        has_hap = has_hap or "hap" in contig_name
        contig_name = "_".join(contig_name.split("_")[:-1]).replace("hap", "")
        contigs[contig_name] = contig
    return contigs, ctg_suffix, has_hap

def nth_index(haystack, needle, n):
    start = haystack.find(needle)
    # print("n_th_index 1", file=sys.stderr)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+1)
        n -= 1
    # print("n_th_index", file=sys.stderr)
    return start

def get_n(haystack, needle, start):
    n = 0
    # print("get_n 1", file=sys.stderr)
    last_start = 0
    while haystack.find(needle, last_start) != start - 1:
        # print("last_start", last_start, file=sys.stderr)
        n += 1
        last_start = haystack.find(needle, last_start) + 1
    # print("get_n", file=sys.stderr)

    return n

def transfer_annotation_exact_match(old_genome, new_genome, annotation, annos_not_found):
    # print("loading", file=sys.stderr)
    old_contigs, _, _ = load_contigs(old_genome)
    new_contigs, new_suffix, new_has_hap = load_contigs(new_genome)
    def fix_ctg_name(ctg):
        if new_has_hap:
            ctg = ctg.replace("_A", "_hapA").replace("_B", "_hapB")
        return ctg + "_" + new_suffix
    # print(old_contigs.keys(), file=sys.stderr)
    # print(new_contigs.keys(), file=sys.stderr)
    # print("loaded", file=sys.stderr)
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
                ctg = "_".join(ctg.split("_")[:-1]).replace("hap", "")
                # print(ctg, a, b, start, end, sep="\t", file=sys.stderr)
                if ctg in old_contigs and ctg in new_contigs:
                    seq = old_contigs[ctg][int(start) - 1:int(end)]
                    was_transferred = False
                    for cut_ends in range(0, max(1, len(seq) - 100) , 100):
                        if cut_ends == 0:
                            seq_curr = seq
                        else:
                            seq_curr = seq[cut_ends:]
                        # print("extracted", file=sys.stderr)
                        if seq_curr in new_contigs[ctg]:
                            # print("found", file=sys.stderr)
                            if new_contigs[ctg].count(seq_curr) != old_contigs[ctg].count(seq_curr):
                                # print("different counts", file=sys.stderr)
                                while len(extra) <= 3:
                                    extra.append("")
                                if extra[3] != "":
                                    extra[3] += ";"
                                extra[3] += "cannot_transfer_reason=different_number_of_occurrences," + str(old_contigs[ctg].count(seq_curr)) + "," + str(new_contigs[ctg].count(seq_curr)) + ",len(seq)=" + str(len(seq_curr))
                                out_file.write("\t".join([fix_ctg_name(ctg), a, b, start, end, *extra]) + "\n")
                                was_transferred = True
                            else:
                                # print("same counts", file=sys.stderr)
                                new_start = nth_index(new_contigs[ctg], seq_curr, 
                                                      get_n(old_contigs[ctg], seq_curr, int(start) + cut_ends))
                                print(fix_ctg_name(ctg), a, b, new_start + 1 - cut_ends, 
                                                    new_start + len(seq_curr) + cut_ends, *extra, sep="\t")
                                was_transferred = True
                                if cut_ends != 0:
                                    print("##", ctg, a, b, start, end, "needed to be cut down to be transferred",
                                          cut_ends, sep="\t", file=sys.stderr)
                            break
                        # print("##", ctg, a, b, start, end, "could not be transferred with cut_ends=",
                        #         cut_ends, sep="\t", file=sys.stderr)
                    if not was_transferred:
                        # print("not found", file=sys.stderr)
                        while len(extra) <= 3:
                            extra.append("")
                        if extra[3] != "":
                            extra[3] += ";"
                        extra[3] += "cannot_transfer_reason=no_match"
                        out_file.write("\t".join([fix_ctg_name(ctg), a, b, start, end, *extra]) + "\n")
                else:
                    while len(extra) <= 3:
                        extra.append("")
                    if extra[3] != "":
                        extra[3] += ";"
                    extra[3] += "cannot_transfer_reason=contig_not_found"
                    out_file.write("\t".join([fix_ctg_name(ctg), a, b, start, end, *extra]) + "\n")



if __name__ == '__main__':
    transfer_annotation_exact_match(*sys.argv[1:])