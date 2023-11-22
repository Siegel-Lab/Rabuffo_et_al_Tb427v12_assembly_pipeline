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
        scaffolds = contig.split("N" * 1000)
        contigs[contig_name] = scaffolds
    return contigs

def rev_comp(seq):
    return "".join({"A": "T", "T": "A", "C": "G", "G": "C"}[s] for s in reversed(seq))

def load_gff(filename):
    ret = {}
    with fileinput.input(filename) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            contig, _, _, start, end, _, _, _, extra = line[:-1].split("\t")
            extra = extra.split(";")
            idx = None
            for e in extra:
                if "ID=" in e:
                    idx = e.split("=")[1]
            start = int(start)
            end = int(end)
            if contig not in ret:
                ret[contig] = []
            ret[contig].append([start, end, idx])
    for l in ret.values():
        l.sort()
    return ret

def annotate_closed_gaps(old_genome, new_genome, old_gaps, overhang=1000, gap_size=1000):
    old_contigs = load_contigs(old_genome)
    new_contigs = load_contigs(new_genome)
    gaps = load_gff(old_gaps)

    print("##gff-version 3")
    for contig_name, new_scaffolds in sorted(new_contigs.items()):
        print("##sequence-region", contig_name, 1, 
              sum(len(new_scaffold) + gap_size for new_scaffold in new_scaffolds) - gap_size, sep="\t")

    for contig_name, old_scaffolds in old_contigs.items():
        if contig_name in new_contigs and len(old_scaffolds) > 1:
            new_scaffolds = new_contigs[contig_name]
            scaffold_assignment = []
            for jdx, old_scaffold in enumerate(old_scaffolds):
                cropped_scaffold = old_scaffold[overhang:-overhang]
                old_scaffold_rev_comp = rev_comp(cropped_scaffold)
                curr_assignment = []
                for idx, new_scaffold in enumerate(new_scaffolds):
                    if cropped_scaffold in new_scaffold:
                        curr_assignment.append(idx)
                    if old_scaffold_rev_comp in new_scaffold:
                        curr_assignment.append(idx)
                if not len(curr_assignment) == 1:
                    print("ERROR: could not find assignment for:", contig_name, "scaffold", jdx, len(curr_assignment),
                          file=sys.stderr)
                    # print(contig_name)
                    # print(curr_assignment)
                    # print(jdx)
                    # print(len(old_scaffolds))
                    # print(len(new_scaffolds))
                    # print(len(old_scaffold))
                    # print(old_scaffold)
                    # print(old_scaffold_rev_comp)
                    assert False
                else:
                    scaffold_assignment.append(curr_assignment[0])
            assert all(None in [a, b] or a <= b for a, b in zip(scaffold_assignment[:-1], scaffold_assignment[1:]))

            for (old_sc_a, new_sc_a), (old_sc_b, new_sc_b) in zip(enumerate(scaffold_assignment[:-1]), 
                                                                  enumerate(scaffold_assignment[1:])):
                idx = "ID=?"
                old_prev_scaf_len = sum(len(s) + gap_size for s in old_scaffolds[:old_sc_a])
                old_idx_a = old_prev_scaf_len + len(old_scaffolds[old_sc_a])
                old_idx_b = old_idx_a + gap_size
                for start, end, idy in gaps[contig_name]:
                    if start == old_idx_a + 1 and end == old_idx_b:
                        idx = "ID="+str(idy)
                        break
                if new_sc_a == new_sc_b:
                    cropped_scaffold_a = old_scaffolds[old_sc_a][overhang:-overhang]
                    cropped_scaffold_b = old_scaffolds[old_sc_b + 1][overhang:-overhang]
                    prev_scaf_len = sum(len(s) + gap_size for s in new_scaffolds[:new_sc_a])
                    idx_a = new_scaffolds[new_sc_a].index(cropped_scaffold_a) + len(cropped_scaffold_a) + prev_scaf_len
                    idx_b = new_scaffolds[new_sc_a].index(cropped_scaffold_b) + prev_scaf_len
                    print(contig_name, ".", "filledgap", idx_a + 1, idx_b, ".", ".", ".", idx, sep="\t")
                else:
                    prev_scaf_len = sum(len(s) + gap_size for s in new_scaffolds[:new_sc_a])
                    idx_a = prev_scaf_len + len(new_scaffolds[new_sc_a])
                    idx_b = idx_a + gap_size
                    print(contig_name, ".", "gap", idx_a + 1, idx_b, ".", ".", ".", 
                          "estimated_length=1000;gap_type=within scaffold;" + idx, sep="\t")


if __name__ == '__main__':
    annotate_closed_gaps(*sys.argv[1:])