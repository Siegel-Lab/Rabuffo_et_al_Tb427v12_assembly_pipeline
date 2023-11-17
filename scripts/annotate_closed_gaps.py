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

def annotate_closed_gaps(old_genome, new_genome, overhang=1000):
    old_contigs = load_contigs(old_genome)
    new_contigs = load_contigs(new_genome)

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
                    print("ERROR: could not find assignment for:", contig_name, "scaffold", jdx, file=sys.stderr)
                    # print(contig_name)
                    # print(curr_assignment)
                    # print(jdx)
                    # print(len(old_scaffolds))
                    # print(len(new_scaffolds))
                    # print(old_scaffold)
                    # print(old_scaffold_rev_comp)
                    assert False
                    scaffold_assignment.append(None)
                else:
                    scaffold_assignment.append(curr_assignment[0])
            assert all(None in [a, b] or a <= b for a, b in zip(scaffold_assignment[:-1], scaffold_assignment[1:]))

            for (old_sc_a, new_sc_a), (old_sc_b, new_sc_b) in zip(enumerate(scaffold_assignment[:-1]), 
                                                                  enumerate(scaffold_assignment[1:])):
                if None in [new_sc_a, new_sc_b]:
                    continue
                if new_sc_a == new_sc_b:
                    cropped_scaffold_a = old_scaffolds[old_sc_a][overhang:-overhang]
                    cropped_scaffold_b = old_scaffolds[old_sc_b + 1][overhang:-overhang]
                    prev_scaf_len = sum(len(s) for s in new_scaffolds[:new_sc_a])
                    idx_a = new_scaffolds[new_sc_a].index(cropped_scaffold_a) + len(cropped_scaffold_a) + prev_scaf_len
                    idx_b = new_scaffolds[new_sc_a].index(cropped_scaffold_b) + prev_scaf_len
                    print(contig_name, ".", "filledgap", idx_a, idx_b, ".", ".", ".", sep="\t")


if __name__ == '__main__':
    annotate_closed_gaps(*sys.argv[1:])