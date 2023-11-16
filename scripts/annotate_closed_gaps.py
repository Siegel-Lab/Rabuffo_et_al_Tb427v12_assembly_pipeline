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

def annotate_closed_gaps(old_genome, new_genome):
    old_contigs = load_contigs(old_genome)
    new_contigs = load_contigs(new_genome)

    for contig_name, old_scaffolds in old_contigs:
        if contig_name in new_contigs and len(old_scaffolds) > 1:
            new_scaffolds = new_contigs[contig_name]
            scaffold_assignment = []
            for old_scaffold in old_scaffolds:
                curr_assignment = []
                for idx, new_scaffold in enumerate(new_scaffolds):
                    if old_scaffold in new_scaffold:
                        curr_assignment.append(idx)
                assert len(curr_assignment) == 1
                scaffold_assignment.append(curr_assignment[0])
            assert all(a <= b for a, b in zip(scaffold_assignment[:-1], scaffold_assignment[1:]))

            for (old_sc_a, new_sc_a), (old_sc_b, new_sc_b) in zip(enumerate(scaffold_assignment[:-1]), 
                                                                  enumerate(scaffold_assignment[1:])):
                if new_sc_a == new_sc_b:
                    prev_scaf_len = sum(len(s) for s in new_scaffolds[:new_sc_a])
                    idx_a = new_scaffolds[new_sc_a].index(old_scaffolds[old_sc_a]) + len(old_scaffolds[old_sc_a]) + prev_scaf_len
                    idx_b = new_scaffolds[new_sc_a].index(old_scaffolds[old_sc_b]) + prev_scaf_len
                    print(contig_name, ".", "filledgap", idx_a, idx_b, ".", ".", ".", sep="\t")


if __name__ == '__main__':
    annotate_closed_gaps(*sys.argv[1:])