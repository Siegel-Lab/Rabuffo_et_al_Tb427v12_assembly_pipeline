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
            fields = line[:-1].split("\t")
            if fields[0] not in ret:
                ret[fields[0]] = []
            ret[fields[0]].append(fields)
    return ret

def annotate_closed_gaps(old_genome, new_genome, annotation_file, transfer_failed_file, overhang=1000, gap_size=1000):
    old_contigs = load_contigs(old_genome)
    new_contigs = load_contigs(new_genome)
    annotations = load_gff(annotation_file)

    with open(transfer_failed_file, "w") as out_file:
        print("##gff-version 3")
        print("##gff-version 3", file=out_file)
        for contig_name, new_scaffolds in sorted(new_contigs.items()):
            print("##sequence-region", contig_name, 1, 
                sum(len(new_scaffold) + gap_size for new_scaffold in new_scaffolds) - gap_size, sep="\t")
            print("##sequence-region", contig_name, 1, 
                sum(len(new_scaffold) + gap_size for new_scaffold in new_scaffolds) - gap_size, sep="\t", file=out_file)

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
                        # elif old_scaffold_rev_comp in new_scaffold:
                        #     curr_assignment.append(idx)
                    if not len(curr_assignment) == 1:
                        print("WARNING: could not find assignment for:", contig_name, "scaffold", jdx, curr_assignment,
                            file=sys.stderr)
                        # if len(scaffold_assignment) > 0:
                        #     print("last assignment:", scaffold_assignment[-1], file=sys.stderr)
                        # print(contig_name)
                        # print(curr_assignment)
                        # print(jdx)
                        # print(len(old_scaffolds), file=sys.stderr)
                        # print(len(new_scaffolds), file=sys.stderr)
                        # print(len(old_scaffold), file=sys.stderr)
                        # print(old_scaffold, file=sys.stderr)
                        # print(old_scaffold_rev_comp)
                        # assert False
                        scaffold_assignment.append(None)
                    else:
                        scaffold_assignment.append(curr_assignment[0])

                assert all(not None in [a, b] or a <= b for a, b in zip(scaffold_assignment[:-1], scaffold_assignment[1:]))

                for old_sc, new_sc in enumerate(scaffold_assignment):
                    old_sc_start = sum(len(s) + gap_size for s in old_scaffolds[:old_sc]) + overhang
                    old_sc_end = old_sc_start + len(old_scaffolds[old_sc]) - overhang*2
                    
                    cropped_scaffold_a = old_scaffolds[old_sc][overhang:-overhang]
                    assert new_scaffolds[new_sc].count(cropped_scaffold_a) == 1
                    idx_new = new_scaffolds[new_sc].index(cropped_scaffold_a)
                    # idx_new = 0
                    # if cropped_scaffold_a in new_scaffolds[new_sc]:
                    #     idx_new = new_scaffolds[new_sc].index(cropped_scaffold_a)
                    # else:
                    #     assert rev_comp(cropped_scaffold_a) in new_scaffolds[new_sc]
                    #     idx_new = new_scaffolds[new_sc].index(rev_comp(cropped_scaffold_a))

                    new_sc_start = sum(len(s) + gap_size for s in new_scaffolds[:new_sc]) + idx_new
                    new_sc_end = new_sc_start + len(cropped_scaffold_a)

                    remaining_annos = []
                    for annotation in annotations[contig_name]:
                        start = int(annotation[3])
                        end = int(annotation[4])

                        if start > old_sc_start and end <= old_sc_end:
                            new_start = start + new_sc_start - old_sc_start
                            new_end = end + new_sc_start - old_sc_start

                            seq_a = cropped_scaffold_a[start - old_sc_start - 1:end - old_sc_start]
                            seq_b = new_scaffolds[new_sc][new_start - new_sc_start - 1 + idx_new:new_end - new_sc_start + idx_new]
                            if seq_a != seq_b:
                                print("WARNING: annotation sequence does not match", contig_name, annotation, 
                                      old_sc, new_sc, old_sc_start, old_sc_end, new_sc_start, new_sc_end, "\n", 
                                      seq_a, "\n",
                                      seq_b, "\n", new_scaffolds[new_sc].index(seq_a), start, end, new_start, new_end, 
                                      idx_new, file=sys.stderr)
                                assert False

                            print(*annotation[:3], new_start, new_end, *annotation[5:], sep="\t")
                        else:
                            remaining_annos.append(annotation)
                    annotations[contig_name] = remaining_annos
            elif contig_name in new_contigs and len(old_scaffolds) == 1:
                for annotation in annotations[contig_name]:
                    print(*annotation, sep="\t")


        for annotation in annotations.values():
            for a in annotation:
                print(*a, sep="\t", file=out_file)




if __name__ == '__main__':
    annotate_closed_gaps(*sys.argv[1:])