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

def annotate_gaps(file_in):
    next_gap_id=1
    for contig_name, contig in iterate_contigs(file_in):
        gap_start = None
        gap_end = None
        last_base_was_gap = False
        for idx, contig_base in enumerate(contig):
            is_gap = contig_base.upper() == "N"
            if not last_base_was_gap and is_gap:
                gap_start = idx + 1
            if last_base_was_gap and not is_gap:
                gap_end = idx
                print(contig_name, ".", "gap", gap_start, gap_end, ".", ".", ".", 
                      "estimated_length=1000;gap_type=within scaffold;ID=" + str(next_gap_id), sep="\t")
                next_gap_id += 1

            last_base_was_gap = is_gap
        if last_base_was_gap:
            gap_end = idx
            print(contig_name, ".", "gap", gap_start, gap_end, ".", ".", ".", 
                    "estimated_length=1000;gap_type=within scaffold;ID=" + str(next_gap_id), sep="\t")
            next_gap_id += 1



if __name__ == "__main__":
    annotate_gaps(*sys.argv[1:])