import fileinput
import sys
from paste_old_sequence_into_remaining_gaps import load_masked


def distinguish_masked_from_unmasked(gaps_before_closing, gaps_after_masking_undone, masked_in):
    masked = load_masked(masked_in)
    old_gff = {}
    with fileinput.input(gaps_before_closing) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            _, _, _, start, end, _, _, _, extra = line.split("\t")
            for e in extra.split(";"):
                if "ID=" in e:
                    idx = int(e.split("=")[1])
                    old_gff[idx] = (int(start), int(end))
                    break

    with fileinput.input(gaps_after_masking_undone) as in_file:
        for line in in_file:
            if line[0] == "#":
                print(line[:-1])
                continue
            contig, x1, what, start, end, x2, x3, x4, extra = line[:-1].split("\t")
            idx = None
            for e in extra.split(";"):
                if "ID=" in e:
                    idx = int(e.split("=")[1])
            if what == "filledgap":
                old_start, old_end = old_gff[idx]
                key = contig + " " + str(old_start) + " " + str(old_end)
                if key in masked:
                    es = extra.split(";")
                    for idx, e in enumerate(es):
                        if "previous_size" in e:
                            es[idx] = "previous_size=" + str(len(masked[key]))
                    extra = ";".join(es)
                print(contig, x1, "filledmasked" if key in masked else "filledgap", start, end, x2, x3, x4, extra, 
                      sep="\t")
            else:
                print(line[:-1])

if __name__ == "__main__":
    distinguish_masked_from_unmasked(*sys.argv[1:])