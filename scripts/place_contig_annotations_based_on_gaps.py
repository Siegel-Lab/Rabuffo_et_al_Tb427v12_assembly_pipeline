import fileinput
import sys

def place_contig_annotations_based_on_gaps(old_gaps, new_gaps_in, contigs_in, tb427Ver, overlap_region=1000):
    old_gaps_by_contig = {}
    with fileinput.input(old_gaps) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            contig, _, _, start, end, _, _, _, *extra = line.split("\t")
            start=int(start)
            end=int(end)
            if contig not in old_gaps_by_contig:
                old_gaps_by_contig[contig] = []
            old_gaps_by_contig[contig].append([start, end])
    for l in old_gaps_by_contig.values():
        l.sort()
    # print(old_gaps_by_contig.keys(), file=sys.stderr)


    contig_lengths = {}
    new_gap_by_contig = {}
    with fileinput.input(new_gaps_in) as in_file:
        for line in in_file:
            if line[0] == "#":
                if line.startswith("##sequence-region"):
                    _, contig, _, end = line.split()
                    contig_lengths[contig] = int(end)
                continue
            contig, _, gap_type, start, end, _, _, _, *extra = line.split("\t")
            start=int(start)
            end=int(end)
            if contig not in new_gap_by_contig:
                new_gap_by_contig[contig] = []
            new_gap_by_contig[contig].append([start, end, gap_type])
    for contig, l in new_gap_by_contig.items():
        if len(l) != len(old_gaps_by_contig[contig]):
            print(contig, len(l), len(old_gaps_by_contig[contig]), file=sys.stderr)
            assert False
        l.sort()
    # print(new_gap_by_id.keys(), file=sys.stderr)


    contigs = {}
    with fileinput.input(contigs_in) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            contig, _, _, start, end, _, _, _, extra = line.split("\t")
            start=int(start)
            end=int(end)
            name = None
            for field in extra.split(";"):
                if field.startswith("Name="):
                    name = field[5:]
                    break
            name = name[:name.index("_CTG")]
            if name in contigs:
                start = min(start, contigs[name][0])
                end = max(end, contigs[name][1])
                if contigs[name][2] != contig:
                    print("WARNING: contig ", name, "is annotated to both alleles of the fully phased reference. Picking allele A.",
                        file=sys.stderr)
                    continue
            contigs[name] = [start, end, contig]

    for name, (start, end, contig) in contigs.items():
        start_gap = None
        end_gap = None
        if "unitig" in name:
            print(contig, ".", "unitig", 1, end, ".", "+", ".", "Name=" + name, sep="\t")
        elif contig in old_gaps_by_contig and contig in new_gap_by_contig:
            for idx, gap in enumerate(old_gaps_by_contig[contig]):
                # if "_core_" in name or "_3A_" in name or "_3B_" in name:
                if gap[1] + 1 == start:
                    assert start_gap is None
                    start_gap = idx

                # if "_core_" in name or "_5A_" in name or "_5B_" in name:
                if end + 1 == gap[0]:
                    assert end_gap is None
                    end_gap = idx

            start_gap_type = "-/-"
            if start_gap is not None:
                start_gap_type = new_gap_by_contig[contig][start_gap][2]
                if new_gap_by_contig[contig][start_gap][2] != "gap":
                    start = new_gap_by_contig[contig][start_gap][0] + overlap_region
                else:
                    start = new_gap_by_contig[contig][start_gap][0]
            else:
                start = 1

            end_gap_type = "-/-"
            if end_gap is not None:
                end_gap_type = new_gap_by_contig[contig][end_gap][2]
                if new_gap_by_contig[contig][end_gap][2] != "gap":
                    end = new_gap_by_contig[contig][end_gap][0] + overlap_region - 1
                else:
                    end = new_gap_by_contig[contig][end_gap][0] - 1
            else:
                end = contig_lengths[contig]

            # exception for BES15
            if name == "BES15_Tb427v11":
                start = 2309280 - 91008
            if name == "Chr3_3A_Tb427v11":
                end = 2309280 - 91008

            if end <= start:
                print("WARNING: contig ", name, "has negative length", file=sys.stderr)
            print(contig, ".", "contig_core" if "_core_" in name else "contig_subt", start, end, ".", 
                  "+" if "_core_" + tb427Ver + "_A" in name or "_3A_" in name or "_5A_" in name else "-", 
                  ".", "Name=" + name + ";Position_Adjusted_by_Gaps=True;start_gap_id=" + str(start_gap) + ";end_gap_id=" + str(end_gap) + ";start_gap_type=" + start_gap_type + ";end_gap_type=" + end_gap_type , sep="\t")
        else:
            print(contig, ".", "contig_core" if "_core_" in name else "contig_subt", 1, contig_lengths[contig], ".", 
                  "+", ".", "Name=" + name + ";Position_Adjusted_by_Gaps=False", sep="\t")




if __name__ == "__main__":
    place_contig_annotations_based_on_gaps(*sys.argv[1:])