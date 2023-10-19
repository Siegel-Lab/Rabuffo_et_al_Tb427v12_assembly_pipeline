import fileinput
import sys

def close_annotation_in_gff(annotation, remaining_gaps, fixed_gaps=None):
    gaps = {}
    with fileinput.input(annotation) as in_file:
        for line in in_file:
            if line[0] == "#":
                print(line[:-1])
                continue
            contig, source, anno_type, start, end, *extra = line.strip().split()
            if anno_type == "gap":
                if not contig in gaps:
                    gaps[contig] = []
                gaps[contig].append([int(start), "closedgap", source, end, extra])
            else:
                print(line[:-1])

    with fileinput.input(remaining_gaps) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            contig, source, anno_type, start, end, *extra = line.strip().split()

            closest_gap = None
            closest_distance = None

            if contig in gaps:
                for idx, gap in enumerate(gaps[contig]):
                    o_start = gap[0]
                    dist = abs(o_start - int(start))

                    # @todo @fixme i cannot just pick the second closest gap here
                    # I have to check properly, why multiple gaps in the outcome assembly are close to the same gap on the reference assembly
                    if (closest_distance is None or dist < closest_distance) and gap[1] == "closedgap":
                        closest_gap = idx
                        closest_distance = dist
            else:
                print("trying to state that gap on", contig, " is still open. however there is no gap to be found.", file=sys.stderr)

            if not closest_gap is None:
                if gaps[contig][closest_gap][1] == "gap":
                    print("Gap '", contig, closest_gap, "' already designated as unfixed.", file=sys.stderr)
                gaps[contig][closest_gap][1] = "gap"
    if not fixed_gaps is None:
        with fileinput.input(fixed_gaps) as in_file:
            for line in in_file:
                if line[0] == "#":
                    continue
                contig, source, anno_type, start, end, *extra = line.strip().split()

                closest_gap = None
                closest_distance = None

                if contig in gaps:
                    for idx, gap in enumerate(gaps[contig]):
                        o_start = gap[0]
                        dist = abs(o_start - int(start))

                        # @todo @fixme same as above
                        if (closest_distance is None or dist < closest_distance) and gap[1] == "closedgap":
                            closest_gap = idx
                            closest_distance = dist

                    if not closest_gap is None:
                        if gaps[contig][closest_gap][1] != "closedgap":
                            print("Gap", contig, closest_gap, gaps[contig][closest_gap][1], 
                                  "not closed before making it", anno_type, file=sys.stderr)
                        gaps[contig][closest_gap][1] = anno_type
                else:
                    print(line[:-1], "not found", file=sys.stderr)
    for contig, gap_list in gaps.items():
        for start, gap_type, source, end, extra in gap_list:
            print(contig, source, gap_type, start, end, *extra, sep="\t")


if __name__ == "__main__":
    close_annotation_in_gff(*sys.argv[1:])