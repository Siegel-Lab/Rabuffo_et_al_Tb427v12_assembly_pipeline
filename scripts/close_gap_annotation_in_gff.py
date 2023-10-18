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
                gaps[contig].append([int(start), True, False, source, end, extra])
            else:
                print(line[:-1])

    with fileinput.input(remaining_gaps) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            contig, source, anno_type, start, end, *extra = line.strip().split()

            closest_gap = None
            closest_distance = None

            for idx, gap in enumerate(gaps[contig]):
                o_start = gap[0]
                dist = abs(o_start - int(start))

                if closest_distance is None or dist < closest_distance:
                    closest_gap = idx
                    closest_distance = dist

            if not closest_gap is None:
                assert gaps[contig][closest_gap][1] == True
                gaps[contig][closest_gap][1] = False
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

                        if closest_distance is None or dist < closest_distance:
                            closest_gap = idx
                            closest_distance = dist

                    if not closest_gap is None:
                        assert gaps[contig][closest_gap][2] == False
                        gaps[contig][closest_gap][2] = True
                else:
                    print(line[:-1], "not found", file=sys.stderr)
    for contig, gap_list in gaps.items():
        for start, is_closed, was_fixed, source, end, extra in gap_list:
            if was_fixed and not is_closed:
                print(contig, ":", int(start)//1000, "kbp claims to be fixed but is not closed", file=sys.stderr)
            print(contig, source, "fixedgap" if was_fixed else ("closedgap" if is_closed else "gap"), start, end, *extra, sep="\t")


if __name__ == "__main__":
    close_annotation_in_gff(*sys.argv[1:])