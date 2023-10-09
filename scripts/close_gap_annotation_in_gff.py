import fileinput
import sys

def close_annotation_in_gff(annotation, remaining_gaps):
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
                gaps[contig].append([int(start), True, source, end, extra])
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
    for contig, gap_list in gaps.items():
        for start, is_closed, source, end, extra in gap_list:
            print(contig, source, "closedgap" if is_closed else "gap", start, end, *extra, sep="\t")


if __name__ == "__main__":
    close_annotation_in_gff(*sys.argv[1:])