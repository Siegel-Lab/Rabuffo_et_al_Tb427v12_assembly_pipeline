import fileinput
import sys

def sort_gff(gff_file, order_in):
    order = {}
    with fileinput.input(order_in) as in_file:
        for idx, line in enumerate(in_file):
            order[line.strip()] = idx
    gff_header = []
    gff_body = []
    with fileinput.input(gff_file) as in_file:
        for line in in_file:
            if line.startswith("##sequence-region"):
                _, contig, _, end = line.strip().split()
                if not contig in order:
                    print("WARNING: contig", contig, "not found in order file. Annotation of this contig is dropped", gff_file, file=sys.stderr)
                    continue
                gff_header.append((order[contig], contig, int(end)))
            if line[0] == "#":
                continue
            contig_name, _, anno_type, start, *_ = line.strip().split()
            if not contig_name in order:
                print("WARNING: contig", contig_name, "not found in order file. Annotation of this contig is dropped", gff_file, file=sys.stderr)
                continue
            gff_body.append((order[contig_name], int(start), anno_type, line.strip()))

    gff_header.sort()
    gff_body.sort()

    print("##gff-version 3")
    for _, contig, end in gff_header:
        print("##sequence-region " + contig + " 1 " + str(end))
    for _, _, _, line in gff_body:
        print(line)


if __name__ == "__main__":
    sort_gff(*sys.argv[1:])