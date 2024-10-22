import fileinput
import sys


def main(gff_file, expand_region_size=1000):
    print("#columns: contig, type, start, idx, size_change, size_before, size_after")
    with fileinput.input(gff_file) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            contig, _, waht, start, end, _, _, _, extra = line[:-1].split("\t")
            extra = extra.split(";")
            idx = None
            prev_size=1000
            for e in extra:
                if "ID=" in e:
                    idx = e.split("=")[1]
                if "previous_size=" in e:
                    prev_size = int(e.split("=")[1])
            start = int(start)
            end = int(end)
            size = 1 + end - start - 2*expand_region_size
            print(contig, waht, start, idx, size-prev_size, prev_size + 2*expand_region_size, size + 2*expand_region_size, sep="\t")

if __name__ == "__main__":
    main(*sys.argv[1:])