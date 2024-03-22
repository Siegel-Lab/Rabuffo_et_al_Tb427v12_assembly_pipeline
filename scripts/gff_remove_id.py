import fileinput
import sys

def gff_remove_id(gff_file):
    with fileinput.input(gff_file) as in_file:
        for line in in_file:
            eles = line.strip().split("\t")
            if not line.startswith("#"):
                eles[-1] = ";".join([x for x in eles[-1].split(";") if not x.startswith("ID=")])
            print("\t".join(eles))


if __name__ == "__main__":
    gff_remove_id(*sys.argv[1:])