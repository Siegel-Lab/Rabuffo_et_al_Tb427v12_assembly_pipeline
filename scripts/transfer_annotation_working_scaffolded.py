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
        contigs[contig_name] = contig
    return contigs

def transfer_annotation_working_scaffold(to_transfer, anno_scaffolded):
    contig_translations = {}
    with fileinput.input(anno_scaffolded) as in_file:
        for line in in_file:
            if line[0] == "#" or len(line) <= 1:
                continue
            cols = line.strip().split("\t")
            if len(cols) == 9:
                ctg, a, anno_type, start, e, b, c, d, extra = cols
            elif len(cols) == 8:
                ctg, a, anno_type, start, e, b, c, d = cols
                extra = ""
            else:
                raise RuntimeError("Unknown format " + str(len(cols)) + " " + line)
            name=None
            if "contig" in anno_type:
                for n in extra.split(";"):
                    if "Name=" in n:
                        name = n.split("=")[1]
                        break
                if name == None:
                    print("no name", line, file=sys.stderr)
                    assert name != None
                # print(ctg, name, file=sys.stderr)
                if "_B" in name and "core" in name:
                    continue
                name = "_".join(name.split("_")[:3])
                if name in contig_translations:
                    contig_translations[name][1] = min(contig_translations[name][1], int(start))
                else:
                    contig_translations[name] = [ctg, int(start)]



    with fileinput.input(to_transfer) as in_file:
        for line in in_file:
            if line[0] == "#":
                print(line, end="")
                continue
            ctg, x, anno_type, start, end, *extra = line.strip().split("\t")
            if ctg not in contig_translations:
                print("not found", ctg, file=sys.stderr)
                print(contig_translations.keys(), file=sys.stderr)
            move_by = contig_translations[ctg][1] - 1
            start = int(start) + move_by
            end = int(end) + move_by
            ctg = contig_translations[ctg][0]
            print(ctg, x, anno_type, start, end, *extra, sep="\t")



if __name__ == '__main__':
    transfer_annotation_working_scaffold(*sys.argv[1:])