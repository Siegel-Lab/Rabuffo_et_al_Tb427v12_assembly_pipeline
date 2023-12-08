import fileinput
import sys

def bin_genome_based_on_gff(gff):
    contig_sizes = {}
    gaps_per_contig = {}
    with fileinput.input(gff) as in_file:
        for line in in_file:
            if line.startswith("##sequence-region"):
                _, contig, _, end = line.split()
                contig_sizes[contig] = int(end)
                continue
            if line[0] == "#":
                continue
            contig, _, _, start, end, _, _, _, *extra = line.split("\t")
            start=int(start)
            end=int(end)
            if contig not in gaps_per_contig:
                gaps_per_contig[contig] = []
            gaps_per_contig[contig].append([start, end])
    for k in contig_sizes.keys():
        v = gaps_per_contig[k] if k in gaps_per_contig else []
        v.sort()
        if len(v) == 0:
            print(k, 1, contig_sizes[k], sep="\t")
            continue
        print(k, 1, v[0][0] - 1, sep="\t")
        for s, e in v:
            print(k, s, e, sep="\t")
        print(k, v[-1][1] + 1, contig_sizes[k], sep="\t")



if __name__ == "__main__":
    bin_genome_based_on_gff(*sys.argv[1:])