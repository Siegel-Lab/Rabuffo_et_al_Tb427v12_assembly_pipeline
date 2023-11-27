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

def extract_region(genome_in, gff_in):
    to_mask = {}
    with fileinput.input(gff_in) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            contig, _, _, start, end, _, _, _, extra = line.split("\t")
            start = int(start)
            end = int(end)
            name = None
            for e in extra.split(";"):
                if "Name=" in e:
                    name = e.split("=")[1]
            if contig not in to_mask:
                to_mask[contig] = []
            to_mask[contig].append((start, end, name))
    for contig_name, contig in iterate_contigs(genome_in):
        if contig_name in to_mask:
            to_mask[contig_name].sort()
            for start, end, name in to_mask[contig_name]:
                print(">" + name)
                for i in range(start, end, 80):
                    print(contig[i:i+80])

if __name__ == "__main__":
    extract_region(*sys.argv[1:])