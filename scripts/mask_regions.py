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

def mask_regions(genome_in, gff_in, num_n=1000):
    to_mask = {}
    with fileinput.input(gff_in) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            contig, _, _, start, end, *_ = line.split("\t")
            start = int(start)
            end = int(end)
            if contig not in to_mask:
                to_mask[contig] = []
            to_mask[contig].append((start, end))
    for contig_name, contig in iterate_contigs(genome_in):
        if contig_name in to_mask:
            to_mask[contig_name].sort(reverse=True)
            for start, end in to_mask[contig_name]:
                contig = contig[:start] + "N" * num_n + contig[end:]
        print(">" + contig_name)
        for i in range(0, len(contig), 80):
            print(contig[i:i+80])

if __name__ == "__main__":
    mask_regions(*sys.argv[1:])