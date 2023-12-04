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

def extract_region(genome_in, gff_in, annotation_in, annotation_out):
    to_extract = {}
    annos_by_contig = {}
    with fileinput.input(annotation_in) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            contig = line.split("\t")[0]
            if contig not in annos_by_contig:
                annos_by_contig[contig] = []
            annos_by_contig[contig].append(line[:-1].split("\t"))
    with fileinput.input(gff_in) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            contig, _, _, start, end, _, _, _, extra = line[:-1].split("\t")
            start = int(start)
            end = int(end)
            name = None
            for e in extra.split(";"):
                if "Name=" in e:
                    name = e.split("=")[1]
            if contig not in to_extract:
                to_extract[contig] = []
            to_extract[contig].append((start, end, name))
    with open(annotation_out, "w") as out_file:
        out_file.write("##gff-version 3\n")
        for contig_name, contig in iterate_contigs(genome_in):
            if contig_name in to_extract:
                to_extract[contig_name].sort()
                for start, end, name in to_extract[contig_name]:
                    out_file.write("##sequence-region " + name + " 1 " + str(1 + end - start) + "\n")
                    print(">" + name)
                    for i in range(start, end, 80):
                        print(contig[i:i+80])
        for contig_name, contig in iterate_contigs(genome_in):
            if contig_name in to_extract:
                for start, end, name in to_extract[contig_name]:
                    for _, x1, x2, anno_start, anno_end, x3, x4, x5, extra in annos_by_contig[contig_name]:
                        anno_start = int(anno_start)
                        anno_end = int(anno_end)
                        if anno_start >= start and anno_end <= end:
                            anno_start -= start
                            anno_end -= start
                            out_file.write("\t".join([name, x1, x2, str(anno_start + 1), str(anno_end + 1), x3, x4, x5, extra]) + "\n")


if __name__ == "__main__":
    extract_region(*sys.argv[1:])