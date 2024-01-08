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


def load_gff(filename):
    ret = {}
    with fileinput.input(filename) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            fields = line[:-1].split("\t")
            if fields[0] not in ret:
                ret[fields[0]] = []
            ret[fields[0]].append(fields)
    return ret

def mask_regions(genome_in, gff_in, masked_out, annotation_in, annotation_out, num_n=1000):
    annotations = load_gff(annotation_in)
    contig_names = []
    to_mask = {}
    with open(masked_out, "w") as out_file:
        with open(annotation_out, "w") as out_annotation:
            print("##gff-version 3", file=out_annotation)
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
                contig_names.append(contig_name)
                if contig_name in to_mask:
                    to_mask[contig_name].sort(reverse=True)
                    for idx, (start, end) in enumerate(to_mask[contig_name]):
                        sequence_to_remove_before = sum( (end - start) - num_n for start, end in to_mask[contig_name][idx+1:])
                        out_pos = start - sequence_to_remove_before
                        out_file.write(">" + contig_name + " " + str(out_pos + 1) + " " + str(out_pos + num_n) + "\n")
                        out_file.write(contig[start:end] + "\n")
                        contig = contig[:start] + "N" * num_n + contig[end:]
                print(">" + contig_name)
                print("##sequence-region", contig_name, 1, len(contig), sep="\t", file=out_annotation)
                for i in range(0, len(contig), 80):
                    print(contig[i:i+80])
                
            for contig_name in contig_names:
                if contig_name in to_mask:
                    for anno in annotations[contig_name]:
                        anno[3] = int(anno[3]) # start
                        anno[4] = int(anno[4]) # end
                        if all(anno[4] < sg or anno[3] > eg for sg, eg in to_mask[contig_name]):
                            to_move = sum((eg - sg) - num_n for sg, eg in to_mask[contig_name] if eg < anno[3])
                            anno[3] -= to_move
                            anno[4] -= to_move
                            print(*anno, sep="\t", file=out_annotation)
                else:
                    for anno in annotations[contig_name]:
                        print(*anno, sep="\t", file=out_annotation)


if __name__ == "__main__":
    mask_regions(*sys.argv[1:])