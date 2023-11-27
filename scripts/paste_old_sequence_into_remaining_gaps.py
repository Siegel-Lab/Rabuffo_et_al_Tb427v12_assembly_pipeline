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

def load_masked(fasta_in):
    ret = {}
    with fileinput.input(fasta_in) as in_file:
        for header, sequence in zip(*[in_file]*2):
            sequence = sequence[:-1]
            ret[header[1:-1]] = sequence
    return ret


def paste_old_sequences(genome_in, old_gff_in, gff_in, masked_in):
    masked = load_masked(masked_in)
    to_undo = {}
    old_gff = {}
    with fileinput.input(old_gff_in) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            _, _, _, start, end, _, _, _, extra = line.split("\t")
            for e in extra.split(";"):
                if "ID=" in e:
                    idx = int(e.split("=")[1])
                    old_gff[idx] = (int(start), int(end))
                    break

    with fileinput.input(gff_in) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            contig, _, what, start, end, _, _, _, extra = line.split("\t")
            idx = None
            for e in extra.split(";"):
                if "ID=" in e:
                    idx = int(e.split("=")[1])
            if what == "gap":
                start, end = old_gff[idx]
                if contig not in to_undo:
                    to_undo[contig] = []
                to_undo[contig].append((start, end))

    for contig_name, contig in iterate_contigs(genome_in):
        if contig_name in to_undo:
            to_undo[contig_name].sort(reverse=True)
            for start, end in to_undo[contig_name]:
                key = contig_name + " " + str(start) + " " + str(end)
                if key in masked:
                    paste_sequence = masked[contig_name + " " + str(start) + " " + str(end)]
                    # print("pasting a", len(paste_sequence), "bp sequence into", contig_name, "at", start, "-", end, 
                    #       file=sys.stderr)
                    # print("replaced sequence:", contig[start:end], file=sys.stderr)
                    contig = contig[:start] + paste_sequence + contig[end:]
        print(">" + contig_name)
        for i in range(0, len(contig), 80):
            print(contig[i:i+80])

if __name__ == "__main__":
    paste_old_sequences(*sys.argv[1:])