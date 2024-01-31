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


def paste_old_sequences(genome_in, gaps_before_closing, gaps_after_closing, masked_in, reversed_out):
    masked = load_masked(masked_in)
    to_undo = {}
    old_gff = {}
    with fileinput.input(gaps_before_closing) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            _, _, _, start, end, _, _, _, extra = line.split("\t")
            for e in extra.split(";"):
                if "ID=" in e:
                    idx = int(e.split("=")[1])
                    old_gff[idx] = (int(start), int(end))
                    break

    with fileinput.input(gaps_after_closing) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            contig, _, what, start, end, _, _, _, extra = line.split("\t")
            idx = None
            for e in extra.split(";"):
                if "ID=" in e:
                    idx = int(e.split("=")[1])
            if what == "gap":
                old_start, old_end = old_gff[idx]
                if contig not in to_undo:
                    to_undo[contig] = []
                to_undo[contig].append((int(start), int(end), contig + " " + str(old_start) + " " + str(old_end)))


    with open(reversed_out, "w") as out_file:
        for contig_name, contig in iterate_contigs(genome_in):
            if contig_name in to_undo:
                to_undo[contig_name].sort(reverse=True)
                for start, end, key in to_undo[contig_name]:
                    if key in masked:
                        paste_sequence = masked[key]
                        # print("pasting a", len(paste_sequence), "bp sequence into", contig_name, "at", start, "-", end, 
                        #       file=sys.stderr)
                        # print("replaced sequence:", contig[start:end], file=sys.stderr)
                        contig = contig[:start-1] + paste_sequence + contig[end:]

                        offset_new_genome = sum([len(masked[key]) - (end - (start - 1)) \
                                                    for start, end, key in to_undo[contig_name] if key in masked])
                        out_file.write("\t".join([contig_name, ".", "undone_masking", 
                                                  str(start + offset_new_genome), 
                                                  str(end + offset_new_genome), ".", ".", "."]) + "\n")
                    else:
                        # print("could not find", contig_name, "at", start, "-", end, "in masked", file=sys.stderr)
                        pass
            print(">" + contig_name)
            for i in range(0, len(contig), 80):
                print(contig[i:i+80])

if __name__ == "__main__":
    paste_old_sequences(*sys.argv[1:])