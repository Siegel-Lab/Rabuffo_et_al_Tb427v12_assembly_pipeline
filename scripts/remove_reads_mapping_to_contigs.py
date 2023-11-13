import fileinput
import sys

def remove_reads_mapping_to_contigs(reads_in, alignments_in, contigs_in):
    remove_contigs = set()
    with fileinput.input(files=(contigs_in)) as contigs:
        for line in contigs:
            remove_contigs.add(line[:-1].strip())

    reads_to_remove = set()
    with fileinput.input(files=(alignments_in)) as alignments:
        for line in alignments:
            line = line[:-1].strip().split('\t')
            if line[2] in remove_contigs:
                reads_to_remove.add(line[0])

    with fileinput.input(files=(reads_in)) as reads:
        for lines in zip(*[reads]*4):
            lines = [l.strip() for l in lines]
            name, sequence, _, qual = lines

            if name[1:] not in reads_to_remove:
                print(name)
                print(sequence)
                print("+")
                print(qual)


if __name__ == '__main__':
    remove_reads_mapping_to_contigs(*sys.argv[1:])