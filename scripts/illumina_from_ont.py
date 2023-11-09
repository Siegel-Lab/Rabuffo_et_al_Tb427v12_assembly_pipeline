import fileinput
import sys

def illumina_from_ont(fasta_in, fasta_out, mate_out, stats_out, min_ont_read_len, illumina_read_len, num_illumina_reads):
    min_ont_read_len = int(min_ont_read_len)
    illumina_read_len = int(illumina_read_len)
    num_illumina_reads = int(num_illumina_reads)

    with fileinput.input(fasta_in) as in_file:
        with open(fasta_out, 'w') as out_file:
            with open(mate_out, 'w') as mate_file:
                with open(stats_out, 'w') as stats_file:
                    for lines in zip(*[in_file]*4):
                        lines = [l.strip() for l in lines]
                        name, sequence, _, qual = lines
                        name = name.strip().split()[0]
                        if len(sequence) < min_ont_read_len:
                            continue

                        for idx in range(num_illumina_reads):
                            out_file.write(name + "_" + str(idx) + "\n" + 
                                           sequence[illumina_read_len*idx:illumina_read_len*(idx+1)] + "\n+\n" + 
                                        qual[illumina_read_len*idx:illumina_read_len*(idx+1)] + "\n")
                        for idx in range(num_illumina_reads):
                            mate_file.write(name + "_" + str(idx) + "\n" + 
                                            sequence[-illumina_read_len*(idx+1):-illumina_read_len*(idx)] + "\n+\n" + 
                                            qual[-illumina_read_len*(idx+1):-illumina_read_len*(idx)] + "\n")
                        stats_file.write(name + "\t" + str(len(sequence) - illumina_read_len) + "\n")


if __name__ == "__main__":
    illumina_from_ont(*sys.argv[1:])