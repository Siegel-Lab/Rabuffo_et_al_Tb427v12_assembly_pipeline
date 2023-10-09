import fileinput
import sys


def fixup_number_of_n(fasta_in, num_in, num_out):
    str_in = ""
    for _ in range(num_in):
        str_in += "n"
    str_out = ""
    for _ in range(num_out):
        str_out += "n"
    with fileinput.input(fasta_in) as in_file:
        for line in in_file:
            if line[0] == ">":
                print(line[:-1])
                continue
            
            print(line.lower().replace(str_in, str_out)[:-1])
