import fileinput
import sys

def get_clipped_bases(cigar):
    first_symbol_idx = 0
    while first_symbol_idx < len(cigar) and cigar[first_symbol_idx].isdigit():
        first_symbol_idx += 1
    if first_symbol_idx == len(cigar):
        return 0
    if not cigar[first_symbol_idx] in ["S", "H"]:
        return 0
    return int(cigar[:first_symbol_idx])



def get_average_distance_deviation(reads_in, mates_in, expected_distances_file):
    expected_distances = {}
    with fileinput.input(expected_distances_file) as in_file:
        for line in in_file:
            if len(line.strip().split()) == 2:
                name, distance = line.strip().split()
                expected_distances[name[1:]] = int(distance)

    stored_reads = {}
    with fileinput.input(reads_in) as in_file2:
        for line in in_file2:
            qname, flag, rname, pos, _, cigar, *_ = line.strip().split()
            clipped = get_clipped_bases(cigar)
            pos = int(pos)
            pos -= clipped
            flag = int(flag)
            #print(qname, flag, rname, pos, sep='\t', file=sys.stderr)
            assert not qname in stored_reads
            stored_reads[qname] = [rname, pos]

    print("#columns:read_name distance_deviation")
    with fileinput.input(mates_in) as in_file3:
        for line in in_file3:
            qname, flag, rname, pos, _, cigar, *_ = line.strip().split()
            clipped = get_clipped_bases(cigar)
            pos = int(pos)
            pos -= clipped
            flag = int(flag)
            if qname in stored_reads and qname in expected_distances:
                o_rname, o_rpos = stored_reads[qname]

                if o_rname != rname:
                    continue

                dist = pos - o_rpos
                print(qname, dist - expected_distances[qname], sep='\t')


if __name__ == "__main__":
    get_average_distance_deviation(*sys.argv[1:])