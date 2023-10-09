import fileinput
import sys

def get_average_distance_deviation(reads_in, mates_in, expected_distances):
    expected_distances = {}
    with fileinput.input(expected_distances) as in_file:
        for line in in_file:
            if line.strip().split() == 2:
                name, distance = line.strip().split()
                expected_distances[name] = int(distance)

    stored_reads = {}
    with fileinput.input(reads_in) as in_file:
        for line in in_file:
            qname, flag, rname, pos, *_ = line.strip().split()
            pos = int(pos)
            flag = int(flag)
            is_mate = flag % 64 != 0
            assert qname not in stored_reads
            stored_reads[qname] = [rname, pos, is_mate]
    print("#columns:read_name distance_deviation")
    with fileinput.input(mates_in) as in_file:
        for line in in_file:
            qname, flag, rname, pos, *_ = line.strip().split()
            pos = int(pos)
            flag = int(flag)
            is_mate = flag % 64 != 0
            if qname in stored_reads:
                o_rname, o_rpos, o_is_mate = stored_reads[qname]

                if o_is_mate == is_mate:
                    continue
                if o_rname != rname:
                    continue
                
                dist = pos - o_rpos
                if is_mate:
                    dist = -dist
                print(qname, dist - expected_distances[qname], sep='\t')


if __name__ == "__main__":
    get_average_distance_deviation(*sys.argv[1:])