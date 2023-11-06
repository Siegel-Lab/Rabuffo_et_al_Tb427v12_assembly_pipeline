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



def get_read_pos_and_strand(reads_in, mates_in):
    stored_reads = {}
    with fileinput.input(reads_in) as in_file2:
        for line in in_file2:
            qname, flag, rname, pos, _, cigar, *_ = line.strip().split()
            clipped = get_clipped_bases(cigar)
            pos = int(pos)
            pos -= clipped
            flag = int(flag)
            rev_stnd = flag & 16
            #print(qname, flag, rname, pos, sep='\t', file=sys.stderr)
            assert not qname in stored_reads
            stored_reads[qname] = [rname, pos, rev_stnd]

    print("#columns:read_name distance_deviation expected_distance chr rpos1 rpos2 strand")
    with fileinput.input(mates_in) as in_file3:
        for line in in_file3:
            qname, flag, rname, pos, _, cigar, *_ = line.strip().split()
            clipped = get_clipped_bases(cigar)
            pos = int(pos)
            pos -= clipped
            flag = int(flag)
            rev_stnd = flag & 16
            if qname in stored_reads:
                o_rname, o_rpos, o_rev_stnd = stored_reads[qname]

                if o_rname != rname:
                    continue
                print(qname, rname, pos, o_rpos, rev_stnd, o_rev_stnd, sep='\t')


if __name__ == "__main__":
    get_read_pos_and_strand(*sys.argv[1:])