import sys
from get_average_distance_deviation import store_reads

def get_read_pos_and_strand(reads_in, mates_in, readlen):
    stored_reads = store_reads(reads_in, int(readlen))
    stored_mates = store_reads(mates_in, -int(readlen))

    print("#columns:read_name distance_deviation expected_distance chr rpos1 rpos2 strand")
    for qname, (rname, pos, rev_stnd, map_q) in stored_mates.items():
        if qname in stored_reads:
            o_rname, o_rpos, o_rev_stnd, o_map_q = stored_reads[qname]

            if o_rname != rname:
                continue
            print(qname, rname, pos, o_rpos, rev_stnd, o_rev_stnd, min(map_q, o_map_q), sep='\t')


if __name__ == "__main__":
    get_read_pos_and_strand(*sys.argv[1:])