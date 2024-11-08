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

def store_reads(reads_in, readlen):
    stored_reads = {}
    with fileinput.input(reads_in) as in_file2:
        for line in in_file2:
            if line[0] == "@":
                continue
            fragname, flag, rname, pos, map_q, cigar, *_ = line.strip().split()
            qname, fragidx = fragname.split("_")
            clipped = get_clipped_bases(cigar)
            fragidx = int(fragidx)
            pos = int(pos)
            map_q = int(map_q)
            pos -= clipped
            flag = int(flag)
            rev_stnd = (flag & 16) == 16
            # sec_or_supp = flag & 2048 or flag & 256 # already filtered out before...
            # if sec_or_supp:
            #     continue
            if not rev_stnd:
                pos -= fragidx * readlen
            else:
                pos += fragidx * readlen
            #print(qname, flag, rname, pos, sep='\t', file=sys.stderr)
            if qname in stored_reads and stored_reads[qname][3] >= map_q:
                continue
            stored_reads[qname] = [rname, pos, rev_stnd, map_q]
    return stored_reads


def get_contig_connecting_reads(reads_in, mates_in, readlen):

    stored_reads = store_reads(reads_in, int(readlen))
    stored_mates = store_reads(mates_in, -int(readlen))

    print("#columns:rname, o_rname, pos, o_rpos, rev_stnd, o_rev_stnd, map_q, o_map_q")
    for qname, (rname, pos, rev_stnd, map_q) in stored_mates.items():
        if qname in stored_reads:
            o_rname, o_rpos, o_rev_stnd, o_map_q = stored_reads[qname]

            if o_rname == rname or "*" in [o_rname, rname]:
                continue

            print(rname, o_rname, pos, o_rpos, rev_stnd, o_rev_stnd, map_q, o_map_q, sep='\t')


if __name__ == "__main__":
    get_contig_connecting_reads(*sys.argv[1:])