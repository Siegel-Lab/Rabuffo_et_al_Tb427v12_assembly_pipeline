import fileinput
import sys

def check_spans_gap(reads_in, mates_in, gff_file):
    gaps = {}
    with fileinput.input(gff_file) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            if len(line.strip().split()) >= 5:
                contig_name, _, anno_type, gap_start, gap_end, *_ = line.strip().split()
                if anno_type == "gap":
                    if not contig_name in gaps:
                        gaps[contig_name] = []
                    gaps[contig_name].append((int(gap_start), int(gap_end)))

    for contig_name in gaps.keys():
        gaps[contig_name].sort()

    stored_reads = {}
    with fileinput.input(reads_in) as in_file2:
        for line in in_file2:
            qname, flag, rname, pos, *_ = line.strip().split()
            pos = int(pos)
            flag = int(flag)
            assert not qname in stored_reads
            stored_reads[qname] = [rname, pos]

    with fileinput.input(mates_in) as in_file3:
        for line in in_file3:
            qname, flag, rname, pos, *_ = line.strip().split()
            pos = int(pos)
            flag = int(flag)
            if qname in stored_reads:
                o_rname, o_rpos = stored_reads[qname]

                if o_rname != rname:
                    continue

                spans_gap = []
                if rname in gaps:
                    for gap_start, gap_end in gaps[rname]:
                        if (o_rpos < gap_start and pos > gap_end) or (pos < gap_start and o_rpos > gap_end):
                            spans_gap.append(rname[:-len("_Tb427v10")] + ":" + str(int(gap_start/1000)) + "kbp")

                if len(spans_gap) > 0:
                    print(qname, *spans_gap, sep='\t')


if __name__ == "__main__":
    check_spans_gap(*sys.argv[1:])