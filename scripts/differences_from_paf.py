import fileinput
import sys

def turn_paf_into_differences(paf_file):
    matches_by_contig = {}
    with fileinput.input(paf_file) as in_file:
        for line in in_file:
            q_name, q_len, q_start, q_end, strand, r_name, r_len, r_start, r_end, num_matches, align_len, mapq, *extra = line.strip().split()

            if q_name != r_name:
                continue
            if strand != '+':
                continue
            if q_name not in matches_by_contig:
                matches_by_contig[q_name] = []
            r_start = int(r_start)
            r_end = int(r_end)
            q_start = int(q_start)
            q_end = int(q_end)
            num_matches = int(num_matches)
            align_len = int(align_len)
            mapq = int(mapq)
            matches_by_contig[q_name].append((q_start, q_end, r_start, r_end, num_matches, align_len, mapq))

    for contig in matches_by_contig:
        matches_by_contig[contig].sort(key=lambda x: x[0])
        translocation_detected = False
        for a, b in zip(matches_by_contig[contig][:-1], matches_by_contig[contig][1:]):
            a_q_start, a_q_end, a_r_start, a_r_end, a_num_matches, a_align_len, a_mapq = a
            b_q_start, b_q_end, b_r_start, b_r_end, b_num_matches, b_align_len, b_mapq = b

            # make sure no translocations ocurred
            if a_q_start >= b_q_start:
                translocation_detected = True
                break

        if not translocation_detected:
            for a, b in zip(matches_by_contig[contig][:-1], matches_by_contig[contig][1:]):
                a_q_start, a_q_end, a_r_start, a_r_end, a_num_matches, a_align_len, a_mapq = a
                b_q_start, b_q_end, b_r_start, b_r_end, b_num_matches, b_align_len, b_mapq = b

                if a_q_end > b_q_start:
                    # the alignments are overlapping on the reference.
                    # fix it by including the duplicated region in the change
                    dup_reg_size = a_q_end - b_q_start
                    print("diff", dup_reg_size)
                    tmp = b_q_start
                    b_q_start = a_q_end
                    a_q_end = tmp
                    a_r_end -= dup_reg_size
                    b_r_start += dup_reg_size
                    print(a, b)
                    
                # print change in paf format
                print(contig, b_q_start - a_q_end, a_q_end, b_q_start, "+", contig, b_r_start - a_r_end, a_r_end, b_r_start, 0, (b_q_start - a_q_end) + (b_r_start - a_r_end), 0, sep='\t')


if __name__ == "__main__":
    turn_paf_into_differences(sys.argv[1])