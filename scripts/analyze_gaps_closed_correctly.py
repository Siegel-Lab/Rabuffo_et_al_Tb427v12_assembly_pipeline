from identify_collapsed_regions import load_dist_dev
import sys
import re

def load_gap_spanning_reads(gap_spanning_reads_file):
    gap_spanning_reads = {}

    with open(gap_spanning_reads_file, "r") as file_in:
        for line in file_in:
            readname, *gaps = line[:-1].strip().split()
            readname = readname.split("_")[0]
            gap_names = []
            for gap in gaps:
                chrom, start, _ = re.split(":|-", gap)
                # chrom = chrom[:-len("_Tb427v10")]
                gap_name = chrom + ":" + str(int(start)//1000) + "kbp"
                gap_names.append(gap_name)
            gap_spanning_reads[readname] = " ".join(gap_names)
    return gap_spanning_reads

def load_gaps(file_name):
    gap_pos = {}
    contig_sizes = {}

    with open(file_name, "r") as in_file:
        for line in in_file:
            if line[0] == "#":
                if line.startswith("##sequence-region"):
                    contig, _, end = line.strip().split()[1:]
                    contig_sizes[contig] = int(end)
                continue
            contig, source, anno_type, start, end, *extra = line.strip().split()
            contig = "_".join(contig.split("_")[:-1])
            
            gap_name = contig + ":" + str(int(start)//1000) + "kbp"
            gap_pos[gap_name] = [contig, int(start), int(end), anno_type]
    return gap_pos, contig_sizes

def extract_reads_for_gap(gap_name, gap_spanning_reads, readnames):
    return [n for n in readnames if n in gap_spanning_reads and gap_spanning_reads[n] == gap_name]

def cluster(l, key, max_dif):
    l.sort(key=key)

    clustered = [[l[0]]]
    if len(l) > 1:
        for read in l[1:]:
            if abs(key(read) - key(clustered[-1][-1])) > max_dif:
                clustered.append([])
            clustered[-1].append(read)
    return clustered

def cluster_reads(readnames, ref_dict, fixed_dict, max_ref_fixed_diff, max_ref_fixed_sum):
    return [c for l in cluster(readnames, lambda x: ref_dict[x] - fixed_dict[x], max_ref_fixed_diff) 
               for c in cluster(l, lambda x: ref_dict[x] + fixed_dict[x], max_ref_fixed_sum)]

def filter_clusters(clusters, min_reads):
    return [c for c in clusters if len(c) > min_reads]

def get_mean_deviation_in_clusters(clusters, ref_dict, fixed_dict, in_fixed=True):
    return [sum([fixed_dict[n] if in_fixed else ref_dict[n] for n in cluster]) / len(cluster) for cluster in clusters]

def make_read_dicts(x_names, x_dev, min_map_q):
    x_dict = {x: y[0] for x, y in zip(x_names, x_dev) if y[-1] >= min_map_q}

    return x_dict

def filter_reads(fixed_names, fixed_dev, gap_spanning_reads, ref_dict, fixed_dict):
    # filter reads
    # remove those where the disrance has not changed
    filtered = set()
    for read_name, distance in zip(fixed_names, fixed_dev):
        if not read_name in gap_spanning_reads:
            filtered.add(read_name)
            # if read_name in ref_dict and ref_dict[read_name] == distance:
            #     filtered.add(read_name)
        pass
    readnames = [n for n in fixed_names if n in ref_dict and n in fixed_dict and not n in filtered]
    return readnames




def analyze_gaps_closed_correctly(dist_dev_ref_file, dist_dev_fixed_file, gap_spanning_reads_file, gaps_new_genome, 
                                  gap_closed_if_fixed_dev_smaller_than=5000,
                                  max_ref_fixed_diff = 10, max_ref_fixed_sum = 1000, min_map_q=30,
                                  min_reads_in_cluster=2 # @todo change back to 5, 2 is just a test
                            ):
    ref_names, ref_dev = load_dist_dev(dist_dev_ref_file)
    fixed_names, fixed_dev = load_dist_dev(dist_dev_fixed_file)
    gap_spanning_reads = load_gap_spanning_reads(gap_spanning_reads_file)
    # print(gap_spanning_reads.values(), file=sys.stderr)

    ref_dict = make_read_dicts(ref_names, ref_dev, min_map_q)
    fixed_dict = make_read_dicts(fixed_names, fixed_dev, min_map_q)

    readnames = filter_reads(fixed_names, fixed_dev, gap_spanning_reads, ref_dict, fixed_dict)
    #print(readnames, file=sys.stderr)

    gap_pos_new_genome, _ = load_gaps(gaps_new_genome)

    gap_names = gap_pos_new_genome.keys()
    # print(gap_names, file=sys.stderr)

    for gap in sorted(gap_names):
        # print(gap, file=sys.stderr)
        read_names = extract_reads_for_gap(gap, gap_spanning_reads, readnames)
        chrom, start, end, gap_type = gap_pos_new_genome[gap]
        if gap_type == "gap":
            print(chrom + "_Tb427v11", ".", "gap", str(start), str(end), ".", ".", ".", sep="\t")
        else:
            if len(read_names) >= min_reads_in_cluster:
                read_clusters = filter_clusters(cluster_reads(read_names, ref_dict, fixed_dict,
                                                            max_ref_fixed_diff=max_ref_fixed_diff, 
                                                                max_ref_fixed_sum=max_ref_fixed_sum),
                                                min_reads=min_reads_in_cluster)
                cluster_fixed = get_mean_deviation_in_clusters(read_clusters, ref_dict, fixed_dict)
                gap_closed = False
                min_fixed = float("inf")
                for x in cluster_fixed:
                    if abs(x) < gap_closed_if_fixed_dev_smaller_than and abs(x) < min_fixed:
                        gap_closed = True
                        min_fixed = abs(x)
                if gap_closed:
                    print(chrom + "_Tb427v11", ".", "fixed", str(start), str(end), ".", ".", ".", sep="\t")
                else:
                    print(chrom + "_Tb427v11", ".", "failed", str(start), str(end), ".", ".", ".", sep="\t")
            else:
                print(chrom + "_Tb427v11", ".", "no_data", str(start), str(end), ".", ".", ".", sep="\t")


if __name__ == "__main__":
    analyze_gaps_closed_correctly(*sys.argv[1:])