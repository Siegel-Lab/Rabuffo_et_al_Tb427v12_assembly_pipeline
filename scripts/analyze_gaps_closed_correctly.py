from identify_collapsed_regions import load_dist_dev, load_gaps
import sys
import re

def load_gap_spanning_reads(gap_spanning_reads_file):
    gap_spanning_reads = {}

    with open(gap_spanning_reads_file, "r") as file_in:
        for line in file_in:
            readname, *gaps = line[:-1].strip().split()
            gap_names = []
            for gap in gaps:
                chrom, start, _ = re.split(":|-", gap)
                #chrom = chrom[:-len("_Tb427v10")]
                gap_name = chrom + ":" + str(int(start)//1000) + "kbp"
                gap_names.append(gap_name)
            gap_spanning_reads[readname] = " ".join(gap_names)
    return gap_spanning_reads


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

def filter_clusters(clusters):
    return [c for c in clusters if len(c) > 3]

def get_mean_deviation_in_clusters(clusters, ref_dict, fixed_dict, in_fixed=True):
    return [sum([fixed_dict[n] if in_fixed else ref_dict[n] for n in cluster]) / len(cluster) for cluster in clusters]

def make_read_dicts(x_names, x_dev):
    x_dict = {x: y[0] for x, y in zip(x_names, x_dev)}

    return x_dict

def filter_reads(fixed_names, fixed_dev, gap_spanning_reads, ref_dict):
    # filter reads
    # remove those where the disrance has not changed
    filtered = set()
    for read_name, distance in zip(fixed_names, fixed_dev):
        if not read_name in gap_spanning_reads:
            filtered.add(read_name)
            # if read_name in ref_dict and ref_dict[read_name] == distance:
            #     filtered.add(read_name)
        pass
    readnames = [n for n in fixed_names if n in ref_dict and not n in filtered]
    return readnames




def analyze_gaps_closed_correctly(dist_dev_ref_file, dist_dev_fixed_file, gap_spanning_reads_file, gaps_old_genome, 
                                  gap_closed_if_fixed_dev_smaller_than=5000, 
                                  max_ref_fixed_diff = 10, max_ref_fixed_sum = 1000):
    ref_names, ref_dev = load_dist_dev(dist_dev_ref_file)
    fixed_names, fixed_dev = load_dist_dev(dist_dev_fixed_file)
    gap_spanning_reads = load_gap_spanning_reads(gap_spanning_reads_file)

    ref_dict = make_read_dicts(ref_names, ref_dev)
    fixed_dict = make_read_dicts(fixed_names, fixed_dev)

    readnames = filter_reads(fixed_names, fixed_dev, gap_spanning_reads, ref_dict)

    gap_pos_old_genome = load_gaps(gaps_old_genome)

    gap_names = gap_pos_old_genome.keys()

    for gap in sorted(gap_names):
        read_names = extract_reads_for_gap(gap, gap_spanning_reads, readnames)
        chrom, start, end = gap_pos_old_genome[gap]
        if len(read_names) > 0:
            read_clusters = filter_clusters(cluster_reads(read_names, ref_dict, fixed_dict,
                                                          max_ref_fixed_diff=max_ref_fixed_diff, 
                                                            max_ref_fixed_sum=max_ref_fixed_sum))
            cluster_fixed = get_mean_deviation_in_clusters(read_clusters, ref_dict, fixed_dict)
            gap_closed = False
            min_fixed = float("inf")
            for x in cluster_fixed:
                if abs(x) < gap_closed_if_fixed_dev_smaller_than and abs(x) < min_fixed:
                    gap_closed = True
                    min_fixed = abs(x)
            if gap_closed:
                print(chrom + "_Tb427v10", ".", "fixedgap", str(start), str(end), ".", ".", ".", 
                      "estimated_length=1000;gap_type=within scaffold;closed_correctly=true", sep="\t")
            else:
                print(chrom + "_Tb427v10", ".", "gap", str(start), str(end), ".", ".", ".", 
                      "estimated_length=1000;gap_type=within scaffold;not_enough_data=true", sep="\t")
        else:
            print(chrom + "_Tb427v10", ".", "gap", str(start), str(end), ".", ".", ".", 
                  "estimated_length=1000;gap_type=within scaffold;not_enough_data=true", sep="\t")


if __name__ == "__main__":
    analyze_gaps_closed_correctly(*sys.argv[1:])