import sys

def load_dist_dev(file_name, max_dist=None):
    names = []
    ret = []
    with open(file_name, "r") as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            else:
                # read_name distance_deviation expected_distance chr rpos1 rpos2 strand
                readname, distance, expected, chr, pos1, pos2, strnd, map_q = line.strip().split()
                if max_dist is None or int(distance) < max_dist:
                    names.append(readname)
                    chr = chr[:-len("_Tb427v10")]
                    ret.append([int(distance), int(expected), chr, int(pos1), int(pos2), strnd, int(map_q)])
    return names, ret


def load_gaps(file_name):
    gap_pos = {}
    contig_sizes = {}

    with open(file_name, "r") as in_file:
        for line in in_file:
            if line[0] == "#":
                if line.startswith("##sequence-region"):
                    contig, _, end = line.strip().split()[1:]
                    contig = contig[:-len("_Tb427v10")]
                    contig_sizes[contig] = int(end)
                continue
            contig, source, anno_type, start, end, *extra = line.strip().split()
            contig = contig[:-len("_Tb427v10")]
            
            gap_name = contig + ":" + str(int(start)//1000) + "kbp"
            gap_pos[gap_name] = [contig, int(start), int(end)]
    return gap_pos, contig_sizes



def post_process(ref_names, ref_dev, min_dev=-100, min_map_q=30, max_dev=-float("inf")):
    data = []
    for r_name, (distance, expected, chr, pos1, pos2, strnd, map_q) in zip(ref_names, ref_dev):
        if distance < min_dev and map_q >= min_map_q and max_dev < distance:
            data.append([chr, min(pos1, pos2), max(pos1, pos2), distance, r_name])

    return data


def percentile(l, p=0.05):
    l.sort()
    return l[int(len(l) * p)]

def cluster(data, distance_y=500, min_reads=5, max_cluster_size=50000):
    data.sort(key=lambda x: (x[0], x[1], x[3]))
    clusters = []
    for chr, start, end, deviation, r_name in data:
        fitting_clusters = []
        for idx, cluster in enumerate(clusters):
            fits_in_cluster = False
            for cluster_chr, cluster_start, cluster_end, cluster_deviation, _ in cluster:
                if cluster_chr == chr and end >= cluster_start and start <= cluster_end and abs(deviation - cluster_deviation) <= distance_y:
                    fitting_clusters.append(idx)
                    break
        if len(fitting_clusters) == 0:
            clusters.append([[chr, start, end, deviation, r_name]])
        else:
            clusters[fitting_clusters[0]].append([chr, start, end, deviation, r_name])
            for idx in fitting_clusters[1:]:
                clusters[fitting_clusters[0]] += clusters[idx]
                clusters[idx] = []

    processed_clusters = []
    for cluster in clusters:
        if len(cluster) > min_reads and not "unitig" in cluster[0][0]:
            cluster_start = percentile([x[1] for x in cluster], 0.95)
            cluster_end = percentile([x[2] for x in cluster], 0.05)
            if abs(cluster_start - cluster_end) > max_cluster_size: # filter out extremely large clusters
                continue
            if cluster_start > cluster_end:
                continue
            #assert cluster_start < cluster_end
            cluster_deviation = sum([x[3] for x in cluster]) / len(cluster)
            #print("cluster from", cluster_start, "to", cluster_end, "with deviation", cluster_deviation, "and", 
            #    len(cluster), "reads")
            processed_clusters.append([cluster[0][0], cluster_start, cluster_end, cluster_deviation, cluster])

    return processed_clusters

def filter_clusters_with_counter_indication(clusters, data, min_indication=5):
    ret = []
    for cluster_chr, cluster_start, cluster_end, cluster_deviation, c in clusters:
        counter_indication = 0
        for chr, start, end, deviation, r_name in data:
            if chr == cluster_chr and deviation >= 0 and deviation <= 200:
                if start <= cluster_start and end >= cluster_end:
                    counter_indication += 1
        if counter_indication < min_indication:
            ret.append([cluster_chr, cluster_start, cluster_end, cluster_deviation, c])


    return ret

def filter_clusters_that_overlap_gap(clusters, gap_pos, min_distance_to_gap=10000):
    # filter out clusters that overlap with gaps
    cluster_overlapping_gap = []
    new_cluster = []
    for cluster_chr, cluster_start, cluster_end, cluster_deviation, c in clusters:
        found_gap = False
        for gap_name, (gap_chr, gap_start, gap_end) in gap_pos.items():
            if cluster_chr == gap_chr and cluster_start <= gap_end + min_distance_to_gap and cluster_end + min_distance_to_gap >= gap_start:
                cluster_overlapping_gap.append([cluster_chr, cluster_start, cluster_end, cluster_deviation, c])
                found_gap = True
                break
        if not found_gap:
            new_cluster.append([cluster_chr, cluster_start, cluster_end, cluster_deviation, c])
    return new_cluster, cluster_overlapping_gap


def has_overlapping_clusters(clusters, cut_back_distance=1000):
    for a, b in zip(clusters[:-1], clusters[1:]):
        if a[0] == b[0] and a[2] + cut_back_distance*3 >= b[1]:
            return True
    return False

def merge_overlapping_clusters_helper(new_cluster, contig_sizes, cut_back_distance=1000):
    ret = []
    for cluster_chr, cluster_start, cluster_end, cluster_deviation, c in new_cluster:
        if cluster_start < cut_back_distance*2:
            continue
        if cluster_end > contig_sizes[cluster_chr] - cut_back_distance*2:
            continue
        if len(ret) == 0 or ret[-1][0] != cluster_chr or ret[-1][2] + cut_back_distance*3 < cluster_start:
            ret.append([cluster_chr, cluster_start, cluster_end, cluster_deviation, c])
        else:
            ret[-1][2] = max(ret[-1][2], cluster_end)
            ret[-1][4] = ret[-1][4] + c
            ret[-1][3] = sum([x[3] for x in ret[-1][4]]) / len(ret[-1][4])
    return ret

def merge_overlapping_clusters(new_cluster, contig_sizes, cut_back_distance=1000):
    new_cluster.sort(key=lambda x: (x[0], x[1]))

    while has_overlapping_clusters(new_cluster, cut_back_distance):
        new_cluster = merge_overlapping_clusters_helper(new_cluster, contig_sizes, cut_back_distance=cut_back_distance)
    return new_cluster


def main(distance_deviation_filename, reference_gaps_filename, missassemblies_gff_filename, min_dev=-100, max_dev=-float("inf"), filter_overlap=True, filter_non_overlap=False, region_name="misassembly", cut_back_distance=1000):
    if isinstance(min_dev, str):
        min_dev = float(min_dev)
    if isinstance(max_dev, str):
        max_dev = float(max_dev)
    if isinstance(filter_overlap, str):
        filter_overlap = filter_overlap == "True"
    if isinstance(filter_non_overlap, str):
        filter_non_overlap = filter_non_overlap == "True"
    ref_names, ref_dev = load_dist_dev(distance_deviation_filename)
    gap_pos, contig_sizes = load_gaps(reference_gaps_filename)
    data = post_process(ref_names, ref_dev, min_dev=min_dev, max_dev=max_dev)

    clusters = cluster(data)
    clusters = filter_clusters_with_counter_indication(clusters, data)
    if filter_overlap:
        kept_clusters, filtered_clusters = filter_clusters_that_overlap_gap(clusters, gap_pos)
        # print("clusters overlapping gaps")
        # for cluster_chr, cluster_start, cluster_end, cluster_deviation, c in filtered_clusters:
        #     print(cluster_chr, int(cluster_start/1000), "-", int(cluster_end/1000), "k")
        # print()
    elif filter_non_overlap:
        filtered_clusters, kept_clusters = filter_clusters_that_overlap_gap(clusters, gap_pos)
        # print("clusters not overlapping gaps")
        # for cluster_chr, cluster_start, cluster_end, cluster_deviation, c in filtered_clusters:
        #     print(cluster_chr, int(cluster_start/1000), "-", int(cluster_end/1000), "k")
        # print()
    else:
        kept_clusters = clusters
    kept_clusters = merge_overlapping_clusters(kept_clusters, contig_sizes, cut_back_distance=cut_back_distance)


    # print("kept clusters")
    # for cluster_chr, cluster_start, cluster_end, cluster_deviation, c in kept_clusters:
    #     print(cluster_chr, int(cluster_start/1000), "-", int(cluster_end/1000), "k")
    # print()

    with open(missassemblies_gff_filename, "w") as file_out:
        for cluster_chr, cluster_start, cluster_end, cluster_deviation, c in kept_clusters:
            file_out.write("\t".join([cluster_chr + "_Tb427v10", ".", region_name, str(cluster_start), str(cluster_end), ".", ".", ".", ""]) + "\n" )

    # if filter_overlap:
    #     gap_without_cluster = {}
    #     for gap_name, (gap_chr, gap_start, gap_end) in gap_pos.items():
    #         found_gap = False
    #         for cluster_chr, cluster_start, cluster_end, cluster_deviation, c in filtered_clusters:
    #             if cluster_chr == gap_chr and cluster_start <= gap_end and cluster_end >= gap_start:
    #                 found_gap = True
    #                 break
    #         if not found_gap:
    #             gap_without_cluster[gap_name] = [gap_chr, gap_start, gap_end]

    #     print("gaps without cluster")
    #     for g in gap_without_cluster.keys():
    #         print(g)
    #     print()

if __name__ == "__main__":
    main(*sys.argv[1:])