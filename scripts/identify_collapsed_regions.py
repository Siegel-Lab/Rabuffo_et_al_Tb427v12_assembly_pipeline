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

    with open(file_name, "r") as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            contig, source, anno_type, start, end, *extra = line.strip().split()
            contig = contig[:-len("_Tb427v10")]
            
            gap_name = contig + ":" + str(int(start)//1000) + "kbp"
            gap_pos[gap_name] = [contig, int(start), int(end)]
    return gap_pos



def post_process(ref_names, ref_dev, min_dev = -100):
    data = []
    for r_name, (distance, expected, chr, pos1, pos2, strnd, mapq_q) in zip(ref_names, ref_dev):
        if distance < -100:
            data.append([chr, min(pos1, pos2), max(pos1, pos2), distance, r_name])

    return data


def percentile(l, p=0.05):
    l.sort()
    return l[int(len(l) * p)]

def cluster(data, distance_y=500, min_reads=5, max_cluster_size=50000):
    data.sort(key=lambda x: (x[0], x[1], x[3]))
    clusters = []
    for chr, start, end, deviation, r_name in data:
        fits_in_cluster = False
        for cluster in clusters:
            for cluster_chr, cluster_start, cluster_end, cluster_deviation, _ in cluster:
                if cluster_chr == chr and end >= cluster_start and start <= cluster_end and abs(deviation - cluster_deviation) <= distance_y:
                    fits_in_cluster = True
                    break
            if fits_in_cluster:
                cluster.append([chr, start, end, deviation, r_name])
                break
        if not fits_in_cluster:
            clusters.append([])
        clusters[-1].append([chr, start, end, deviation, r_name])

    processed_clusters = []
    for cluster in clusters:
        if len(cluster) > min_reads and not "unitig" in cluster[0][0]:
            cluster_start = percentile([x[1] for x in cluster], 0.05)
            cluster_end = percentile([x[1] for x in cluster], 0.95)
            if abs(cluster_start - cluster_end) > max_cluster_size: # filter out extremely large clusters
                continue
            assert cluster_start < cluster_end
            cluster_deviation = sum([x[3] for x in cluster]) / len(cluster)
            #print("cluster from", cluster_start, "to", cluster_end, "with deviation", cluster_deviation, "and", 
            #    len(cluster), "reads")
            processed_clusters.append([cluster[0][0], cluster_start, cluster_end, cluster_deviation, cluster])

    return processed_clusters



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


def has_overlapping_clusters(clusters):
    for a, b in zip(clusters[:-1], clusters[1:]):
        if a[0] == b[0] and a[2] >= b[1]:
            return True
    return False

def merge_overlapping_clusters_helper(new_cluster):
    ret = []
    for cluster_chr, cluster_start, cluster_end, cluster_deviation, c in new_cluster:
        if len(ret) == 0 or ret[-1][0] != cluster_chr or ret[-1][2] < cluster_start:
            ret.append([cluster_chr, cluster_start, cluster_end, cluster_deviation, c])
        else:
            ret[-1][2] = max(ret[-1][2], cluster_end)
            ret[-1][4] = ret[-1][4] + c
            ret[-1][3] = sum([x[3] for x in ret[-1][4]]) / len(ret[-1][4])
    return ret

def merge_overlapping_clusters(new_cluster):
    new_cluster.sort(key=lambda x: (x[0], x[1]))

    while has_overlapping_clusters(new_cluster):
        new_cluster = merge_overlapping_clusters_helper(new_cluster)
    return new_cluster

if __name__ == "__main__":
    distance_deviation_filename = sys.argv[1] # "../data/out/virtual_paired_read_dist/referece.distance_deviation"
    reference_gaps_filename = sys.argv[2] # "../data/out/samba_out_1/reference.gaps.gff3"
    missassemblies_gff_filename = sys.argv[3] # "../data/in/analysis_in/mis_assemblies.gff"

    ref_names, ref_dev = load_dist_dev(distance_deviation_filename)
    gap_pos = load_gaps(reference_gaps_filename)
    data = post_process(ref_names, ref_dev)

    clusters = cluster(data)
    new_cluster, cluster_overlapping_gap = filter_clusters_that_overlap_gap(clusters, gap_pos)
    new_cluster = merge_overlapping_clusters(new_cluster)

    print("clusters overlapping gaps")
    for cluster_chr, cluster_start, cluster_end, cluster_deviation, c in cluster_overlapping_gap:
        print(cluster_chr, int(cluster_start/1000), "-", int(cluster_end/1000), "k")
    print()

    print("new clusters")
    for cluster_chr, cluster_start, cluster_end, cluster_deviation, c in new_cluster:
        print(cluster_chr, int(cluster_start/1000), "-", int(cluster_end/1000), "k")
    print()

    with open(missassemblies_gff_filename, "w") as file_out:
        for cluster_chr, cluster_start, cluster_end, cluster_deviation, c in new_cluster:
            file_out.write("\t".join([cluster_chr + "_Tb427v10", ".", "misassembly", str(cluster_start), str(cluster_end), 
                                    ".", ".", ".", ""]) + "\n" )

    gap_without_cluster = {}
    for gap_name, (gap_chr, gap_start, gap_end) in gap_pos.items():
        found_gap = False
        for cluster_chr, cluster_start, cluster_end, cluster_deviation, c in cluster_overlapping_gap:
            if cluster_chr == gap_chr and cluster_start <= gap_end and cluster_end >= gap_start:
                found_gap = True
                break
        if not found_gap:
            gap_without_cluster[gap_name] = [gap_chr, gap_start, gap_end]

    print("gaps without cluster")
    for g in gap_without_cluster.keys():
        print(g)
    print()
