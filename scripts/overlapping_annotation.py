# usage: python overlapping_annotation.py <input_file> <from_annotation> <to_annotation>

# you need python 3 & the intervaltree package to run this script

import sys
from intervaltree import Interval, IntervalTree

from_anno_list = {}
to_anno_list = {}

from_annos = sys.argv[2].split(",")
to_annos = sys.argv[3].split(",")
add_distance = 0
if len(sys.argv) > 4:
    add_distance = int(sys.argv[4])

with open(sys.argv[1], "r") as in_file:
    for line in in_file:
        if len(line) == 0 or line[0] == "#":
            continue
        else:
            chrom, _, anno, start, end, *extra = line[:-1].split()
            if anno in from_annos:
                if not anno in from_anno_list:
                    from_anno_list[anno] = []
                from_anno_list[anno].append([chrom, int(start), int(end), line[:-1]])
            if anno in to_annos:
                if anno not in to_anno_list:
                    to_anno_list[anno] = {}
                if chrom not in to_anno_list[anno]:
                    to_anno_list[anno][chrom] = IntervalTree()
                start = int(start)
                end = int(end)
                if end >= start:
                    to_anno_list[anno][chrom].add(Interval(start, end))

print("add distance:", add_distance)

to_anno_keys = list(to_anno_list.keys())

print("", end="\t")
for to_anno in to_anno_keys:
    print(to_anno, end="\t")
print("total")

for from_anno, l in from_anno_list.items():
    print(from_anno, end="\t")
    for to_anno in to_anno_keys:
        overlap = 0
        for chrom, start, end, _ in l:
            if chrom in to_anno_list[to_anno]:
                if start - add_distance < end + add_distance:
                    overlap += min(1, len(to_anno_list[to_anno][chrom].overlap(start - add_distance, 
                                                                               end + add_distance)))
        print(overlap, end="\t")
    print(len(l))

    


# py ../scripts/overlapping_annotation.py ../data/out/21_overview_of_remaining_gaps/gff_last_column_removed_2.gff Centromere,rRNA,tRNA contig_subt,contig_core,unitig,gap,closedgap_full,closedgap_core,closedgap_a,closedgap_b,closedgap_masked,expanded_region,unexpanded_reg 10000
# add distance: 10000
#   gap     closedgap_full  closedgap_a     closedgap_b     expanded_region unexpanded_reg  closedgap_masked        contig_subt     contig_core     unitig  total
#rRNA            0       34      20      0       0       0       0       12      64      44      120
#tRNA            0       0       0       0       3       7       0       1       130     0       131
#Centromere      5       10      0       0       0       0       1       5       16      0       21



# min distance: 10kb for gap & repeat. 0kb for contig
#            open_gap  closed_gap  collapsed_repeat  fixed_repeat  contig_subt  contig_core  unitig  total
#rRNA        0         54          0                 0             12           64           44      120
#tRNA        0         0           7                 3             1            130          0       131
#Centromere  5         10          0                 0             5            16           0       21
#improvement:          ^^^^^^^^^^                    ^^^^^^^^^^^^