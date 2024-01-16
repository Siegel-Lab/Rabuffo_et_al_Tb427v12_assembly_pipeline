import fileinput
import sys

def extract_ids(gff_line, actor):    
    if "Tb427_" in gff_line:
        split = gff_line.split("Tb427_")
        for kdx, ele in enumerate(split[1:]):
            a = ele[0]
            b = ele[1:]
            idx = len(b)
            for c in ";,:.":
                if c in b:
                    jdx = b.index(c)
                    if jdx != -1 and jdx < idx:
                        idx = jdx
            b2 = b[:idx]
            c = b[idx:]
            # print(ele, a, b2, c, file=sys.stderr)
            split[kdx+1] = a + str(actor(int(b2))) + c

        return "Tb427_".join(split)
    else:
        return gff_line

def load_annos(file):
    ret = {}
    with fileinput.input(file) as in_file:
        for line in in_file:
            if line[0] == "#":
                continue
            ctg, _, _, start, end, *_ = line.strip().split()
            if ctg not in ret:
                ret[ctg] = []
            ret[ctg].append((int(start), int(end), line.strip()))
    return ret

def extract_annotation(annotation_file, gap_file, merge_annos_file):
    annos = load_annos(annotation_file)
    gaps = load_annos(gap_file)
    merge_annos = load_annos(merge_annos_file)
    ids = set()
    def add_to_ids(idx):
        ids.add(idx)
        return idx
    for x in merge_annos.values():
        for _, _, line in x:
            extract_ids(line, add_to_ids)
    
    max_id = max(ids)
    previous_id = {}

    idx_counter = max_id
    def fix_idx(idx):
        if idx in previous_id:
            return previous_id[idx]
        else:
            nonlocal idx_counter
            idx_counter += 100
            previous_id[idx] = idx_counter
            assert idx_counter not in ids
            return idx_counter
    for ctg, x in annos.items():
        for start, end, line in x:
            keep = False
            if ctg in gaps:
                for a, b, _ in gaps[ctg]:
                    if a <= end and start <= b and (a < start or end < b):
                        keep = True
                        break
                if keep:
                    print(extract_ids(line, fix_idx))


if __name__ == "__main__":
    extract_annotation(*sys.argv[1:])
