import fileinput
import sys


class TXX:
    def __init__(self, idx, start, end, anno_type, strnd, line):
        self.outgoing_right = False
        self.outgoing_left = False
        self.incoming_right = False
        self.incoming_left = False
        self.idx = idx
        self.strnd = strnd
        self.anno_type = anno_type
        self.start = start
        self.end = end
        self.line = line

    def center(self):
        return (self.start + self.end) // 2
    
    def should_go_right(self):
        return self.anno_type == "dTSS" or (self.strnd == "+" and self.anno_type == "sTSS")
    def should_go_left(self):
        return self.anno_type == "dTSS" or (self.strnd == "-" and self.anno_type == "sTSS")
    def should_come_right(self):
        return self.anno_type == "cTTS" or (self.strnd == "-" and self.anno_type == "sTTS")
    def should_come_left(self):
        return self.anno_type == "cTTS" or (self.strnd == "+" and self.anno_type == "sTTS")
    
    def forw_strand(self):
        return self.anno_type == "dTSS" or self.anno_type == "cTTS" or self.strnd == "+"
    def rev_strand(self):
        return self.anno_type == "dTSS" or self.anno_type == "cTTS" or self.strnd == "-"
    
    def is_as_expected(self):
        return (self.should_go_right() and self.outgoing_right) or \
               (self.should_go_left() and self.outgoing_left) or \
               (self.should_come_right() and self.incoming_right) or \
               (self.should_come_left() and self.incoming_left)


def main(file_in):
    txx_per_contig = {}
    txx_by_idx = {}
    with fileinput.input(file_in) as in_file:
        for line in in_file:
            if line[0] == "#":
                print(line[:-1])
                continue
            contig, _, anno_type, start, end, _, strnd, _, *extra = line[:-1].split("\t")
            if not anno_type in ["dTSS", "sTSS", "cTTS", "sTTS"]:
                continue
            if contig not in txx_per_contig:
                txx_per_contig[contig] = ([], [])
            idx = "None"
            for e in extra:
                if "ID=" in e:
                    idx = e.split("=")[1]
                    break
            txx = TXX(idx, int(start), int(end), anno_type, strnd, line[:-1])
            txx_by_idx[idx] = txx
            if txx.forw_strand():
                txx_per_contig[contig][0].append(txx)
            if txx.rev_strand():
                txx_per_contig[contig][1].append(txx)

    for contig, (txxs_f, txxs_r) in txx_per_contig.items():
        for strnd, txxs in [("+", txxs_f), ("-", txxs_r)]:
            txxs.sort(key=lambda txx: txx.start)
            for txx_a, txx_b in zip(txxs[:-1], txxs[1:]):
                # if txx_a.should_go_right() != txx_b.should_come_left():
                #     print("Warning: unexpected TXX:\n", txx_a.line, "\n", txx_b.line, file=sys.stderr)
                #     # assert txx_a.should_go_right() == txx_b.should_come_left()
                # if txx_a.should_come_right() != txx_b.should_go_left():
                #     print("Warning: unexpected TXX:\n", txx_a.line, "\n", txx_b.line, file=sys.stderr)
                #     # assert txx_a.should_come_right() == txx_b.should_go_left()

                if txx_a.should_go_right() and txx_b.should_come_left() and strnd == "+":
                    print(contig, "Siegel_Group", "PTU", txx_a.center(), txx_b.center(), ".", "+", ".", 
                        "TSS_ID=" + txx_a.idx + ";TSS_type=" + txx_a.anno_type + ";TTS_ID=" + txx_b.idx + ";TTS_type=" + txx_b.anno_type, sep="\t")
                    txx_by_idx[txx_a.idx].outgoing_right = True
                    txx_by_idx[txx_b.idx].incoming_left = True
                if txx_a.should_come_right() and txx_b.should_go_left() and strnd == "-":
                    print(contig, "Siegel_Group", "PTU", txx_a.center(), txx_b.center(), ".", "-", ".", 
                        "TSS_ID=" + txx_b.idx + ";TSS_type=" + txx_b.anno_type + ";TTS_ID=" + txx_a.idx + ";TTS_type=" + txx_a.anno_type, sep="\t")
                    txx_by_idx[txx_a.idx].incoming_right = True
                    txx_by_idx[txx_b.idx].outgoing_left = True

    for _, txx in txx_by_idx.items():
        if not txx.is_as_expected():
            print("Warning: unexpected TXX:", txx.line, 
                txx.outgoing_right, txx.outgoing_left, txx.incoming_right, txx.incoming_left,
                file=sys.stderr)




if __name__ == "__main__":
    main(*sys.argv[1:])