import csv
import os
import os.path
import re
import sys
from Bio.Blast.Applications import NcbiblastnCommandline
import subprocess
from Bio import SeqIO
from Bio.Emboss.Applications import WaterCommandline
from Bio.Blast import NCBIXML
#increase the field size limit to 500K (so that a whole PTU can fit in)
csv.field_size_limit(500000)

def main():
    repl_and_seq = {}
    for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
        repl_and_seq[seq_record.id] = str(seq_record.seq)
    with open(sys.argv[4], "w") as output_fh:
        #        output_fh.write("##gff-version 3\n")
        for row in csv.reader(open(sys.argv[2]), delimiter="\t"):
            search_seq(row, repl_and_seq, output_fh)


def search_seq(row, repl_and_seq, output_fh):
    features_seq = row[10]
    # replicon_name_start = "_".join(row[9].split("_")[:2]) + "_"
    # E.g. search for "Chr10" and "core"
    replicon_name_reg_ex = "{}_.*{}".format(
        row[9].split("_")[0], row[9].split("_")[1])
    new_replicons_with_same_start = list(filter(
        lambda replicon: re.match(replicon_name_reg_ex, replicon),
        repl_and_seq.keys()))
    assert len(new_replicons_with_same_start) == 1
    new_replicon_name = new_replicons_with_same_start[0]
    ref_seq = repl_and_seq[new_replicon_name]
    # This one is located in a gap won't be found
    #    if row[9] == "Chr7_core_Tb427v3:1909284-1910284":
    #        return
    # First try a simple string match ...
    start_pos, number_of_hits = search_by_exact_string_match(
        ref_seq, features_seq)
    if start_pos != -1 and number_of_hits == 1:
        end_pos = start_pos + len(features_seq)
        _write_entry(row, start_pos, end_pos, new_replicon_name,
                     output_fh)
        return

    # ... if the string match does not work search with blast
    start_pos, end_pos = search_by_blast(
        ref_seq, features_seq, row[9], new_replicon_name, sys.argv[3])
    _write_entry(row, start_pos, end_pos, new_replicon_name,
                 output_fh)


def search_by_blast(ref_seq, query_seq, seq_id_query, seq_id_ref,
                    alignment_folder):
    ref_fasta_file_name = "ref_{}.fa".format(seq_id_query.replace(":", "_"))
    query_fasta_file_name = "query_{}.fa".format(
        seq_id_query.replace(":", "_"))
    blast_db_prefix = "{}/tmp_blast_db".format(alignment_folder)
    blast_output_file = "{}/blast_output_{}.xml".format(
        alignment_folder, seq_id_query)
    with open(ref_fasta_file_name, "w") as ref_fh:
        ref_fh.write(">{}\n{}\n".format(seq_id_ref, ref_seq))
    with open(query_fasta_file_name, "w") as query_fh:
        query_fh.write(">query\n{}\n".format(query_seq))
    subprocess.call(["makeblastdb", "-in", ref_fasta_file_name,
                     "-dbtype", "nucl", "-out", blast_db_prefix])
    blastn_cline = NcbiblastnCommandline(
        query=query_fasta_file_name, db=blast_db_prefix,
        evalue=0.001, out=blast_output_file, outfmt="5")
    stdout, stderr = blastn_cline()
    os.remove("{}.nhr".format(blast_db_prefix))
    os.remove("{}.nin".format(blast_db_prefix))
    os.remove("{}.nsq".format(blast_db_prefix))
    os.remove(ref_fasta_file_name)
    os.remove(query_fasta_file_name)
    blast_record = NCBIXML.read(open(blast_output_file))
    assert len(list(blast_record.alignments)) == 1
    alignment = blast_record.alignments[0]
    if (float(blast_record.query_length) / float(
            alignment.hsps[0].align_length) < 0.95):
        sys.stderr.write("Very low similarity for {} and {}. Stopped! See "
                         "file {}\n".format(seq_id_query, seq_id_ref,
                                            blast_output_file))
        sys.exit(2)
    return alignment.hsps[0].sbjct_start, alignment.hsps[0].sbjct_end
    
            
def _write_entry(row, start_pos, end_pos, new_replicon_name, output_fh):
    row = row[:9]
    row[0] = new_replicon_name
    row[3] = str(start_pos)
    row[4] = str(end_pos)
    output_fh.write("\t".join(row) + "\n")
    

def search_by_exact_string_match(ref_seq, query_seq):
    start_pos = ref_seq.find(query_seq) + 1
    number_of_hits = ref_seq.count(query_seq)
    return start_pos, number_of_hits


def search_by_waterman(ref_seq, query_seq, seq_id, alignment_folder):
    # Warning - using this approach can take a while.
    ref_fasta = "ref_{}.fa".format(seq_id.replace(":", "_"))
    query_fasta = "query_{}.fa".format(seq_id.replace(":", "_"))

    alignment_file_path = "{}/waterman_{}.txt".format(
            alignment_folder, seq_id)
    if not os.path.exists(alignment_file_path):
        with open(ref_fasta, "w") as ref_fh:
            ref_fh.write(">ref\n{}\n".format(ref_seq))
        with open(query_fasta, "w") as query_fh:
            query_fh.write(">query\n{}\n".format(query_seq))
        waterman_cline = WaterCommandline(
            asequence=ref_fasta, bsequence=query_fasta,
            gapopen=10, gapextend=0.5, outfile=alignment_file_path)
        stdout, stderr = waterman_cline()
        os.remove(ref_fasta)
        os.remove(query_fasta)
    alignmet_results = open(alignment_file_path).read()
    gap_perc = None
    start_pos = None
    for line in alignmet_results.split("\n"):
        if line.startswith("# Gaps:"):
            gap_perc = float(line[:-2].split("(")[-1].strip())
        if line.startswith("ref"):
            start_pos = int(line.split()[1])
            break
    if gap_perc > 5.0:
        sys.stderr.write("High gap percentage for {}.\n".format(
            alignment_file_path))
    return start_pos


main()
