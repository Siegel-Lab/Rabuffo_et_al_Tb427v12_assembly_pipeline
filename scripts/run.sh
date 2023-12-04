#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH --mem=500G
#SBATCH --time=7-00:00:00
#SBATCH --partition=fat
#SBATCH --job-name=ont_try_assembly
#SBATCH -o slurm_out/ont_assembly-%j.out


setup() {
    source /home/mschmidt/.miniconda3/etc/profile.d/conda.sh

    conda activate ont_assembly

    set -e

    module load ngs/minimap2/2.10
    module load ngs/samtools/1.9
    module load ngs/deeptools/3.5.0
    module load ngs/bedtools2/2.28.0

    GENOME_FOLDER_IN=$(realpath ../data/in/genome_in/HGAP3_Tb427v10_diploid)
    GENOME_FILENAME_IN="HGAP3_Tb427v10_diploid_scaffolded"
    GFF_IN_XXX=$(realpath ../data/in/genome_in/HGAP3_Tb427v10_diploid/HGAP3_Tb427v10_diploid_scaffolded.gff3)
    ANA_LYSIS_IN=$(realpath ../data/in/analysis_in)
    REF_CENTRO=$(realpath ../data/in/genome_in/HGAP3_Tb427v10/HGAP3_Tb427v10.fasta)
    GFF_CENTRO_IN=$(realpath ../data/in/genome_in/HGAP3_Tb427v10/HGAP3_Tb427v10_manual.gff3)
    ONT_READS_IN=$(realpath ../data/in/ont_reads_in/merged.nanopore.gz)
    
    BIN_DIR=$(realpath ../bin/)
    DATA_DIR=$(realpath ../data)
    SCRIPTS_DIR=$(realpath .)


    mkdir -p ${DATA_DIR}/out
}

#
# @todo section:
#
# - check for reads that connect different contigs
#
# - remember to ask raul for the final v11 assembly
# - pick data from anna and what raul sent
#



main(){
    setup

    # move the Centromere annotation from the assembly where the cores are merged to the fully phased assembly
    # @todo this does not work with an exact match...
    # move_annotation ${DATA_DIR}/out/1_move_centro_anno \
    #                 "Centromere" \
    #                 ${GFF_CENTRO_IN} \
    #                 ${REF_CENTRO} \
    #                 ${GENOME_FOLDER_IN}/${GENOME_FILENAME_IN}.fasta \
    #                 ${GFF_IN_XXX}
    #                 # -> ${DATA_DIR}/out/1_move_centro_anno/annotation_combined.gff

    remove_annotated_gaps ${DATA_DIR}/out/1_remove_gap_annotation \
                          ${GFF_IN_XXX}
                          # -> ${DATA_DIR}/out/1_remove_gap_annotation/annotation.gapless.gff3

    # reannotate the gaps in the fully phased assembly
    annotate_gaps ${DATA_DIR}/out/2_ref_reannotated_gaps \
                  ${GENOME_FOLDER_IN}/${GENOME_FILENAME_IN}.fasta
                  # -> ${DATA_DIR}/out/2_ref_reannotated_gaps/gaps.gff3

    # create a picture
    generate_overview_pic ${DATA_DIR}/out/3_overview_of_gaps \
                          ${DATA_DIR}/out/2_ref_reannotated_gaps/gaps.gff3 \
                          "gap" \
                          "gap=purple"

    mask_and_close ${DATA_DIR}/out/4_close_gaps_full_genome \
                   ${GENOME_FOLDER_IN}/${GENOME_FILENAME_IN}.fasta \
                   ${ONT_READS_IN} \
                   ${DATA_DIR}/out/1_remove_gap_annotation/annotation.gapless.gff3 \
                   ${GENOME_FOLDER_IN}/${GENOME_FILENAME_IN}.fasta 
                    # -> ${DATA_DIR}/out/4_close_gaps_full_genome/7.1_undo_failed_masking/masking_undone.fasta

    split_genome_in_a_and_b ${DATA_DIR}/out/5_split_genome \
                            ${DATA_DIR}/out/4_close_gaps_full_genome/7.1_undo_failed_masking/masking_undone.fasta \
                            ${ONT_READS_IN}
                            # -> ${DATA_DIR}/out/5_split_genome/A.fasta
                            # -> ${DATA_DIR}/out/5_split_genome/B.fasta
                            # -> ${DATA_DIR}/out/5_split_genome/remainder.fasta
                            # -> ${DATA_DIR}/out/5_split_genome/reads.A.fasta
                            # -> ${DATA_DIR}/out/5_split_genome/reads.B.fasta

    mask_and_close ${DATA_DIR}/out/6_closed_gaps_a \
                   ${DATA_DIR}/out/5_split_genome/A.fasta \
                   ${DATA_DIR}/out/5_split_genome/reads.A.fasta \
                   ${DATA_DIR}/out/1_remove_gap_annotation/annotation.gapless.gff3 \
                   ${GENOME_FOLDER_IN}/${GENOME_FILENAME_IN}.fasta 
                    # -> ${DATA_DIR}/out/6_closed_gaps_a/7.1_undo_failed_masking/masking_undone.fasta

    mask_and_close ${DATA_DIR}/out/7_closed_gaps_b \
                   ${DATA_DIR}/out/5_split_genome/B.fasta \
                   ${DATA_DIR}/out/5_split_genome/reads.B.fasta \
                    ${DATA_DIR}/out/1_remove_gap_annotation/annotation.gapless.gff3 \
                   ${GENOME_FOLDER_IN}/${GENOME_FILENAME_IN}.fasta
                    # -> ${DATA_DIR}/out/7_closed_gaps_b/7.1_undo_failed_masking/masking_undone.fasta


    merge_genomes ${DATA_DIR}/out/8_merged_genomes \
                  ${DATA_DIR}/out/6_closed_gaps_a/7.1_undo_failed_masking/masking_undone.fasta \
                  ${DATA_DIR}/out/7_closed_gaps_b/7.1_undo_failed_masking/masking_undone.fasta \
                  ${DATA_DIR}/out/5_split_genome/remainder.fasta
                 # -> ${DATA_DIR}/out/8_merged_genomes/assembly.fasta
                 # -> ${DATA_DIR}/out/8_merged_genomes/annotation.gff

    transfer_annotation ${DATA_DIR}/out/8.1_transfer_annotation \
                        ${DATA_DIR}/out/8_merged_genomes/assembly.fasta \
                        ${DATA_DIR}/out/1_remove_gap_annotation/annotation.gapless.gff3 \
                        ${GENOME_FOLDER_IN}/${GENOME_FILENAME_IN}.fasta \
                        ${DATA_DIR}/out/8_merged_genomes/gaps.gff3
                        # -> ${DATA_DIR}/out/8.1_transfer_annotation/annotation_combined.gff
                        # -> ${DATA_DIR}/out/8.1_transfer_annotation/annotation.transfered.gff

    transfer_fixed_regions ${DATA_DIR}/out/8.2_transfer_fixed_regions \
                        "filledgap\|filledmasked;closedgap_full;${DATA_DIR}/out/4_close_gaps_full_genome/7.1_undo_failed_masking \
                         filledgap\|filledmasked;closedgap_a;${DATA_DIR}/out/6_closed_gaps_a/7.1_undo_failed_masking \
                         filledgap\|filledmasked;closedgap_b;${DATA_DIR}/out/7_closed_gaps_b/7.1_undo_failed_masking" \
                        ${DATA_DIR}/out/8.1_transfer_annotation/annotation_combined.gff \
                        ${DATA_DIR}/out/8_merged_genomes/assembly.fasta
                        # -> ${DATA_DIR}/out/8.2_transfer_fixed_regions/annotation_combined.gff
                        # -> ${DATA_DIR}/out/8.2_transfer_fixed_regions/combined.transfered.gff

    generate_overview_pic ${DATA_DIR}/out/9_overview_of_remaining_gaps \
                          ${DATA_DIR}/out/8.2_transfer_fixed_regions/annotation_combined.gff \
                          "gene closedgap_full closedgap_a closedgap_b gap" \
                          "gene=lightgrey;closedgap_full=green;closedgap_a=green;closedgap_b=green;gap=purple"

    gap_spanning_reads ${DATA_DIR}/out/10_gap_spanning_reads \
                    ${DATA_DIR}/out/8_merged_genomes/assembly.fasta \
                    ${ONT_READS_IN} \
                    ${DATA_DIR}/out/8_merged_genomes/gaps.gff3
                    # -> ${DATA_DIR}/out/10_gap_spanning_reads/distance_deviation.tsv


    # mask_region ${DATA_DIR}/out/10.1_test_masked_repeats \
    #             ${DATA_DIR}/out/8_merged_genomes/assembly.fasta \
    #             ${DATA_DIR}/in/mask_repeats/manual_mask.gff
    #             # -> ${DATA_DIR}/out/10.1_test_masked_repeats/masked.fasta
                
    # gap_spanning_reads ${DATA_DIR}/out/10.2_test_gap_spanning_reads \
    #                 ${DATA_DIR}/out/10.1_test_masked_repeats/joined.fasta \
    #                 ${ONT_READS_IN} \
    #                 ${DATA_DIR}/in/mask_repeats/manual_mask.gff

    # annotate_gaps ${DATA_DIR}/out/10.3_test_annotate_gaps \
    #             ${DATA_DIR}/out/10.1_test_masked_repeats/joined.fasta
    #             # -> ${DATA_DIR}/out/16_reannotated_gaps/gaps.gff3

    # align_reads_to_genome ${DATA_DIR}/out/10.4_test_align_reads \
    #     ${DATA_DIR}/out/10.1_test_masked_repeats/joined.fasta \
    #     ${ONT_READS_IN}


    identify_collapsed_regions ${DATA_DIR}/out/13_identify_collapsed_regions \
        ${DATA_DIR}/out/10_gap_spanning_reads/distance_deviation.tsv \
        "#\|closedgap_full\|closedgap_a\|closedgap_b\|gap" \
        ${DATA_DIR}/out/8.2_transfer_fixed_regions/annotation_combined.gff
        # -> ${DATA_DIR}/out/13_identify_collapsed_regions/annotation.gff
        # -> ${DATA_DIR}/out/13_identify_collapsed_regions/collapsed_regions.gff

    generate_overview_pic ${DATA_DIR}/out/14_overview_collapsed_repeats \
                          ${DATA_DIR}/out/13_identify_collapsed_regions/annotation.gff \
                          "gene misassembly closedgap_full closedgap_a closedgap_b gap" \
                          "gene=lightgrey;misassembly=pink;closedgap_full=green;closedgap_a=green;closedgap_b=green;gap=purple"

    mask_region ${DATA_DIR}/out/15_masked_repeats \
                ${DATA_DIR}/out/8_merged_genomes/assembly.fasta \
                ${DATA_DIR}/out/13_identify_collapsed_regions/collapsed_regions.gff
                # -> ${DATA_DIR}/out/15_masked_repeats/masked.fasta

    annotate_gaps ${DATA_DIR}/out/16_reannotated_gaps \
                ${DATA_DIR}/out/15_masked_repeats/masked.fasta
                # -> ${DATA_DIR}/out/16_reannotated_gaps/gaps.gff3

    generate_overview_pic ${DATA_DIR}/out/17_overview_of_masked_regions \
                        ${DATA_DIR}/out/16_reannotated_gaps/gaps.gff3 \
                        "gene gap" \
                        "gene=lightgrey;gap=purple"

    close_gaps ${DATA_DIR}/out/18_closed_gaps \
            ${DATA_DIR}/out/15_masked_repeats \
            masked \
            ${ONT_READS_IN} \
            ${DATA_DIR}/out/16_reannotated_gaps/gaps.gff3
            # -> ${DATA_DIR}/out/18_closed_gaps/gaps.gff3
            # -> ${DATA_DIR}/out/18_closed_gaps/assembly.fasta

    # GENOME_AFTER_CLOSING=$2
    # GAPS_BEFORE_CLOSING=$3
    # GAPS_AFTER_CLOSING=$4
    # MASKED_SEQUENCES=$5
    # GENOME_BEFORE_CLOSING=$6
    undo_masking ${DATA_DIR}/out/18.1_undo_failed_masking \
                ${DATA_DIR}/out/18_closed_gaps/assembly.fasta \
                ${DATA_DIR}/out/16_reannotated_gaps/gaps.gff3 \
                ${DATA_DIR}/out/18_closed_gaps/gaps.gff3 \
                ${DATA_DIR}/out/15_masked_repeats/removed_sequences.fasta \
                ${DATA_DIR}/out/15_masked_repeats/masked.fasta
                # -> ${DATA_DIR}/out/18.1_undo_failed_masking/masking_undone.fasta
                # -> ${DATA_DIR}/out/18.1_undo_failed_masking/gaps.gff3
                # -> ${DATA_DIR}/out/18.1_undo_failed_masking/fixed_gaps_and_masked.gff3

    annotate_gaps ${DATA_DIR}/out/18.2_reannotated_gaps \
                ${DATA_DIR}/out/18.1_undo_failed_masking/masking_undone.fasta
                # -> ${DATA_DIR}/out/18.2_reannotated_gaps/gaps.gff3


    transfer_annotation ${DATA_DIR}/out/19_transfer_annotation \
                        ${DATA_DIR}/out/18_closed_gaps/assembly.fasta \
                        ${DATA_DIR}/out/1_remove_gap_annotation/annotation.gapless.gff3 \
                        ${GENOME_FOLDER_IN}/${GENOME_FILENAME_IN}.fasta \
                        ${DATA_DIR}/out/18.2_reannotated_gaps/gaps.gff3
                        # -> ${MaC_OUT_FOLDER}/19_transfer_annotation/annotation.transfered.gff

    transfer_fixed_regions ${DATA_DIR}/out/20_transfer_fixed_regions \
                        "filledgap\|filledmasked;closedgap_full;${DATA_DIR}/out/4_close_gaps_full_genome/7.1_undo_failed_masking \
                         filledgap\|filledmasked;closedgap_a;${DATA_DIR}/out/6_closed_gaps_a/7.1_undo_failed_masking \
                         filledgap\|filledmasked;closedgap_b;${DATA_DIR}/out/7_closed_gaps_b/7.1_undo_failed_masking \
                         filledmasked;expanded_region;${DATA_DIR}/out/18.1_undo_failed_masking \
                         filledgap;closedgap_masked;${DATA_DIR}/out/18.1_undo_failed_masking" \
                        ${DATA_DIR}/out/19_transfer_annotation/annotation_combined.gff \
                        ${DATA_DIR}/out/18.1_undo_failed_masking/masking_undone.fasta
                        # -> ${DATA_DIR}/out/20_transfer_fixed_regions/annotation_combined.gff


    generate_overview_pic ${DATA_DIR}/out/21_overview_of_remaining_gaps \
                        ${DATA_DIR}/out/20_transfer_fixed_regions/annotation_combined.gff \
                        "gene filledgap closedgap_full closedgap_a closedgap_b expanded_region closedgap_masked gap" \
                        "gene=lightgrey;closedgap_full=green;closedgap_a=green;closedgap_b=green;closedgap_masked=green;expanded_region=blue;gap=purple"

    generate_overview_pic ${DATA_DIR}/out/21.1_overview_of_untransferred_annotations \
                        ${DATA_DIR}/out/19_transfer_annotation/annotation.failed.gff \
                        "gene" \
                        "gene=lightgrey"

    gap_spanning_reads ${DATA_DIR}/out/22_vpr_new_genome \
                                 ${DATA_DIR}/out/18.1_undo_failed_masking/masking_undone.fasta \
                                 ${ONT_READS_IN} \
                                 ${DATA_DIR}/out/18.2_reannotated_gaps/gaps.gff3
                                 # -> ${DATA_DIR}/out/22_vpr_new_genome/distance_deviation.tsv

    # assembly_gaps_individually ${DATA_DIR}/out/22.1_individually_assembled_gaps \
    #         ${DATA_DIR}/out/22_vpr_new_genome/gap_spanning_reads.tsv \
    #         ${ONT_READS_IN}
    #         # -> ....


    annotate_cores_and_subtelomeric_contigs ${DATA_DIR}/out/23_annotate_cores_and_subt \
            ${DATA_DIR}/out/2_ref_reannotated_gaps/gaps.gff3 \
            ${DATA_DIR}/out/1_remove_gap_annotation/annotation.gapless.gff3 \
            ${DATA_DIR}/out/20_transfer_fixed_regions/annotation_combined.gff \
            "#\|closedgap_full\|closedgap_a\|closedgap_b\|gap" 
            # -> ${DATA_DIR}/out/23_annotate_cores_and_subt/annotation.gff

    extract_regions ${DATA_DIR}/out/24_extract_haploid_genome \
                ${DATA_DIR}/out/18.1_undo_failed_masking/masking_undone.fasta \
                ${DATA_DIR}/out/23_annotate_cores_and_subt/contig_and_subt.gff \
                ${DATA_DIR}/out/23_annotate_cores_and_subt/annotation.gff
                # -> ${DATA_DIR}/out/24_extract_haploid_genome/masked.fasta

    
    generate_overview_pic ${DATA_DIR}/out/25_overview_of_final_assembly \
                        ${DATA_DIR}/out/24_extract_haploid_genome/annotation.gff \
                        "gene filledgap closedgap_full closedgap_a closedgap_b expanded_region closedgap_masked gap" \
                        "gene=lightgrey;closedgap_full=green;closedgap_a=green;closedgap_b=green;closedgap_masked=green;expanded_region=blue;gap=purple"

    # here comes the analysis part !

    align_reads_to_genome ${DATA_DIR}/out/26_aligned_reads_on_new_genome \
        ${DATA_DIR}/out/18_closed_gaps/assembly.fasta \
        ${ONT_READS_IN}


    # analyze_gaps_closed_correctly ${DATA_DIR}/out/11_analyze_gaps_closed_correctly \
    #                               ${DATA_DIR}/out/3.1_gap_spanning_reads_old_genome/distance_deviation.tsv \
    #                               ${DATA_DIR}/out/10_vpr_new_genome/distance_deviation.tsv \
    #                               ${DATA_DIR}/out/3.1_gap_spanning_reads_old_genome/gap_spanning_reads.tsv \
    #                               ${DATA_DIR}/out/2_ref_reannotated_gaps/gaps.gff3 \
    #                               ${DATA_DIR}/out/1_move_centro_anno/annotation_combined.gff
    #                               # -> ${DATA_DIR}/out/11_analyze_gaps_closed_correctly/gaps_fixed.gff3

    # generate_overview_pic ${DATA_DIR}/out/12_overview_of_fixed_gaps \
    #                       ${DATA_DIR}/out/11_analyze_gaps_closed_correctly/gaps_fixed.gff3 \
    #                       "gap fixedgap" \
    #                       "gap=purple;fixedgap=green"

}

remove_annotated_gaps(){
    OUT_FOLDER=$1
    GFF_IN=$2

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/remove_annotated_gaps.done ]; then
        echo running remove_annotated_gaps in ${OUT_FOLDER}

        grep -vR "\tgap\t" ${GFF_IN} > ${OUT_FOLDER}/annotation.gapless.gff3

        echo "OK" > ${OUT_FOLDER}/remove_annotated_gaps.done
    fi
}

# Output:
# ${MaC_OUT_FOLDER}/7.1_undo_failed_masking/masking_undone.fasta
mask_and_close(){
    MaC_OUT_FOLDER=$1
    MaC_GENOME=$2
    MaC_READS_IN=$3
    MaC_ORIGINAL_GFF_IN=$4
    MaC_ORIGINAL_GENOME_IN=$5

    mkdir -p ${MaC_OUT_FOLDER}


    annotate_gaps ${MaC_OUT_FOLDER}/0_reannotated_gaps \
                ${MaC_GENOME}
                # -> ${MaC_OUT_FOLDER}/0_reannotated_gaps/gaps.gff3

    gap_spanning_reads ${MaC_OUT_FOLDER}/1_gap_spanning_reads \
                    ${MaC_GENOME} \
                    ${MaC_READS_IN} \
                    ${MaC_OUT_FOLDER}/0_reannotated_gaps/gaps.gff3
                    # -> ${MaC_OUT_FOLDER}/1_gap_spanning_reads/distance_deviation.tsv

    cut_superfluous_regions ${MaC_OUT_FOLDER}/2_cut_superfluous_regions \
                            ${MaC_GENOME} \
                            ${MaC_READS_IN} \
                            ${MaC_OUT_FOLDER}/1_gap_spanning_reads/distance_deviation.tsv \
                            ${MaC_OUT_FOLDER}/0_reannotated_gaps/gaps.gff3 \
                            # -> ${MaC_OUT_FOLDER}/2_cut_superfluous_regions/all_annotation.gff
                            # -> ${MaC_OUT_FOLDER}/2_cut_superfluous_regions/filtered_superfluous_regions.gff

    generate_overview_pic ${MaC_OUT_FOLDER}/3_overview_of_superfluous_regions \
                        ${MaC_OUT_FOLDER}/2_cut_superfluous_regions/all_annotation.gff \
                        "gene superfluous gap" \
                        "gene=lightgrey;gap=purple;superfluous=pink"

    mask_region ${MaC_OUT_FOLDER}/4_masked_superfluous_regions \
                ${MaC_GENOME} \
                ${MaC_OUT_FOLDER}/2_cut_superfluous_regions/filtered_superfluous_regions.gff
                # -> ${MaC_OUT_FOLDER}/4_masked_superfluous_regions/masked.fasta
                # -> ${MaC_OUT_FOLDER}/4_masked_superfluous_regions/removed_sequences.fasta

    annotate_gaps ${MaC_OUT_FOLDER}/5_annotate_new_gaps \
                ${MaC_OUT_FOLDER}/4_masked_superfluous_regions/masked.fasta
                # -> ${MaC_OUT_FOLDER}/5_annotate_new_gaps/gaps.gff3

    generate_overview_pic ${MaC_OUT_FOLDER}/6_overview_of_masked_regions \
                        ${MaC_OUT_FOLDER}/5_annotate_new_gaps/gaps.gff3 \
                        "gap" \
                        "gap=purple"

    close_gaps ${MaC_OUT_FOLDER}/7_closed_gaps \
            ${MaC_OUT_FOLDER}/4_masked_superfluous_regions \
            masked \
            ${MaC_READS_IN} \
            ${MaC_OUT_FOLDER}/5_annotate_new_gaps/gaps.gff3
            # -> ${MaC_OUT_FOLDER}/7_closed_gaps/gaps.gff3
            # -> ${MaC_OUT_FOLDER}/7_closed_gaps/assembly.fasta

    undo_masking ${MaC_OUT_FOLDER}/7.1_undo_failed_masking \
                ${MaC_OUT_FOLDER}/7_closed_gaps/assembly.fasta \
                ${MaC_OUT_FOLDER}/5_annotate_new_gaps/gaps.gff3 \
                ${MaC_OUT_FOLDER}/7_closed_gaps/gaps.gff3 \
                ${MaC_OUT_FOLDER}/4_masked_superfluous_regions/removed_sequences.fasta \
                ${MaC_OUT_FOLDER}/4_masked_superfluous_regions/masked.fasta
                # -> ${MaC_OUT_FOLDER}/7.1_undo_failed_masking/masking_undone.fasta
                # -> ${MaC_OUT_FOLDER}/7.1_undo_failed_masking/gaps.gff3

    transfer_annotation ${MaC_OUT_FOLDER}/8_transfer_annotation \
                        ${MaC_OUT_FOLDER}/7_closed_gaps/assembly.fasta \
                        ${MaC_ORIGINAL_GFF_IN} \
                        ${MaC_ORIGINAL_GENOME_IN} \
                        ${MaC_OUT_FOLDER}/7.1_undo_failed_masking/gaps.gff3
                        # -> ${MaC_OUT_FOLDER}/8_transfer_annotation/annotation_combined.gff

    generate_overview_pic ${MaC_OUT_FOLDER}/9_overview_of_remaining_gaps \
                        ${MaC_OUT_FOLDER}/8_transfer_annotation/annotation_combined.gff \
                        "gene filledgap gap" \
                        "gene=lightgrey;filledgap=green;gap=purple"

    # align_reads_to_genome ${MaC_OUT_FOLDER}/10_aligned_reads \
    #     ${MaC_OUT_FOLDER}/7_closed_gaps/assembly.fasta \
    #     ${MaC_READS_IN}



}

# Output:
# ${OUT_FOLDER}/reads.fasta
# ${OUT_FOLDER}/mates.fasta
#
# ${OUT_FOLDER}/reads.sam
# ${OUT_FOLDER}/mates.sam
#
# ${OUT_FOLDER}/reads.filtered.sam
# ${OUT_FOLDER}/mates.filtered.sam
#
# ${OUT_FOLDER}/distance_deviation.tsv
# ${OUT_FOLDER}/expected_distances.tsv
# ${OUT_FOLDER}/read_pos_and_strnd.tsv
#
virtual_paired_read_distance(){
    OUT_FOLDER=$1
    GENOME=$2
    READS_IN=$3

    mkdir -p ${OUT_FOLDER}

    VPR_LEN=500
    VPR_TRIALS=3
    MIN_ONT_LEN=$( expr "${VPR_LEN}" "*" "(" "${VPR_TRIALS}" "*" "2" "+" "2" ")" )

    if [ ! -e ${OUT_FOLDER}/virtual_paired_read_distance.done ]; then
        echo running virtual_paired_read_distance in ${OUT_FOLDER}

        zcat ${READS_IN} | python3 ${SCRIPTS_DIR}/illumina_from_ont.py - ${OUT_FOLDER}/reads.fasta ${OUT_FOLDER}/mates.fasta ${OUT_FOLDER}/expected_distances.tsv ${MIN_ONT_LEN} ${VPR_LEN} ${VPR_TRIALS}

        minimap2 -t 8 -ax map-ont ${GENOME} ${OUT_FOLDER}/reads.fasta > ${OUT_FOLDER}/reads.sam 2> ${OUT_FOLDER}/reads.minimap.errlog
        minimap2 -t 8 -ax map-ont ${GENOME} ${OUT_FOLDER}/mates.fasta > ${OUT_FOLDER}/mates.sam 2> ${OUT_FOLDER}/mates.minimap.errlog

        # # filter out low mapping quality reads
        samtools view -S -F 2304 ${OUT_FOLDER}/reads.sam > ${OUT_FOLDER}/reads.filtered.sam 
        samtools view -S -F 2304 ${OUT_FOLDER}/mates.sam > ${OUT_FOLDER}/mates.filtered.sam 

        # compute distances
        python3 ${SCRIPTS_DIR}/get_average_distance_deviation.py ${OUT_FOLDER}/reads.filtered.sam ${OUT_FOLDER}/mates.filtered.sam ${OUT_FOLDER}/expected_distances.tsv ${VPR_LEN} > ${OUT_FOLDER}/distance_deviation.tsv

        # compute distances
        python3 ${SCRIPTS_DIR}/get_read_pos_and_strand.py ${OUT_FOLDER}/reads.filtered.sam ${OUT_FOLDER}/mates.filtered.sam ${VPR_LEN} > ${OUT_FOLDER}/read_pos_and_strnd.tsv

        python3 ${SCRIPTS_DIR}/contig_connecting_reads.py ${OUT_FOLDER}/reads.sam ${OUT_FOLDER}/mates.sam ${VPR_LEN} > ${OUT_FOLDER}/contig_connecting_reads.tsv

        echo "OK" > ${OUT_FOLDER}/virtual_paired_read_distance.done
    fi
}


# Output:
# ${OUT_FOLDER}/gaps.gff3
annotate_gaps(){
    OUT_FOLDER=$1
    GENOME_IN=$2

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/annotate_gaps.done ]; then
        echo running annotate_gaps in ${OUT_FOLDER}

        python3 ${SCRIPTS_DIR}/annotate_gaps.py ${GENOME_IN} > ${OUT_FOLDER}/gaps.only.gff3

        echo "##gff-version 3" > ${OUT_FOLDER}/empty.annotation.gff3
        faidx ${GENOME_IN} -i chromsizes | awk '{OFS = "\t"; print "##sequence-region", $1, "1", $2}' >> ${OUT_FOLDER}/empty.annotation.gff3
        cat ${OUT_FOLDER}/empty.annotation.gff3 ${OUT_FOLDER}/gaps.only.gff3 > ${OUT_FOLDER}/gaps.gff3

        echo "OK" > ${OUT_FOLDER}/annotate_gaps.done
    fi
}

# Output:
# ${OUT_FOLDER}/assembly.fasta
# ${OUT_FOLDER}/gaps.gff3
# ${OUT_FOLDER}/stats.txt
#
close_gaps()
{
    OUT_FOLDER=$1
    REFERENCE_FOLDER=$2
    REFERENCE_NAME=$3
    READS_IN=$4
    GAPS_IN=$5

    mkdir -p ${OUT_FOLDER}


    if [ ! -e ${OUT_FOLDER}/close_gaps.done ]; then
        echo running close_gaps in ${OUT_FOLDER}

        cd ${OUT_FOLDER}
        ${BIN_DIR}/MaSuRCA-4.1.0/bin/close_scaffold_gaps.sh -r ${REFERENCE_FOLDER}/${REFERENCE_NAME}.fasta -q ${READS_IN} -t 18 -m 500 -i 95 -d ont
        cd ${SCRIPTS_DIR}

        python3 ${SCRIPTS_DIR}/fixup_number_of_n.py ${OUT_FOLDER}/${REFERENCE_NAME}.fasta.split.joined.fa 100 1000 > ${OUT_FOLDER}/assembly.fasta

        python3 ${SCRIPTS_DIR}/annotate_closed_gaps.py \
                ${REFERENCE_FOLDER}/${REFERENCE_NAME}.fasta \
                ${OUT_FOLDER}/assembly.fasta \
                ${GAPS_IN} \
            > ${OUT_FOLDER}/gaps.gff3


        echo "N's before running mascura (here, a gap consists of 1,000 Ns):" > ${OUT_FOLDER}/stats.txt
        grep -o "N" ${REFERENCE_FOLDER}/${REFERENCE_NAME}.fasta | wc -l >> ${OUT_FOLDER}/stats.txt
        echo "N's after running mascura (here, a gap consists of 1,000 Ns):" >> ${OUT_FOLDER}/stats.txt
        grep -o "N" ${OUT_FOLDER}/assembly.fasta | wc -l >> ${OUT_FOLDER}/stats.txt

        echo "Number of contigs before running samba:" >> ${OUT_FOLDER}/stats.txt
        grep -c ">" ${REFERENCE_FOLDER}/${REFERENCE_NAME}.fasta >> ${OUT_FOLDER}/stats.txt

        echo "Number of contigs after running samba:" >> ${OUT_FOLDER}/stats.txt
        grep -c ">" ${OUT_FOLDER}/assembly.fasta >> ${OUT_FOLDER}/stats.txt

        echo "OK" > ${OUT_FOLDER}/close_gaps.done
    fi
}

# Output:
# ${OUT_FOLDER}/assembly.fasta
#
merge_contigs()
{
    OUT_FOLDER=$1
    REFERENCE_FOLDER=$2
    REFERENCE_NAME=$3
    READS_IN=$4

    mkdir -p ${OUT_FOLDER}


    if [ ! -e ${OUT_FOLDER}/merge_contigs.done ]; then
        echo running merge_contigs in ${OUT_FOLDER}

        cd ${OUT_FOLDER}
        ${BIN_DIR}/MaSuRCA-4.1.0/bin/samba.sh -r ${REFERENCE_FOLDER}/${REFERENCE_NAME}.fasta -q ${READS_IN} -t 18 -m 500 -d ont
        cd ${SCRIPTS_DIR}

        cp ${OUT_FOLDER}/${REFERENCE_NAME}.fasta.scaffolds.fa ${OUT_FOLDER}/assembly.fasta 


        echo "Number of contigs before running samba:" >> ${OUT_FOLDER}/stats.txt
        grep -c ">" ${REFERENCE_FOLDER}/${REFERENCE_NAME}.fasta >> ${OUT_FOLDER}/stats.txt

        echo "Number of contigs after running samba:" >> ${OUT_FOLDER}/stats.txt
        grep -c ">" ${OUT_FOLDER}/assembly.fasta >> ${OUT_FOLDER}/stats.txt

        echo "OK" > ${OUT_FOLDER}/merge_contigs.done
    fi
}



# Output:
# ${OUT_FOLDER}/reference.sorted.fasta
# ${OUT_FOLDER}/query.sorted.fasta
#
# ${OUT_FOLDER}/comparison.paf
# ${OUT_FOLDER}/differences.paf
#
compare_assemblies(){
    OUT_FOLDER=$1
    REFERENCE=$2
    QUERY=$3

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/compare_assemblies.done ]; then
        echo running compare_assemblies in ${OUT_FOLDER}

        cat ${REFERENCE} | fassort > ${OUT_FOLDER}/reference.sorted.fasta

        cat ${QUERY} | fassort > ${OUT_FOLDER}/query.sorted.fasta

        minimap2 -x asm5 ${OUT_FOLDER}/query.sorted.fasta ${OUT_FOLDER}/reference.sorted.fasta > ${OUT_FOLDER}/comparison.paf

        python3 ${SCRIPTS_DIR}/differences_from_paf.py ${OUT_FOLDER}/comparison.paf > ${OUT_FOLDER}/differences.paf

        echo "OK" > ${OUT_FOLDER}/compare_assemblies.done
    fi 
}

# Output:
# ${OUT_FOLDER}/gap_spanning_reads.tsv
#
gap_spanning_reads(){
    OUT_FOLDER=$1
    GENOME=$2
    READS_IN=$3
    GAPS_IN=$4

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/gap_spanning_reads.done ]; then
        echo running gap_spanning_reads in ${OUT_FOLDER}
        
        virtual_paired_read_distance ${OUT_FOLDER} \
                                     ${GENOME} \
                                     ${READS_IN}
        # -> ${OUT_FOLDER}/reads.filtered.sam
        # -> ${OUT_FOLDER}/mates.filtered.sam

        python3 ${SCRIPTS_DIR}/spans_gap.py ${OUT_FOLDER}/reads.filtered.sam ${OUT_FOLDER}/mates.filtered.sam ${GAPS_IN} > ${OUT_FOLDER}/gap_spanning_reads.tsv

        echo "OK" > ${OUT_FOLDER}/gap_spanning_reads.done
    fi
}

# Output:
# ${OUT_FOLDER}/${ANNOTATION_NAME}.transfered.gff
# ${OUT_FOLDER}/annotation_combined.gff
move_annotation(){
    OUT_FOLDER=$1
    ANNOTATION_NAME=$2
    GFF_IN=$3
    ASSEMBLY_IN=$4
    ASSEMBLY_TRANSFER=$5
    GFF_MERGE=$6

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/move_annotation.done ]; then
        echo running move_annotation in ${OUT_FOLDER}

        grep ${ANNOTATION_NAME} ${GFF_IN} | grep -v "_3B_Tb427v10" > ${OUT_FOLDER}/${ANNOTATION_NAME}.filtered.gff


        transfer_annotation ${OUT_FOLDER} \
                            ${ASSEMBLY_TRANSFER} \
                            ${OUT_FOLDER}/${ANNOTATION_NAME}.filtered.gff \
                            ${ASSEMBLY_IN}
                            # -> ${OUT_FOLDER}/annotation.transfered.gff
        mv ${OUT_FOLDER}/annotation.transfered.gff ${OUT_FOLDER}/${ANNOTATION_NAME}.transfered.gff

        grep -v ${ANNOTATION_NAME} ${GFF_MERGE} > ${OUT_FOLDER}/gff_in.${ANNOTATION_NAME}less.gff3

        cat ${OUT_FOLDER}/gff_in.${ANNOTATION_NAME}less.gff3 ${OUT_FOLDER}/${ANNOTATION_NAME}.transfered.gff \
                > ${OUT_FOLDER}/annotation_combined.gff


        echo "OK" > ${OUT_FOLDER}/move_annotation.done
    fi

}


# Output:
# ${OUT_FOLDER}/annotation.transfered.gff
# ${OUT_FOLDER}/annotation_combined.gff
transfer_annotation(){
    OUT_FOLDER=$1
    ASSEMBLY_TRANSFER=$2
    GFF_IN=$3
    ASSEMBLY_IN=$4
    GFF_TRANSFER=$5

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/transfer_annotation.done ]; then
        echo running transfer_annotation in ${OUT_FOLDER}

        python3 ${SCRIPTS_DIR}/transfer_annotation_exact_match.py \
                ${ASSEMBLY_IN} \
                ${ASSEMBLY_TRANSFER} \
                ${GFF_IN} \
                ${OUT_FOLDER}/annotation.failed.gff \
            > ${OUT_FOLDER}/annotation.transfered.gff

        echo "failed to transfer this many annotations:"
        wc -l ${OUT_FOLDER}/annotation.failed.gff
        echo ""

        if [ ! -z "${GFF_TRANSFER}" ]
        then
            grep -v "#" ${GFF_TRANSFER} > ${OUT_FOLDER}/gff_transfer.no#.gff
            cat ${OUT_FOLDER}/annotation.transfered.gff ${OUT_FOLDER}/gff_transfer.no#.gff > ${OUT_FOLDER}/annotation_combined.gff
        fi
 
        echo "OK" > ${OUT_FOLDER}/transfer_annotation.done
    fi
}

# Output
# ${OUT_FOLDER}/overview.svg
generate_overview_pic(){
    OUT_FOLDER=$1
    GFF_IN=$2
    FEATURES=$3
    FEATURE_COLORS=$4

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/generate_overview_pic.done ]; then
        echo running generate_overview_pic in ${OUT_FOLDER}

        CONTIGS_WITH_GAPS=$(grep "##sequence-region" ${GFF_IN} | grep "Chr\|BES" | awk '{print $2}' | sort | uniq)
        # echo CONTIGS_WITH_GAPS: ${CONTIGS_WITH_GAPS}

        awk '{print $1, $2, $3, $4, $5, $6, $7, $8}' OFS='\t' ${GFF_IN} > ${OUT_FOLDER}/gff_last_column_removed.gff

        conda deactivate
        conda activate GENEastics_env
        python3 ${BIN_DIR}/geneastics.py \
            --replicons ${CONTIGS_WITH_GAPS} \
            --gff_file ${OUT_FOLDER}/gff_last_column_removed.gff \
            --feature_types ${FEATURES} \
            --alpha 0.99 \
            --feature_color_mapping ${FEATURE_COLORS} \
            --x_tick_distance 500000 \
            --font_size 1 \
            --output_file ${OUT_FOLDER}/overview.svg
        conda deactivate
        conda activate ont_assembly

        echo "" > ${OUT_FOLDER}/stats.txt
        for F in ${FEATURES}; do
            echo -n "${F} " >> ${OUT_FOLDER}/stats.txt
            grep -P "\t${F}\t" ${GFF_IN} | wc -l >> ${OUT_FOLDER}/stats.txt
        done

        echo "OK" > ${OUT_FOLDER}/generate_overview_pic.done
    fi
}

# Output:
# ${OUT_FOLDER}/gaps_closed_correctly.gff
# ${OUT_FOLDER}/gaps_fixed.gff3
analyze_gaps_closed_correctly(){
    OUT_FOLDER=$1
    DIST_DEV_OLD_GENOME=$2
    DIST_DEV_NEW_GENOME=$3
    GAP_SPANNING_READS=$4
    GAPS_OLD_GENOME=$5
    GFF_IN=$6

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/analyze_gaps_closed_correctly.done ]; then
        echo running analyze_gaps_closed_correctly in ${OUT_FOLDER}

        python3 ${SCRIPTS_DIR}/analyze_gaps_closed_correctly.py \
                ${DIST_DEV_OLD_GENOME} \
                ${DIST_DEV_NEW_GENOME} \
                ${GAP_SPANNING_READS} \
                ${GAPS_OLD_GENOME} \
            > ${OUT_FOLDER}/gaps_closed_correctly.gff

        grep -v "gap" ${GFF_IN} > ${OUT_FOLDER}/annotation.gapless.gff

        cat ${OUT_FOLDER}/annotation.gapless.gff ${OUT_FOLDER}/gaps_closed_correctly.gff > ${OUT_FOLDER}/gaps_fixed.gff3


        echo "OK" > ${OUT_FOLDER}/analyze_gaps_closed_correctly.done
    fi
}

# Output:
# ${OUT_FOLDER}/collapsed_regions.gff
# ${OUT_FOLDER}/annotation.gff
identify_collapsed_regions(){
    OUT_FOLDER=$1
    DIST_DEV_FILE=$2
    GAP_PATTERN_IN=$3
    GFF_IN=$4

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/identify_collapsed_regions.done ]; then
        echo running identify_collapsed_regions in ${OUT_FOLDER}

        grep ${GAP_PATTERN_IN} ${GFF_IN} > ${OUT_FOLDER}/gaps.gff

        python3 ${SCRIPTS_DIR}/identify_collapsed_regions.py \
                ${DIST_DEV_FILE} \
                ${OUT_FOLDER}/gaps.gff \
                ${OUT_FOLDER}/collapsed_regions.gff


        cat ${GFF_IN} ${GAPS_IN} ${OUT_FOLDER}/collapsed_regions.gff > ${OUT_FOLDER}/annotation.gff


        echo "OK" > ${OUT_FOLDER}/identify_collapsed_regions.done
    fi
}


# Output:
# ${OUT_FOLDER}/A.fasta
# ${OUT_FOLDER}/B.fasta
# ${OUT_FOLDER}/reads.A.fasta
# ${OUT_FOLDER}/reads.B.fasta
split_genome_in_a_and_b(){
    OUT_FOLDER=$1
    GENOME_IN=$2
    READS_IN=$3

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/split_genome_in_a_and_b.done ]; then
        echo running split_genome_in_a_and_b in ${OUT_FOLDER}

        # split genomes
        grep ">" ${GENOME_IN} | grep "A_Tb427v10" | awk '{print substr($1,2)}' > ${OUT_FOLDER}/A.lst
        grep ">" ${GENOME_IN} | grep "B_Tb427v10" | awk '{print substr($1,2)}' > ${OUT_FOLDER}/B.lst
        grep ">" ${GENOME_IN} | grep -v "A_Tb427v10" | grep -v "B_Tb427v10" | awk '{print substr($1,2)}' \
                > ${OUT_FOLDER}/remainder.lst

        ${BIN_DIR}/seqtk/seqtk subseq ${GENOME_IN} ${OUT_FOLDER}/A.lst > ${OUT_FOLDER}/A.fasta
        ${BIN_DIR}/seqtk/seqtk subseq ${GENOME_IN} ${OUT_FOLDER}/B.lst > ${OUT_FOLDER}/B.fasta
        ${BIN_DIR}/seqtk/seqtk subseq ${GENOME_IN} ${OUT_FOLDER}/remainder.lst > ${OUT_FOLDER}/remainder.fasta

        # split reads
        minimap2 -t 8 -ax map-ont ${GENOME} ${READS_IN} 2> ${OUT_FOLDER}/reads.minimap.errlog \
                | samtools view -q 30 -S -F 2304 > ${OUT_FOLDER}/reads.sam

        zcat ${READS_IN} | python3 remove_reads_mapping_to_contigs.py - ${OUT_FOLDER}/reads.sam ${OUT_FOLDER}/B.lst \
                    | python3 remove_reads_mapping_to_contigs.py - ${OUT_FOLDER}/reads.sam ${OUT_FOLDER}/remainder.lst \
                    > ${OUT_FOLDER}/reads.A.fasta
        zcat ${READS_IN} | python3 remove_reads_mapping_to_contigs.py - ${OUT_FOLDER}/reads.sam ${OUT_FOLDER}/A.lst \
                    | python3 remove_reads_mapping_to_contigs.py - ${OUT_FOLDER}/reads.sam ${OUT_FOLDER}/remainder.lst \
                    > ${OUT_FOLDER}/reads.B.fasta

        echo "OK" > ${OUT_FOLDER}/split_genome_in_a_and_b.done
    fi
}

# Output:
# ${OUT_FOLDER}/assembly.fasta
# ${OUT_FOLDER}/annotation.gff
merge_genomes(){
    OUT_FOLDER=$1
    GENOME_A=$2
    GENOME_B=$3
    GENOME_C=$4

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/merge_genomes.done ]; then
        echo running merge_genomes in ${OUT_FOLDER}

        cat ${GENOME_A} ${GENOME_B} ${GENOME_C} > ${OUT_FOLDER}/assembly.fasta

        annotate_gaps ${OUT_FOLDER} \
                      ${OUT_FOLDER}/assembly.fasta
                      # -> ${OUT_FOLDER}/annotation.gff

        # grep -v "#" ${GFF_A} > ${OUT_FOLDER}/a.no#.gff
        # grep -v "#" ${GFF_B} > ${OUT_FOLDER}/b.no#.gff

        # echo "##gff-version 3" > ${OUT_FOLDER}/annotation.gff
        # faidx ${OUT_FOLDER}/assembly.fasta -i chromsizes | awk '{OFS = "\t"; print "##sequence-region", $1, "1", $2}' >> ${OUT_FOLDER}/annotation.gff

        # cat ${OUT_FOLDER}/a.no#.gff ${OUT_FOLDER}/b.no#.gff >> ${OUT_FOLDER}/annotation.gff

        echo "OK" > ${OUT_FOLDER}/merge_genomes.done
    fi

}

cut_superfluous_regions(){
    OUT_FOLDER=$1
    GENOME_IN=$2
    ONT_IN=$3
    DIST_DEV_FILE=$4
    GAPS_IN=$5

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/cut_superfluous_regions.done ]; then
        echo running cut_superfluous_regions in ${OUT_FOLDER}

        # this is to inspect the regions that are cut out in jbrowse

        # but actually we use VPR and the gap annotation to remove the superfluous regions
        python3 ${SCRIPTS_DIR}/identify_collapsed_regions.py \
                ${DIST_DEV_FILE} \
                ${GAPS_IN} \
                ${OUT_FOLDER}/all_superfluous_regions.gff \
                inf \
                500 \
                False \
                False \
                superfluous

        python3 ${SCRIPTS_DIR}/identify_collapsed_regions.py \
                ${DIST_DEV_FILE} \
                ${GAPS_IN} \
                ${OUT_FOLDER}/filtered_superfluous_regions.gff \
                inf \
                500 \
                False \
                True \
                superfluous

        cat ${GAPS_IN} ${OUT_FOLDER}/all_superfluous_regions.gff > ${OUT_FOLDER}/all_annotation.gff
        cat ${GAPS_IN} ${OUT_FOLDER}/filtered_superfluous_regions.gff > ${OUT_FOLDER}/filtered_annotation.gff

        echo "OK" > ${OUT_FOLDER}/cut_superfluous_regions.done
    fi
}

align_reads_to_genome(){
    OUT_FOLDER=$1
    GENOME_IN=$2
    ONT_IN=$3

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/align_reads_to_genome.done ]; then
        echo running align_reads_to_genome in ${OUT_FOLDER}

        # this is to inspect the regions that are cut out in jbrowse

        minimap2 -t 8 -ax map-ont ${GENOME_IN} ${ONT_IN} 2> ${OUT_FOLDER}/reads.minimap.errlog | samtools view -hF 256 -q 30 | samtools sort -@ 8 -o ${OUT_FOLDER}/reads.bam
        samtools index -@ 8 ${OUT_FOLDER}/reads.bam


        echo "OK" > ${OUT_FOLDER}/align_reads_to_genome.done
    fi

}

# Output:
# ${OUT_FOLDER}/masked.fasta
# ${OUT_FOLDER}/removed_sequences.fasta
mask_region(){
    OUT_FOLDER=$1
    GENOME_IN=$2
    GFF_IN=$3

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/mask_region.done ]; then
        echo running mask_region in ${OUT_FOLDER}

        
        python3 ${SCRIPTS_DIR}/mask_regions.py \
            ${GENOME_IN} \
            ${GFF_IN} \
            ${OUT_FOLDER}/removed_sequences.fasta \
            > ${OUT_FOLDER}/masked.fasta

        echo "OK" > ${OUT_FOLDER}/mask_region.done
    fi
}

# Output:
# ${OUT_FOLDER}/new_assembly.fasta
extract_regions(){
    OUT_FOLDER=$1
    GENOME_IN=$2
    GFF_IN=$3
    ANNO_IN=$4

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/extract_regions.done ]; then
        echo running extract_regions in ${OUT_FOLDER}

        python3 ${SCRIPTS_DIR}/extract_region.py \
            ${GENOME_IN} \
            ${GFF_IN} \
            ${ANNO_IN} \
            ${OUT_FOLDER}/annotation.gff \
            > ${OUT_FOLDER}/new_assembly.fasta

        echo "OK" > ${OUT_FOLDER}/extract_regions.done
    fi
}

# Output:
# ${OUT_FOLDER}/gaps.gff
annotate_closed_gaps(){
    OUT_FOLDER=$1
    GENOME_IN=$2
    GENOME_AFTER_CLOSING=$3
    GAPS_BEFORE_CLOSING=$4

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/annotate_closed_gaps.done ]; then
        echo running annotate_closed_gaps in ${OUT_FOLDER}

        python3 ${SCRIPTS_DIR}/annotate_closed_gaps.py \
                ${GENOME_IN} \
                ${GENOME_AFTER_CLOSING} \
                ${GAPS_BEFORE_CLOSING} \
            > ${OUT_FOLDER}/gaps.gff3

        echo "OK" > ${OUT_FOLDER}/annotate_closed_gaps.done
    fi
}

# Output:
# ${OUT_FOLDER}/masking_undone.fasta
# ${OUT_FOLDER}/gaps.gff
undo_masking(){
    OUT_FOLDER=$1
    GENOME_AFTER_CLOSING=$2
    GAPS_BEFORE_CLOSING=$3
    GAPS_AFTER_CLOSING=$4
    MASKED_SEQUENCES=$5
    GENOME_BEFORE_CLOSING=$6

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/undo_masking.done ]; then
        echo running undo_masking in ${OUT_FOLDER}

        python3 ${SCRIPTS_DIR}/paste_old_sequence_into_remaining_gaps.py \
            ${GENOME_AFTER_CLOSING} \
            ${GAPS_BEFORE_CLOSING} \
            ${GAPS_AFTER_CLOSING} \
            ${MASKED_SEQUENCES} \
            > ${OUT_FOLDER}/masking_undone.fasta

        annotate_closed_gaps ${OUT_FOLDER} \
                ${GENOME_BEFORE_CLOSING} \
                ${OUT_FOLDER}/masking_undone.fasta \
                ${GAPS_BEFORE_CLOSING} \
            # -> ${OUT_FOLDER}/gaps.gff3


        python3 ${SCRIPTS_DIR}/distinguish_masked_from_unmasked_closed_gaps.py \
            ${GAPS_BEFORE_CLOSING} \
            ${OUT_FOLDER}/gaps.gff3 \
            ${MASKED_SEQUENCES} \
            > ${OUT_FOLDER}/fixed_gaps_and_masked.gff3

        echo "OK" > ${OUT_FOLDER}/undo_masking.done
    fi
}

# Output:
# ${OUT_FOLDER}/annotation_combined.gff
# ${OUT_FOLDER}/combined.transfered.gff
transfer_fixed_regions(){
    OUT_FOLDER=$1
    GFF_FILES_IN=$2
    BASE_ANNOTATION_IN=$3
    ASSEMBLY_TRANSFER=$4

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/transfer_fixed_regions.done ]; then
        echo running transfer_fixed_regions in ${OUT_FOLDER}

        cp ${BASE_ANNOTATION_IN} ${OUT_FOLDER}/annotation_combined.gff
        grep "#" ${BASE_ANNOTATION_IN} > ${OUT_FOLDER}/combined.transfered.gff

        while read line
        do
            arrIN=(${line//;/ })
            FEATURE_NAME_IN=${arrIN[0]}
            FEATURE_NAME_OUT=${arrIN[1]}
            GFF_FILE=${arrIN[2]}/fixed_gaps_and_masked.gff3
            ASSEMBLY_IN=${arrIN[2]}/masking_undone.fasta

            grep -P "${FEATURE_NAME_IN}" ${GFF_FILE} | awk -v feature_name="$FEATURE_NAME_OUT" '{print $1, $2, feature_name, $4, $5, $6, $7, $8, $9}' OFS='\t' > ${OUT_FOLDER}/${FEATURE_NAME_OUT}.gff
            
            python3 ${SCRIPTS_DIR}/transfer_annotation_exact_match.py \
                    ${ASSEMBLY_IN} \
                    ${ASSEMBLY_TRANSFER} \
                    ${OUT_FOLDER}/${FEATURE_NAME_OUT}.gff \
                    ${OUT_FOLDER}/${FEATURE_NAME_OUT}.failed.gff \
                > ${OUT_FOLDER}/${FEATURE_NAME_OUT}.transfered.gff

            echo "failed to transfer this many ${FEATURE_NAME_OUT}:"
            wc -l ${OUT_FOLDER}/${FEATURE_NAME_OUT}.failed.gff
            echo ""

            grep -v "#" ${OUT_FOLDER}/${FEATURE_NAME_OUT}.transfered.gff >> ${OUT_FOLDER}/annotation_combined.gff
            grep -v "#" ${OUT_FOLDER}/${FEATURE_NAME_OUT}.transfered.gff >> ${OUT_FOLDER}/combined.transfered.gff

        done < <(echo ${GFF_FILES_IN} | tr ' ' '\n')

        echo "OK" > ${OUT_FOLDER}/transfer_fixed_regions.done
    fi
}

# Output:
# -> ${OUT_FOLDER}/annotation.gff
annotate_cores_and_subtelomeric_contigs(){
    OUT_FOLDER=$1
    OLD_GAPS=$2
    CONTIG_GFF_IN=$3
    GFF_IN=$4
    NEW_GAP_PATTERN=$5

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/annotate_cores_and_subtelomeric_contigs.done ]; then
        echo running annotate_cores_and_subtelomeric_contigs in ${OUT_FOLDER}


        grep ${NEW_GAP_PATTERN} ${GFF_IN} > ${OUT_FOLDER}/new_gaps.gff

        grep -P "\tcontig\t" ${CONTIG_GFF_IN} | grep "core\|3A_\|3B_\|5A_\|5B_" | \
        python3 ${SCRIPTS_DIR}/place_contig_annotations_based_on_gaps.py \
                ${OLD_GAPS} \
                ${OUT_FOLDER}/new_gaps.gff \
                - \
            > ${OUT_FOLDER}/contig_and_subt.gff

        cp ${GFF_IN} ${OUT_FOLDER}/annotation.gff
        cat ${OUT_FOLDER}/contig_and_subt.gff >> ${OUT_FOLDER}/annotation.gff

        echo "OK" > ${OUT_FOLDER}/annotate_cores_and_subtelomeric_contigs.done
    fi
}


assembly_gaps_individually(){
    OUT_FOLDER=$1
    GAP_SPANNING_READS=$2
    ONT_READS_IN=$3

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/assembly_gaps_individually.done ]; then
        echo running assembly_gaps_individually in ${OUT_FOLDER}

        awk '{print $2}' ${GAP_SPANNING_READS} | sort | uniq > ${OUT_FOLDER}/gaps

        for GAP in $(cat ${OUT_FOLDER}/gaps)
        do
            mkdir -p ${OUT_FOLDER}/${GAP}
            awk 'NF==2{print}{}' ${GAP_SPANNING_READS} | grep ${GAP} | awk '{print substr($1, 1, length($1) - 2)}' \
                > ${OUT_FOLDER}/${GAP}/reads.lst

            if [ -s ${OUT_FOLDER}/${GAP}/reads.lst ]; then
                ${BIN_DIR}/seqtk/seqtk subseq ${ONT_READS_IN} ${OUT_FOLDER}/${GAP}/reads.lst > ${OUT_FOLDER}/${GAP}/reads.fasta

                cd ${OUT_FOLDER}/${GAP}

                ${BIN_DIR}/smartdenovo/smartdenovo.pl -c 1 ${OUT_FOLDER}/${GAP}/reads.fasta > wtasm.mak
                make -f wtasm.mak

                cd ${SCRIPTS_DIR}
            fi

        done

        echo "OK" > ${OUT_FOLDER}/assembly_gaps_individually.done
    fi
}


main


