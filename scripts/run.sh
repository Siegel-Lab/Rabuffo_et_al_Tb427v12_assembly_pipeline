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
# - in the GFF files the first line of the gap annotations is pasted at the end of the last line of the remaining annotations (a newline is missing somewhere)
# - move all overview pictures to the new genome
#   - for this transfer annotations
# - do not filter out misassemblies that overlap existing gaps -> removing more sequence might be a good thing
#



main(){
    setup

    # move the Centromere annotation from the assembly where the cores are merged to the fully phased assembly
    # @todo this does not work with an exact match...
    move_annotation ${DATA_DIR}/out/1_move_centro_anno \
                    "Centromere" \
                    ${GFF_CENTRO_IN} \
                    ${REF_CENTRO} \
                    ${GENOME_FOLDER_IN}/${GENOME_FILENAME_IN}.fasta \
                    ${GFF_IN_XXX}
                    # -> ${DATA_DIR}/out/1_move_centro_anno/annotation_combined.gff

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
                   ${DATA_DIR}/out/1_move_centro_anno/annotation_combined.gff \
                    # -> ${DATA_DIR}/out/4_close_gaps_full_genome/7_closed_gaps/assembly.fasta
                    # -> ${DATA_DIR}/out/4_close_gaps_full_genome/7_closed_gaps/gaps.gff3
                    # -> ${DATA_DIR}/out/4_close_gaps_full_genome/8_transfer_annotation/annotation.transfered.gff

    split_genome_in_a_and_b ${DATA_DIR}/out/5_split_genome \
                            ${DATA_DIR}/out/4_close_gaps_full_genome/7_closed_gaps/assembly.fasta \
                            ${ONT_READS_IN}
                            # -> ${DATA_DIR}/out/5_split_genome/A.fasta
                            # -> ${DATA_DIR}/out/5_split_genome/B.fasta
                            # -> ${DATA_DIR}/out/5_split_genome/remainder.fasta
                            # -> ${DATA_DIR}/out/5_split_genome/reads.A.fasta
                            # -> ${DATA_DIR}/out/5_split_genome/reads.B.fasta

    
    mask_and_close ${DATA_DIR}/out/6_closed_gaps_a \
                   ${DATA_DIR}/out/5_split_genome/A.fasta \
                   ${DATA_DIR}/out/5_split_genome/reads.A.fasta \
                   ${DATA_DIR}/out/1_move_centro_anno/annotation_combined.gff \
                    # -> ${DATA_DIR}/out/6_closed_gaps_a/7_closed_gaps/assembly.fasta
                    # -> ${DATA_DIR}/out/6_closed_gaps_a/7_closed_gaps/gaps.gff3
                    # -> ${DATA_DIR}/out/6_closed_gaps_a/8_transfer_annotation/annotation.transfered.gff
    
    mask_and_close ${DATA_DIR}/out/7_closed_gaps_b \
                   ${DATA_DIR}/out/5_split_genome/B.fasta \
                   ${DATA_DIR}/out/5_split_genome/reads.B.fasta \
                   ${DATA_DIR}/out/1_move_centro_anno/annotation_combined.gff \
                    # -> ${DATA_DIR}/out/7_closed_gaps_b/7_closed_gaps/assembly.fasta
                    # -> ${DATA_DIR}/out/7_closed_gaps_b/7_closed_gaps/gaps.gff3
                    # -> ${DATA_DIR}/out/7_closed_gaps_b/8_transfer_annotation/annotation.transfered.gff


    merge_genomes ${DATA_DIR}/out/8_merged_genomes \
                  ${DATA_DIR}/out/6_closed_gaps_a/7_closed_gaps/assembly.fasta \
                  ${DATA_DIR}/out/6_closed_gaps_a/7_closed_gaps/gaps.gff3 \
                  ${DATA_DIR}/out/7_closed_gaps_b/7_closed_gaps/assembly.fasta \
                  ${DATA_DIR}/out/7_closed_gaps_b/7_closed_gaps/gaps.gff3 \
                  ${DATA_DIR}/out/5_split_genome/remainder.fasta
                 # -> ${DATA_DIR}/out/8_merged_genomes/assembly.fasta
                 # -> ${DATA_DIR}/out/8_merged_genomes/annotation.gff


    # transfer_annotation ${DATA_DIR}/out/8.1_transfer_annotation \
    #                     ${DATA_DIR}/out/8_merged_genomes/assembly.fasta \
    #                     ${DATA_DIR}/out/1_move_centro_anno/annotation_combined.gff \
    #                     ${REF_CENTRO}
    #                     # -> ${DATA_DIR}/out/8.1_transfer_annotation/annotation.transfered.gff

    generate_overview_pic ${DATA_DIR}/out/9_overview_of_remaining_gaps \
                           ${DATA_DIR}/out/8_merged_genomes/annotation.gff \
                          "gap" \
                          "gap=purple"
                        #   ${DATA_DIR}/out/8.1_transfer_annotation/annotation.transfered.gff \

    # virtual_paired_read_distance ${DATA_DIR}/out/10_vpr_new_genome \
    #                              ${DATA_DIR}/out/8_merged_genomes/assembly.fasta \
    #                              ${ONT_READS_IN}
    #                              # -> ${DATA_DIR}/out/10_vpr_new_genome/distance_deviation.tsv


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


    # identify_collapsed_regions ${DATA_DIR}/out/13_identify_collapsed_regions \
    #     ${DATA_DIR}/out/3.1_gap_spanning_reads_old_genome/distance_deviation.tsv \
    #     ${DATA_DIR}/out/4_closed_gaps/gaps.gff3 \
    #     ${DATA_DIR}/out/4_closed_gaps/empty.annotation.gff3
    #     # -> ${DATA_DIR}/out/13_identify_collapsed_regions/annotation.gff

    # generate_overview_pic ${DATA_DIR}/out/14_overview_collapsed_repeats \
    #                       ${DATA_DIR}/out/13_identify_collapsed_regions/annotation.gff \
    #                       "misassembly gap" \
    #                       "misassembly=pink;gap=purple"



}

# Output:
# ${MaC_OUT_FOLDER}/7_closed_gaps/assembly.fasta
# ${MaC_OUT_FOLDER}/8_transfer_annotation/annotation.transfered.gff
# ${MaC_OUT_FOLDER}/7_closed_gaps/gaps.gff3
mask_and_close(){
    MaC_OUT_FOLDER=$1
    MaC_GENOME=$2
    MaC_READS_IN=$3
    MaC_GFF_IN=$4

    mkdir -p ${MaC_OUT_FOLDER}
    if [ ! -e ${MaC_OUT_FOLDER}/mask_and_close.done ]; then
        echo running mask_and_close in ${MaC_OUT_FOLDER}

        
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
                                ${MaC_GFF_IN}
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
                ${MaC_READS_IN}
                # -> ${MaC_OUT_FOLDER}/7_closed_gaps/gaps.gff3
                # -> ${MaC_OUT_FOLDER}/7_closed_gaps/assembly.fasta

        transfer_annotation ${MaC_OUT_FOLDER}/8_transfer_annotation \
                            ${MaC_OUT_FOLDER}/7_closed_gaps/assembly.fasta \
                            ${MaC_GFF_IN} \
                            ${REF_CENTRO}
                            # -> ${MaC_OUT_FOLDER}/8_transfer_annotation/annotation.transfered.gff

        generate_overview_pic ${MaC_OUT_FOLDER}/9_overview_of_remaining_gaps \
                            ${MaC_OUT_FOLDER}/8_transfer_annotation/annotation.transfered.gff \
                            "gene filledgap gap" \
                            "gene=lightgrey;filledgap=green;gap=purple"


        echo "OK" > ${MaC_OUT_FOLDER}/mask_and_close.done
    fi
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

    mkdir -p ${OUT_FOLDER}


    if [ ! -e ${OUT_FOLDER}/close_gaps.done ]; then
        echo running close_gaps in ${OUT_FOLDER}

        cd ${OUT_FOLDER}
        ${BIN_DIR}/MaSuRCA-4.1.0/bin/close_scaffold_gaps.sh -r ${REFERENCE_FOLDER}/${REFERENCE_NAME}.fasta -q ${READS_IN} -t 18 -m 500 -i 95 -d ont
        cd ${SCRIPTS_DIR}

        python3 ${SCRIPTS_DIR}/fixup_number_of_n.py ${OUT_FOLDER}/${REFERENCE_NAME}.fasta.split.joined.fa 100 1000 > ${OUT_FOLDER}/assembly.fasta

        annotate_gaps ${OUT_FOLDER} \
                      ${OUT_FOLDER}/assembly.fasta
        # -> ${OUT_FOLDER}/gaps.gff3

        python3 ${SCRIPTS_DIR}/annotate_closed_gaps.py \
                ${REFERENCE_FOLDER}/${REFERENCE_NAME}.fasta \
                ${OUT_FOLDER}/assembly.fasta \
            >> ${OUT_FOLDER}/gaps.gff3


        echo "N's before running mascura (here, a gap consists of 1,000 Ns):" > ${OUT_FOLDER}/stats.txt
        grep -o "N" ${REFERENCE_FOLDER}/${REFERENCE_NAME}.fasta | wc -l >> ${OUT_FOLDER}/stats.txt
        echo "N's after running mascura (here, a gap consists of 1,000 Ns):" >> ${OUT_FOLDER}/stats.txt
        grep -o "N" ${OUT_FOLDER}/assembly.fasta | wc -l >> ${OUT_FOLDER}/stats.txt

        echo "Number of contigs before running mascura:" >> ${OUT_FOLDER}/stats.txt
        grep -c ">" ${REFERENCE_FOLDER}/${REFERENCE_NAME}.fasta >> ${OUT_FOLDER}/stats.txt

        echo "Number of contigs after running mascura:" >> ${OUT_FOLDER}/stats.txt
        grep -c ">" ${OUT_FOLDER}/assembly.fasta >> ${OUT_FOLDER}/stats.txt

        echo "OK" > ${OUT_FOLDER}/close_gaps.done
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
transfer_annotation(){
    OUT_FOLDER=$1
    ASSEMBLY_TRANSFER=$2
    GFF_IN=$3
    ASSEMBLY_IN=$4

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

        conda deactivate
        conda activate GENEastics_env
        python3 ${BIN_DIR}/geneastics.py \
            --replicons ${CONTIGS_WITH_GAPS} \
            --gff_file ${GFF_IN} \
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
            grep ${F} ${GFF_IN} | wc -l >> ${OUT_FOLDER}/stats.txt
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


        # python3 ${SCRIPTS_DIR}/close_gap_annotation_in_gff.py \
        #         ${GFF_IN} \
        #         ${OUT_FOLDER}/remaining.gaps.gff \
        #         ${OUT_FOLDER}/gaps_closed_correctly.gff \
        #     > ${OUT_FOLDER}/gaps_fixed.gff3


        echo "OK" > ${OUT_FOLDER}/analyze_gaps_closed_correctly.done
    fi
}

# Output:
# ${OUT_FOLDER}/collapsed_regions.gff
# ${OUT_FOLDER}/annotation.gff
identify_collapsed_regions(){
    OUT_FOLDER=$1
    DIST_DEV_FILE=$2
    GAPS_IN=$3
    GFF_IN=$4

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/identify_collapsed_regions.done ]; then
        echo running identify_collapsed_regions in ${OUT_FOLDER}


        python3 ${SCRIPTS_DIR}/identify_collapsed_regions.py \
                ${DIST_DEV_FILE} \
                ${GAPS_IN} \
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
    GFF_A=$3
    GENOME_B=$4
    GFF_B=$5
    GENOME_C=$6

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/merge_genomes.done ]; then
        echo running merge_genomes in ${OUT_FOLDER}

        cat ${GENOME_A} ${GENOME_B} ${GENOME_C} > ${OUT_FOLDER}/assembly.fasta

        grep -v "#" ${GFF_A} > ${OUT_FOLDER}/a.no#.gff
        grep "#" ${GFF_A} > ${OUT_FOLDER}/a.#.gff
        grep -v "#" ${GFF_B} > ${OUT_FOLDER}/b.no#.gff
        grep "#" ${GFF_B} | grep -v "gff-version" > ${OUT_FOLDER}/b.#.gff

        cat ${OUT_FOLDER}/a.#.gff ${OUT_FOLDER}/b.#.gff ${OUT_FOLDER}/a.no#.gff ${OUT_FOLDER}/b.no#.gff > ${OUT_FOLDER}/annotation.gff

        echo "OK" > ${OUT_FOLDER}/merge_genomes.done
    fi

}

cut_superfluous_regions(){
    OUT_FOLDER=$1
    GENOME_IN=$2
    ONT_IN=$3
    DIST_DEV_FILE=$4
    GAPS_IN=$5
    GFF_IN=$6

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/cut_superfluous_regions.done ]; then
        echo running cut_superfluous_regions in ${OUT_FOLDER}

        # this is to inspect the regions that are cut out in jbrowse

        # minimap2 -t 8 -ax map-ont ${GENOME_IN} ${ONT_IN} 2> ${OUT_FOLDER}/reads.minimap.errlog | samtools sort -o ${OUT_FOLDER}/reads.bam
        # samtools index ${OUT_FOLDER}/reads.bam
        # samtools view -bSF 256 -q 30 ${OUT_FOLDER}/reads.bam > ${OUT_FOLDER}/reads.filtered.bam
        # samtools index ${OUT_FOLDER}/reads.filtered.bam

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

        cat ${GFF_IN} ${GAPS_IN} ${OUT_FOLDER}/all_superfluous_regions.gff > ${OUT_FOLDER}/all_annotation.gff
        cat ${GFF_IN} ${GAPS_IN} ${OUT_FOLDER}/filtered_superfluous_regions.gff > ${OUT_FOLDER}/filtered_annotation.gff

        echo "OK" > ${OUT_FOLDER}/cut_superfluous_regions.done
    fi
}

# Output:
# ${OUT_FOLDER}/masked.fasta
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
            > ${OUT_FOLDER}/masked.fasta

        echo "OK" > ${OUT_FOLDER}/mask_region.done
    fi
}


main


