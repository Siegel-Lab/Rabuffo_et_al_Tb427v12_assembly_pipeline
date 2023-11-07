#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH --mem=500G
#SBATCH --time=7-00:00:00
#SBATCH --partition=fat
#SBATCH --job-name=ont_try_assembly
#SBATCH -o slurm_out/ont_assembly-%j.out

## STEP 1: configure

#
# META:
# Every function represents one stage of the processing.
# The various output files of every function are stored under a common folder (should be a folder in data/out/).
# The inputs of a function are coolectedd from the previous outputs
#



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
    GFF_IN=$(realpath ../data/in/genome_in/HGAP3_Tb427v10_diploid/HGAP3_Tb427v10_diploid_scaffolded.gff3)
    ANA_LYSIS_IN=$(realpath ../data/in/analysis_in)
    REF_CENTRO=$(realpath ../data/in/genome_in/HGAP3_Tb427v10/HGAP3_Tb427v10.fasta)
    GFF_CENTRO_IN=$(realpath ../data/in/genome_in/HGAP3_Tb427v10/HGAP3_Tb427v10_manual.gff3)
    ONT_READS_IN=$(realpath ../data/in/ont_reads_in/merged.nanopore.gz)
    
    BIN_DIR=$(realpath ../bin/)
    DATA_DIR=$(realpath ../data)
    SCRIPTS_DIR=$(realpath .)


    mkdir -p ${DATA_DIR}/out
}




main(){
    setup


    move_annotation "${DATA_DIR}/out/1_move_centro_anno" \
                    "Centromere" \
                    ${GFF_CENTRO_IN} \
                    ${REF_CENTRO} \
                    ${GENOME_FOLDER_IN}/${GENOME_FILENAME_IN}.fasta \
                    ${GFF_IN}
    # -> "${DATA_DIR}/out/1_move_centro_anno/annotation_combined.gff"


    annotate_gaps "${DATA_DIR}/out/2_ref_reannotated_gaps" \
                  ${GENOME_FOLDER_IN}/${GENOME_FILENAME_IN}.fasta \
                  "${DATA_DIR}/out/1_move_centro_anno/annotation_combined.gff"
    # -> "${DATA_DIR}/out/2_ref_reannotated_gaps/reannotated.gff3"


    generate_overview_pic "${DATA_DIR}/out/3_overview_of_gaps" \
                          "${DATA_DIR}/out/2_ref_reannotated_gaps/reannotated.gff3" \
                          "gene gap Centromere" \
                          "Centromere=blue;gap=purple;gene=lightgrey" \
                            ${GENOME_FOLDER_IN}/${GENOME_FILENAME_IN}.fasta


    close_gaps "${DATA_DIR}/out/4_closed_gaps" \
               ${GENOME_FOLDER_IN} \
               ${GENOME_FILENAME_IN} \
               ${ONT_READS_IN}
    # -> "${DATA_DIR}/out/4_closed_gaps/closed_gaps.gff3"
    # -> "${DATA_DIR}/out/4_closed_gaps/assembly.fasta"

    generate_overview_pic "${DATA_DIR}/out/5_overview_of_closed_gaps" \
                          "${DATA_DIR}/out/4_closed_gaps/closed_gaps.gff3" \
                          "gene gap Centromere closedgap" \
                          "Centromere=blue;gap=purple;gene=lightgrey;closedgap=green" \
                          ${GENOME_FOLDER_IN}/${GENOME_FILENAME_IN}.fasta








    # compare_assemblies ${GENOME_FOLDER_IN}/${GENOME_FILENAME_IN}.fasta "reference" ${FIXED_N_ASSEMBLY} "fixed_n"

    # virtual_paired_read_distance ${GENOME_FOLDER_IN}/${GENOME_FILENAME_IN}.fasta "referece"

    # # gap_spanning_reads blub ${VIRT_PAIR_R_DIST}/referece.filtered.reads.sam ${VIRT_PAIR_R_DIST}/referece.filtered.mates.sam ${DATA_DIR}/out/samba_out_1/reference.gaps.gff3

    # # check gap spanning

    # virtual_paired_read_distance ${FIXED_N_ASSEMBLY} "fixed_n"

    # move_annotation "Centromere"

    # generate_overview_pic


    # mask_repeats(){
    #     REFERENCE=$1
    #     REF_NAME=$2

    #     if [ ! -e ${MASK_REPEATS_DIR}/${REF_NAME}.masked.fasta ];then
    #         python3 ${SCRIPTS_DIR}/mask_regions.py ${REFERENCE} ${MASK_REPEATS_DIR}/${REF_NAME}.mis_assemblies.gff > ${MASK_REPEATS_DIR}/${REF_NAME}.masked.fasta
    #     fi

    #     # count the number of Ns in the input genome
    #     echo "N's in input assembly (a gap consists of 1,000 Ns):"
    #     grep -o "N" ${REFERENCE} | wc -l

    #     # count the number of Ns in the output genome
    #     echo "N's in masked assembly (here, a gap consists of 1,000 Ns):"
    #     grep -o "N" ${MASK_REPEATS_DIR}/${REF_NAME}.masked.fasta | wc -l
    # }

    # mask_repeats ${FIXED_N_ASSEMBLY} "fixed_n"


    # close_gaps ${MASK_REPEATS_DIR} ${REF_NAME}.masked fixed_repeats
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

    if [ ! -e ${OUT_FOLDER}/virtual_paired_read_distance.done ]; then
        echo running virtual_paired_read_distance in ${OUT_FOLDER}

        zcat ${READS_IN} | python3 ${SCRIPTS_DIR}/illumina_from_ont.py - ${OUT_FOLDER}/reads.fasta ${OUT_FOLDER}/mates.fasta ${OUT_FOLDER}/expected_distances.tsv 2000 500

        minimap2 -ax map-ont ${GENOME} ${OUT_FOLDER}/reads.fasta > ${OUT_FOLDER}/reads.sam 2> ${OUT_FOLDER}/reads.minimap.errlog
        minimap2 -ax map-ont ${GENOME} ${OUT_FOLDER}/mates.fasta > ${OUT_FOLDER}/mates.sam 2> ${OUT_FOLDER}/mates.minimap.errlog

        # filter out low mapping quality reads
        samtools view -Sq 30 -F 2304 ${OUT_FOLDER}/reads.sam > ${OUT_FOLDER}/reads.filtered.sam 
        samtools view -Sq 30 -F 2304 ${OUT_FOLDER}/mates.sam > ${OUT_FOLDER}/mates.filtered.sam 

        # compute distances
        python3 ${SCRIPTS_DIR}/get_average_distance_deviation.py ${OUT_FOLDER}/reads.filtered.sam ${OUT_FOLDER}/mates.filtered.sam ${OUT_FOLDER}/expected_distances.tsv > ${OUT_FOLDER}/distance_deviation.tsv

        # compute distances
        python3 ${SCRIPTS_DIR}/get_read_pos_and_strand.py ${OUT_FOLDER}/reads.filtered.sam ${OUT_FOLDER}/mates.filtered.sam > ${OUT_FOLDER}/read_pos_and_strnd.tsv

        echo "OK" > ${OUT_FOLDER}/virtual_paired_read_distance.done
    fi
}


# Output:
# ${OUT_FOLDER}/gaps.gff3
# ${OUT_FOLDER}/reannotated.gff3
annotate_gaps(){
    OUT_FOLDER=$1
    GENOME_IN=$2
    GFF_IN=$3

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/annotate_gaps.done ]; then
        echo running annotate_gaps in ${OUT_FOLDER}

        python3 ${SCRIPTS_DIR}/annotate_gaps.py ${GENOME_IN} > ${OUT_FOLDER}/gaps.gff3

        grep -v "gap" ${GFF_IN} > ${OUT_FOLDER}/gff_in.gapless.gff3
        cat ${OUT_FOLDER}/gff_in.gapless.gff3 ${OUT_FOLDER}/gaps.gff3 > ${OUT_FOLDER}/reannotated.gff3

        echo "OK" > ${OUT_FOLDER}/annotate_gaps.done
    fi
}

# Output:
# ${OUT_FOLDER}/assembly.fasta
# ${OUT_FOLDER}/gaps.gff3
# ${OUT_FOLDER}/stats.txt
# ${OUT_FOLDER}/closed_gaps.gff3
#
close_gaps()
{
    OUT_FOLDER=$1
    REFERENCE_FOLDER=$2
    REFERENCE_NAME=$3
    READS_IN=$4
    GFF_IN=$4

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/close_gaps.done ]; then
        echo running close_gaps in ${OUT_FOLDER}

        cd ${OUT_FOLDER}
        ${BIN_DIR}/MaSuRCA-4.1.0/bin/close_scaffold_gaps.sh -r ${REFERENCE_FOLDER}/${REFERENCE_NAME}.fasta -q ${READS_IN} -t 18 -m 500 -i 95 -d ont
        cd ${SCRIPTS_DIR}

        python3 ${SCRIPTS_DIR}/fixup_number_of_n.py ${OUT_FOLDER}/${REFERENCE_NAME}.fasta.split.joined.fa 100 1000 > ${OUT_FOLDER}/assembly.fasta

        annotate_gaps ${OUT_FOLDER} ${OUT_FOLDER}/assembly.fasta

        python3 ${SCRIPTS_DIR}/close_gap_annotation_in_gff.py ${GFF_IN} gaps.gff3 > ${OUT_FOLDER}/closed_gaps.gff3

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
    READS_IN=$2
    MATES_IN=$3
    GAPS_IN=$4

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/gap_spanning_reads.done ]; then
        echo running gap_spanning_reads in ${OUT_FOLDER}

        python3 ${SCRIPTS_DIR}/spans_gap.py ${READS_IN} ${MATES_IN} ${GAPS_IN} > ${OUT_FOLDER}/gap_spanning_reads.tsv

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

        conda deactivate
        conda activate ont_assembly_2

        module load ngs/bedtools2/2.26.0
        module load ncbi-blast/2.7.1+

        grep ${ANNOTATION_NAME} ${GFF_IN} > ${OUT_FOLDER}/${ANNOTATION_NAME}.filtered.gff
        
        python3 ${BIN_DIR}/seq_length.py ${ASSEMBLY_IN} > ${OUT_FOLDER}/genome.sizes

        # Crop coordinates based on chromosome length

        bedtools slop -i ${OUT_FOLDER}/${ANNOTATION_NAME}.filtered.gff -g ${OUT_FOLDER}/genome.sizes -b 0 > ${OUT_FOLDER}/${ANNOTATION_NAME}.cropped.gff

        ## Get fasta sequences for each entry and add it to the annotation file

        bedtools getfasta -tab -fi ${ASSEMBLY_IN} -bed ${OUT_FOLDER}/${ANNOTATION_NAME}.cropped.gff > ${OUT_FOLDER}/${ANNOTATION_NAME}.sequences.fasta

        paste ${OUT_FOLDER}/${ANNOTATION_NAME}.cropped.gff ${OUT_FOLDER}/${ANNOTATION_NAME}.sequences.fasta > ${OUT_FOLDER}/${ANNOTATION_NAME}.sequence_annotation.out

        ## Run a modified version of Konrad's Förstner script to transfer annotation

        python3 ${BIN_DIR}/map_annotation_via_string_match.py ${ASSEMBLY_TRANSFER} ${OUT_FOLDER}/${ANNOTATION_NAME}.sequence_annotation.out ${OUT_FOLDER} ${OUT_FOLDER}/${ANNOTATION_NAME}.transfered.gff

        grep -v ${ANNOTATION_NAME} ${GFF_MERGE} > ${OUT_FOLDER}/gff_in.${ANNOTATION_NAME}less.gff3
        cat ${OUT_FOLDER}/gff_in.${ANNOTATION_NAME}less.gff3 ${OUT_FOLDER}/${ANNOTATION_NAME}.transfered.gff > ${OUT_FOLDER}/annotation_combined.gff

        conda deactivate
        conda activate ont_assembly

        echo "OK" > ${OUT_FOLDER}/move_annotation.done
    fi



}

# Output
# ${OUT_FOLDER}/overview.svg
generate_overview_pic(){
    OUT_FOLDER=$1
    GFF_IN=$2
    FEATURES=$3
    FEATURE_COLORS=$4
    GENOME_IN=$5

    mkdir -p ${OUT_FOLDER}

    if [ ! -e ${OUT_FOLDER}/generate_overview_pic.done ]; then
        echo running generate_overview_pic in ${OUT_FOLDER}

        CONTIGS_WITH_GAPS=$(grep "Chr\|BES" ${GFF_IN} | awk '{print $1}' | grep -v "#" | sort | uniq)
        #echo CONTIGS_WITH_GAPS: ${CONTIGS_WITH_GAPS}

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

        echo "OK" > ${OUT_FOLDER}/generate_overview_pic.done
    fi
}


# generate_overview_pic(){
#     GFF_GAPLESS=${OVERVIEW_DIR}/reference.gapless.gff3
#     GFF_FIXED_GAP=${OVERVIEW_DIR}/reference.gap_annotation_fixed.gff3


#     if [ ! -e ${MASK_REPEATS_DIR}/${REF_NAME}.mis_assemblies.gff ];then
#         python3 ${SCRIPTS_DIR}/identify_collapsed_regions.py ${VIRT_PAIR_R_DIST}/${NAME}.distance_deviation ${DATA_DIR}/out/samba_out_1/reference.gaps.gff3 ${MASK_REPEATS_DIR}/${REF_NAME}.mis_assemblies.gff
#     fi

#     if [ ! -e ${OVERVIEW_DIR}/open_gaps.gff ]; then
#         grep -v "gap" ${GFF_IN} > ${GFF_GAPLESS}
#         cat ${GFF_GAPLESS} ${DATA_DIR}/out/samba_out_1/reference.gaps.gff3 ${MOVE_ANNO_DIR}/Centromere.transfered.curated.gff ${MASK_REPEATS_DIR}/${REF_NAME}.mis_assemblies.gff > ${OVERVIEW_DIR}/open_gaps.gff
#     fi

#     if [ ! -e ${OVERVIEW_DIR}/annotation.gaps_closed.gff3 ]; then
#         grep -v "gap" ${GFF_IN} > ${GFF_GAPLESS}
#         cat ${GFF_GAPLESS} ${DATA_DIR}/out/samba_out_1/reference.gaps.gff3 ${MOVE_ANNO_DIR}/Centromere.transfered.curated.gff ${MASK_REPEATS_DIR}/${REF_NAME}.mis_assemblies.gff > ${GFF_FIXED_GAP}

#         python3 ${SCRIPTS_DIR}/close_gap_annotation_in_gff.py ${GFF_FIXED_GAP} ${FIXED_N_ASSEMBLY}.gaps.gff3 > ${OVERVIEW_DIR}/annotation.gaps_closed.gff3
#     fi
    
#     CONTIGS_WITH_GAPS=$(grep "gap" ${GFF_FIXED_GAP} | awk '{print $1} END { print "Chr2_B_Tb427v10" }' | sort | uniq)

#     if [ ! -e ${OVERVIEW_DIR}/genome_overview_gaps.svg ]; then
#         conda deactivate
#         conda activate GENEastics_env
#         python3 ${BIN_DIR}/geneastics.py \
#             --replicons ${CONTIGS_WITH_GAPS} \
#             --gff_file ${OVERVIEW_DIR}/open_gaps.gff \
#             --feature_types "gene" "misassembly" "gap" "Centromere" \
#             --alpha 0.99 \
#             --feature_color_mapping "Centromere=blue;gap=purple;gene=lightgrey;misassembly=pink" \
#             --x_tick_distance 500000 \
#             --font_size 1 \
#             --output_file ${OVERVIEW_DIR}/genome_overview_gaps.svg
#         conda deactivate
#     fi
#     if [ ! -e ${OVERVIEW_DIR}/genome_overview_closed.svg ]; then
#         conda deactivate
#         conda activate GENEastics_env
#         python3 ${BIN_DIR}/geneastics.py \
#             --replicons ${CONTIGS_WITH_GAPS} \
#             --gff_file ${OVERVIEW_DIR}/annotation.gaps_closed.gff3 \
#             --feature_types "gene" "Centromere" "gap" "closedgap" \
#             --alpha 0.99 \
#             --feature_color_mapping "Centromere=blue;gap=purple;gene=lightgrey;closedgap=green" \
#             --x_tick_distance 500000 \
#             --font_size 1 \
#             --output_file ${OVERVIEW_DIR}/genome_overview_closed.svg
#         conda deactivate
#     fi
#     if [ ! -e ${OVERVIEW_DIR}/genome_overview_fixed.svg ]; then
#         if [ -e ${ANA_LYSIS_IN}/closed_gaps_analysis.gff ]; then
#             if [ ! -e ${OVERVIEW_DIR}/annotation.gaps_fixed.gff3 ]; then
#                 python3 ${SCRIPTS_DIR}/close_gap_annotation_in_gff.py ${GFF_FIXED_GAP} ${FIXED_N_ASSEMBLY}.gaps.gff3 ${ANA_LYSIS_IN}/closed_gaps_analysis.gff > ${OVERVIEW_DIR}/annotation.gaps_fixed.gff3
#             fi
#             conda deactivate
#             conda activate GENEastics_env
#             python3 ${BIN_DIR}/geneastics.py \
#                 --replicons ${CONTIGS_WITH_GAPS} \
#                 --gff_file ${OVERVIEW_DIR}/annotation.gaps_fixed.gff3 \
#                 --feature_types "gene" "Centromere" "gap" "closedgap" "fixedgap" "notenoughdatagap" \
#                 --alpha 0.99 \
#                 --feature_color_mapping "Centromere=blue;gap=purple;gene=lightgrey;closedgap=red;notenoughdatagap=orange;fixedgap=green" \
#                 --x_tick_distance 500000 \
#                 --font_size 1 \
#                 --output_file ${OVERVIEW_DIR}/genome_overview_fixed.svg
#             conda deactivate
#         fi
#     fi
#     conda activate ont_assembly
# }


main


