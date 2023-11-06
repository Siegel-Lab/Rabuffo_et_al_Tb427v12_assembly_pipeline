#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH --mem=500G
#SBATCH --time=7-00:00:00
#SBATCH --partition=fat
#SBATCH --job-name=ont_try_assembly
#SBATCH -o slurm_out/ont_assembly-%j.out

source /home/mschmidt/.miniconda3/etc/profile.d/conda.sh

conda activate ont_assembly

set -e
## STEP 1: configure



GENOME_FOLDER_IN=$(realpath ../data/in/genome_in/HGAP3_Tb427v10_diploid)
GENOME_FILENAME_IN="HGAP3_Tb427v10_diploid_scaffolded"
GFF_IN=$(realpath ../data/in/genome_in/HGAP3_Tb427v10_diploid/HGAP3_Tb427v10_diploid_scaffolded.gff3)
ANA_LYSIS_IN=$(realpath ../data/in/analysis_in)
REF_CENTRO=$(realpath ../data/in/genome_in/HGAP3_Tb427v10/HGAP3_Tb427v10.fasta)
GFF_CENTRO_IN=$(realpath ../data/in/genome_in/HGAP3_Tb427v10/HGAP3_Tb427v10_manual.gff3)
ONT_READS_IN=$(realpath ../data/in/ont_reads_in/merged.nanopore.gz)
MASUCRA_BIN=$(realpath ../bin/MaSuRCA-4.1.0/bin)
COMPARISON_DIR=../data/out/comparison
VIRT_PAIR_R_DIST=../data/out/virtual_paired_read_dist
MASK_REPEATS_DIR=../data/out/mask_repeats
SAMBA_OUT=../data/out/samba_out_1
MOVE_ANNO_DIR=../data/out/move_anno_dir

BIN_DIR=$(realpath ../bin/)
OVERVIEW_DIR=../data/out/overview
DATA_DIR=$(realpath ../data)
SCRIPTS_DIR=$(realpath .)


setup() {
    mkdir -p ${DATA_DIR}/out
    mkdir -p ${COMPARISON_DIR}
    COMPARISON_DIR=$(realpath ${COMPARISON_DIR})
    mkdir -p ${VIRT_PAIR_R_DIST}
    VIRT_PAIR_R_DIST=$(realpath ${VIRT_PAIR_R_DIST})
    mkdir -p ${MASK_REPEATS_DIR}
    MASK_REPEATS_DIR=$(realpath ${MASK_REPEATS_DIR})
    mkdir -p ${OVERVIEW_DIR}
    OVERVIEW_DIR=$(realpath ${OVERVIEW_DIR})
    mkdir -p ${SAMBA_OUT}
    SAMBA_OUT=$(realpath ${SAMBA_OUT})
    mkdir -p ${MOVE_ANNO_DIR}
    MOVE_ANNO_DIR=$(realpath ${MOVE_ANNO_DIR})
}

setup



module load ngs/minimap2/2.10
module load ngs/samtools/1.9
module load ngs/deeptools/3.5.0
module load ngs/bedtools2/2.28.0


virtual_paired_read_distance(){
    GENOME=$1
    NAME=$2
    # create virtual illumina reads from the ont reads
    # align them to both assemblies and check their expected to actual distance
    # this is a measure of how well the assembly is doing

    # lets create the virtual illumina reads

    if [ ! -e ${VIRT_PAIR_R_DIST}/${NAME}.reads.fasta ]; then
        zcat ${ONT_READS_IN} | python3 ${SCRIPTS_DIR}/illumina_from_ont.py - ${VIRT_PAIR_R_DIST}/${NAME}.reads.fasta ${VIRT_PAIR_R_DIST}/${NAME}.mates.fasta ${VIRT_PAIR_R_DIST}/${NAME}.expected_distances 2000 500
    fi 

    
    if [ ! -e ${VIRT_PAIR_R_DIST}/${NAME}.reads.sam ]; then
        minimap2 -ax map-ont ${GENOME} ${VIRT_PAIR_R_DIST}/${NAME}.reads.fasta > ${VIRT_PAIR_R_DIST}/${NAME}.reads.sam 2> ${VIRT_PAIR_R_DIST}/${NAME}.reads.minimap.errlog
        minimap2 -ax map-ont ${GENOME} ${VIRT_PAIR_R_DIST}/${NAME}.mates.fasta > ${VIRT_PAIR_R_DIST}/${NAME}.mates.sam 2> ${VIRT_PAIR_R_DIST}/${NAME}.mates.minimap.errlog
    fi 

    # filter out low mapping quality reads
    if [ ! -e ${VIRT_PAIR_R_DIST}/${NAME}.filtered.reads.sam ]; then
        samtools view -Sq 30 -F 2304 ${VIRT_PAIR_R_DIST}/${NAME}.reads.sam > ${VIRT_PAIR_R_DIST}/${NAME}.filtered.reads.sam 
        samtools view -Sq 30 -F 2304 ${VIRT_PAIR_R_DIST}/${NAME}.mates.sam > ${VIRT_PAIR_R_DIST}/${NAME}.filtered.mates.sam 
    fi 

    # compute distances
    if [ ! -e ${VIRT_PAIR_R_DIST}/${NAME}.distance_deviation ]; then
        python3 ${SCRIPTS_DIR}/get_average_distance_deviation.py ${VIRT_PAIR_R_DIST}/${NAME}.filtered.reads.sam ${VIRT_PAIR_R_DIST}/${NAME}.filtered.mates.sam ${VIRT_PAIR_R_DIST}/${NAME}.expected_distances > ${VIRT_PAIR_R_DIST}/${NAME}.distance_deviation
    fi

    # compute distances
    if [ ! -e ${VIRT_PAIR_R_DIST}/${NAME}.read_pos_and_strnd ]; then
        python3 ${SCRIPTS_DIR}/get_read_pos_and_strand.py ${VIRT_PAIR_R_DIST}/${NAME}.filtered.reads.sam ${VIRT_PAIR_R_DIST}/${NAME}.filtered.mates.sam > ${VIRT_PAIR_R_DIST}/${NAME}.read_pos_and_strnd
    fi

    echo "Distance deviation for ${NAME} is stored in file ${VIRT_PAIR_R_DIST}/${NAME}.distance_deviation"
}




close_gaps()
{
    REFERENCE_FOLDER=$1
    REFERENCE_NAME=$2
    REF_NAME=$3

    cd ${DATA_DIR}/out/samba_out_1

    ${MASUCRA_BIN}/close_scaffold_gaps.sh -r ${REFERENCE_FOLDER}/${REFERENCE_NAME}.fasta -q ${ONT_READS_IN} -t 18 -m 500 -i 95 -d ont

    # .valid_join_pairs.txt contains how scaffolds have been split into contigs
    # .fasta.split.joined.fa is the main outputfile

    # .patches.coords seems interesting

    FILLED_N_ASSEMBLY=${DATA_DIR}/out/samba_out_1/${REFERENCE_NAME}.masked.fasta.split.joined.fa
    FIXED_N_ASSEMBLY=${DATA_DIR}/out/samba_out_1/${REF_NAME}.fixed_n.fasta

    if [ ! -e ${FIXED_N_ASSEMBLY} ];then
        python3 ${SCRIPTS_DIR}/fixup_number_of_n.py ${FILLED_N_ASSEMBLY} 100 1000 > ${FIXED_N_ASSEMBLY}

    fi
    if [ ! -e ${DATA_DIR}/out/samba_out_1/${REF_NAME}.gaps.gff3 ];then
        python3 ${SCRIPTS_DIR}/annotate_gaps.py ${REFERENCE_FOLDER}/${REFERENCE_NAME}.fasta > ${DATA_DIR}/out/samba_out_1/${REF_NAME}.gaps.gff3
    fi
    if [ ! -e ${FIXED_N_ASSEMBLY}.gaps.gff3 ];then
        python3 ${SCRIPTS_DIR}/annotate_gaps.py ${FIXED_N_ASSEMBLY} > ${FIXED_N_ASSEMBLY}.gaps.gff3
    fi

    # count the number of Ns in the output genome
    echo "N's in fixed_n assembly (here, a gap consists of 1,000 Ns):"
    grep -o "N" ${FIXED_N_ASSEMBLY} | wc -l

    ############################################
    ## We have 49 gaps in the original assembly
    ## And 17 in the fixed one
    ############################################


    # count the number of contigs in the original assembly
    echo "Number of contigs in original assembly:"
    grep -c ">" ${REFERENCE}

    # count the number of contigs in the output assembly
    echo "Number of contigs in n_filled assembly:"
    grep -c ">" ${FIXED_N_ASSEMBLY}

    
    ############################################
    ## We have 317 contigs in the original assembly
    ## And 308 in the fixed one
    ## 
    ## In total 32 gaps have been closed?, 9 of them using a unitig present in the original assembly
    ##
    ############################################
}

close_gaps ${GENOME_FOLDER_IN} ${GENOME_FILENAME_IN} reference


compare_assemblies(){
    REFERENCE=$1
    REF_NAME=$2
    QUERY=$3
    Q_NAME=$4
    echo ${REFERENCE} ${REF_NAME} ${QUERY} ${Q_NAME}

    
    if [ ! -e ${DATA_DIR}/out/comparison/${REF_NAME}.sorted.fasta ];then
        cat ${REFERENCE} | fassort > ${DATA_DIR}/out/comparison/${REF_NAME}.sorted.fasta
    fi 
    if [ ! -e ${DATA_DIR}/comparison/${Q_NAME}.sorted.fasta ];then
        cat ${QUERY} | fassort > ${DATA_DIR}/out/comparison/${Q_NAME}.sorted.fasta
    fi 


    if [ ! -e ${COMPARISON_DIR}/${REF_NAME}.${Q_NAME}.comparison.paf ];then
        minimap2 -x asm5 ${DATA_DIR}/out/comparison/${Q_NAME}.sorted.fasta ${DATA_DIR}/out/comparison/${REF_NAME}.sorted.fasta > ${COMPARISON_DIR}/${REF_NAME}.${Q_NAME}.comparison.paf
    fi 

    if [ ! -e ${COMPARISON_DIR}/${REF_NAME}.${Q_NAME}.differences.paf ];then
        python3 ${SCRIPTS_DIR}/differences_from_paf.py ${COMPARISON_DIR}/${REF_NAME}.${Q_NAME}.comparison.paf > ${COMPARISON_DIR}/${REF_NAME}.${Q_NAME}.differences.paf
    fi 
}

compare_assemblies ${GENOME_FOLDER_IN}/${GENOME_FILENAME_IN}.fasta "reference" ${FIXED_N_ASSEMBLY} "fixed_n"

virtual_paired_read_distance ${GENOME_FOLDER_IN}/${GENOME_FILENAME_IN}.fasta "referece"


# check gap spanning
if [ ! -e ${VIRT_PAIR_R_DIST}/gap_spanning_reads ]; then
    python3 ${SCRIPTS_DIR}/spans_gap.py ${VIRT_PAIR_R_DIST}/referece.filtered.reads.sam ${VIRT_PAIR_R_DIST}/referece.filtered.mates.sam ${DATA_DIR}/out/samba_out_1/referece.gaps.gff3 > ${VIRT_PAIR_R_DIST}/gap_spanning_reads
fi

virtual_paired_read_distance ${FIXED_N_ASSEMBLY} "fixed_n"


move_annotation(){
    ANNOTATION=$1

    conda deactivate
    conda activate ont_assembly_2
    
    module load ngs/bedtools2/2.26.0
    module load ncbi-blast/2.7.1+

    ## Get manual annotations from old genome version

    if [ ! -e ${MOVE_ANNO_DIR}/${ANNOTATION}.filtered.gff ]; then
        grep ${ANNOTATION} ${GFF_CENTRO_IN} > ${MOVE_ANNO_DIR}/${ANNOTATION}.filtered.gff
    fi

    ## Crop coordinates on the entries based on chromosome length
    # Get chromosome lengths

    if [ ! -e ${MOVE_ANNO_DIR}/${ANNOTATION}.transfered.gff ]; then
        python3 $BIN_DIR/seq_length.py ${REF_CENTRO} > ${MOVE_ANNO_DIR}/genome.sizes

        # Crop coordinates based on chromosome length

        bedtools slop -i ${MOVE_ANNO_DIR}/${ANNOTATION}.filtered.gff -g ${MOVE_ANNO_DIR}/genome.sizes -b 0 > ${MOVE_ANNO_DIR}/${ANNOTATION}.cropped.gff

        ## Get fasta sequences for each entry and add it to the annotation file

        bedtools getfasta -tab -fi ${REF_CENTRO} -bed ${MOVE_ANNO_DIR}/${ANNOTATION}.cropped.gff > ${MOVE_ANNO_DIR}/${ANNOTATION}.sequences.fasta

        paste ${MOVE_ANNO_DIR}/${ANNOTATION}.cropped.gff ${MOVE_ANNO_DIR}/${ANNOTATION}.sequences.fasta > ${MOVE_ANNO_DIR}/${ANNOTATION}.sequence_annotation.out

        ## Run a modified version of Konrad's Förstner script to transfer annotation

        python3 $BIN_DIR/map_annotation_via_string_match.py ${FIXED_N_ASSEMBLY} ${MOVE_ANNO_DIR}/${ANNOTATION}.sequence_annotation.out ${MOVE_ANNO_DIR} ${MOVE_ANNO_DIR}/${ANNOTATION}.transfered.gff
    fi

    conda deactivate
    conda activate ont_assembly
}

move_annotation "Centromere"


generate_overview_pic(){
    GFF_GAPLESS=${OVERVIEW_DIR}/reference.gapless.gff3
    GFF_FIXED_GAP=${OVERVIEW_DIR}/reference.gap_annotation_fixed.gff3


    if [ ! -e ${OVERVIEW_DIR}/open_gaps.gff ]; then
        grep -v "gap" ${GFF_IN} > ${GFF_GAPLESS}
        cat ${GFF_GAPLESS} ${DATA_DIR}/out/samba_out_1/reference.gaps.gff3 ${MOVE_ANNO_DIR}/Centromere.transfered.curated.gff ${ANA_LYSIS_IN}/mis_assemblies.gff > ${OVERVIEW_DIR}/open_gaps.gff
    fi

    if [ ! -e ${OVERVIEW_DIR}/annotation.gaps_closed.gff3 ]; then
        grep -v "gap" ${GFF_IN} > ${GFF_GAPLESS}
        cat ${GFF_GAPLESS} ${DATA_DIR}/out/samba_out_1/reference.gaps.gff3 ${MOVE_ANNO_DIR}/Centromere.transfered.curated.gff ${ANA_LYSIS_IN}/mis_assemblies.gff > ${GFF_FIXED_GAP}

        python3 ${SCRIPTS_DIR}/close_gap_annotation_in_gff.py ${GFF_FIXED_GAP} ${FIXED_N_ASSEMBLY}.gaps.gff3 > ${OVERVIEW_DIR}/annotation.gaps_closed.gff3
    fi
    
    CONTIGS_WITH_GAPS=$(grep "gap" ${GFF_FIXED_GAP} | awk '{print $1} END { print "Chr2_B_Tb427v10" }' | sort | uniq)

    if [ ! -e ${OVERVIEW_DIR}/genome_overview_gaps.svg ]; then
        conda deactivate
        conda activate GENEastics_env
        python3 ${BIN_DIR}/geneastics.py \
            --replicons ${CONTIGS_WITH_GAPS} \
            --gff_file ${OVERVIEW_DIR}/open_gaps.gff \
            --feature_types "gene" "misassembly" "gap" "Centromere" \
            --alpha 0.99 \
            --feature_color_mapping "Centromere=blue;gap=purple;gene=lightgrey;misassembly=pink" \
            --x_tick_distance 500000 \
            --font_size 1 \
            --output_file ${OVERVIEW_DIR}/genome_overview_gaps.svg
        conda deactivate
    fi
    if [ ! -e ${OVERVIEW_DIR}/genome_overview_closed.svg ]; then
        conda deactivate
        conda activate GENEastics_env
        python3 ${BIN_DIR}/geneastics.py \
            --replicons ${CONTIGS_WITH_GAPS} \
            --gff_file ${OVERVIEW_DIR}/annotation.gaps_closed.gff3 \
            --feature_types "gene" "Centromere" "gap" "closedgap" \
            --alpha 0.99 \
            --feature_color_mapping "Centromere=blue;gap=purple;gene=lightgrey;closedgap=green" \
            --x_tick_distance 500000 \
            --font_size 1 \
            --output_file ${OVERVIEW_DIR}/genome_overview_closed.svg
        conda deactivate
    fi
    if [ ! -e ${OVERVIEW_DIR}/genome_overview_fixed.svg ]; then
        if [ -e ${ANA_LYSIS_IN}/closed_gaps_analysis.gff ]; then
            if [ ! -e ${OVERVIEW_DIR}/annotation.gaps_fixed.gff3 ]; then
                python3 ${SCRIPTS_DIR}/close_gap_annotation_in_gff.py ${GFF_FIXED_GAP} ${FIXED_N_ASSEMBLY}.gaps.gff3 ${ANA_LYSIS_IN}/closed_gaps_analysis.gff > ${OVERVIEW_DIR}/annotation.gaps_fixed.gff3
            fi
            conda deactivate
            conda activate GENEastics_env
            python3 ${BIN_DIR}/geneastics.py \
                --replicons ${CONTIGS_WITH_GAPS} \
                --gff_file ${OVERVIEW_DIR}/annotation.gaps_fixed.gff3 \
                --feature_types "gene" "Centromere" "gap" "closedgap" "fixedgap" "notenoughdatagap" \
                --alpha 0.99 \
                --feature_color_mapping "Centromere=blue;gap=purple;gene=lightgrey;closedgap=red;notenoughdatagap=orange;fixedgap=green" \
                --x_tick_distance 500000 \
                --font_size 1 \
                --output_file ${OVERVIEW_DIR}/genome_overview_fixed.svg
            conda deactivate
        fi
    fi
    conda activate ont_assembly
}


generate_overview_pic


mask_repeats(){
    REFERENCE=$1
    REF_NAME=$2

    if [ ! -e ${MASK_REPEATS_DIR}/${REF_NAME}.mis_assemblies.gff ];then
        python3 ${SCRIPTS_DIR}/identify_collapsed_regions.py ${VIRT_PAIR_R_DIST}/${NAME}.distance_deviation ${DATA_DIR}/out/samba_out_1/reference.gaps.gff3 ${MASK_REPEATS_DIR}/${REF_NAME}.mis_assemblies.gff
    fi

    if [ ! -e ${MASK_REPEATS_DIR}/${REF_NAME}.masked.fasta ];then
        python3 ${SCRIPTS_DIR}/mask_regions.py ${REFERENCE} ${MASK_REPEATS_DIR}/${REF_NAME}.mis_assemblies.gff > ${MASK_REPEATS_DIR}/${REF_NAME}.masked.fasta
    fi

    # count the number of Ns in the input genome
    echo "N's in input assembly (a gap consists of 1,000 Ns):"
    grep -o "N" ${REFERENCE} | wc -l

    # count the number of Ns in the output genome
    echo "N's in masked assembly (here, a gap consists of 1,000 Ns):"
    grep -o "N" ${MASK_REPEATS_DIR}/${REF_NAME}.masked.fasta | wc -l
}

exit
mask_repeats ${FIXED_N_ASSEMBLY} "fixed_n"







