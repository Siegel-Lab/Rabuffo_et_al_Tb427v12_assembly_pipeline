#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH --mem=500G
#SBATCH --time=7-00:00:00
#SBATCH --partition=fat
#SBATCH --job-name=samba

source /home/mschmidt/.miniconda3/etc/profile.d/conda.sh

conda activate ont_assembly

set -e
## STEP 1: configure


GENOME_FASTA_IN=$(realpath ../data/in/genome_in/HGAP3_Tb427v10_diploid/HGAP3_Tb427v10_diploid_scaffolded.fasta)
GFF_IN=$(realpath ../data/in/genome_in/HGAP3_Tb427v10_diploid/HGAP3_Tb427v10_diploid_scaffolded.gff3)
REF_CENTRO=$(realpath ../data/in/genome_in/HGAP3_Tb427v10/HGAP3_Tb427v10.fasta)
GFF_CENTRO_IN=$(realpath ../data/in/genome_in/HGAP3_Tb427v10/HGAP3_Tb427v10_manual.gff3)
ONT_READS_IN=$(realpath ../data/in/ont_reads_in/duplex_reads.fastq.gz)
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
mask_repeats(){
    REFERENCE=$1
    REF_NAME=$2
    REPEATS_IN=$3

    ## first find repeats
    # align the nanopore data onto the reference
    # call peaks in the coverage

    if [ ! -e ${MASK_REPEATS_DIR}/ont_reads_anligntment.sorted.bam ]; then
        minimap2 -a -x map-ont ${REFERENCE} ${ONT_READS_IN} | samtools view -bS - > ${MASK_REPEATS_DIR}/ont_reads_anligntment.bam
    fi 
    samtools sort -@ 18 ${MASK_REPEATS_DIR}/ont_reads_anligntment.bam > ${MASK_REPEATS_DIR}/ont_reads_anligntment.sorted.bam
    samtools index -@ 18 ${MASK_REPEATS_DIR}/ont_reads_anligntment.sorted.bam

    if [ ! -e ${MASK_REPEATS_DIR}/ont_reads_anligntment.bw ]; then
        bamCoverage --bam ${MASK_REPEATS_DIR}/ont_reads_anligntment.sorted.bam -o ${MASK_REPEATS_DIR}/ont_reads_anligntment.bw
    fi 

    if [ ! -e ${MASK_REPEATS_DIR}/ont_reads_anligntment.bed ]; then
        bedtools bamtobed -i ${MASK_REPEATS_DIR}/ont_reads_anligntment.sorted.bam > ${MASK_REPEATS_DIR}/ont_reads_anligntment.bed
    fi
    #awk '{OFS="\t"; print $1, $2, $3, $5}' ${MASK_REPEATS_DIR}/ont_reads_anligntment.bed > ${MASK_REPEATS_DIR}/ont_reads_anligntment.cut.bed
    faidx ${REFERENCE} -i chromsizes > ${MASK_REPEATS_DIR}/${REF_NAME}.sizes
    bedtools makewindows -g ${MASK_REPEATS_DIR}/${REF_NAME}.sizes -w 10000 -s 5000 > ${MASK_REPEATS_DIR}/${REF_NAME}.windows
    bedtools map -a ${MASK_REPEATS_DIR}/${REF_NAME}.windows -b ${MASK_REPEATS_DIR}/ont_reads_anligntment.bed -c 4 -o mean > ${MASK_REPEATS_DIR}/ont_reads_coverage_of_windows.bed

    # if [ ! -e ${MASK_REPEATS_DIR}/${REF_NAME}.masked.fasta ]; then
    #     python3 ${SCRIPTS_DIR}/mask_regions.py ${REFERENCE} ${REPEATS_IN} > ${MASK_REPEATS_DIR}/${REF_NAME}.fasta
    # fi
}

# mask_repeats ${GENOME_FASTA_IN} "reference" ${DATA_DIR}/in/mask_repeats/manual_repeats.gff


close_gaps()
{
    REFERENCE=$1
    REF_NAME=$2

    cd ${DATA_DIR}/out/samba_out_1

    ${MASUCRA_BIN}/close_scaffold_gaps.sh -r ${REFERENCE} -q ${ONT_READS_IN} -t 18 -m 500 -i 95 -d ont

    # .valid_join_pairs.txt contains how scaffolds have been split into contigs
    # .fasta.split.joined.fa is the main outputfile

    # .patches.coords seems interesting

    FILLED_N_ASSEMBLY=${DATA_DIR}/out/samba_out_1/HGAP3_Tb427v10_diploid_scaffolded.fasta.split.joined.fa
    FIXED_N_ASSEMBLY=${DATA_DIR}/out/samba_out_1/${REF_NAME}.fixed_n.fasta

    if [ ! -e ${FIXED_N_ASSEMBLY} ];then
        python3 ${SCRIPTS_DIR}/fixup_number_of_n.py ${FILLED_N_ASSEMBLY} 100 1000 > ${FIXED_N_ASSEMBLY}

    fi
    if [ ! -e ${REFERENCE}.gaps.gff3 ];then
        python3 ${SCRIPTS_DIR}/annotate_gaps.py ${REFERENCE} > ${REFERENCE}.gaps.gff3
    fi
    if [ ! -e ${FIXED_N_ASSEMBLY}.gaps.gff3 ];then
        python3 ${SCRIPTS_DIR}/annotate_gaps.py ${FIXED_N_ASSEMBLY} > ${FIXED_N_ASSEMBLY}.gaps.gff3
    fi

    # count the number of Ns in the input genome
    echo "N's in input assembly (a gap consists of 1,000 Ns):"
    grep -o "N" ${GENOME_FASTA_IN} | wc -l

    # count the number of Ns in the output genome
    echo "N's in fixed_n assembly (here, a gap consists of 1,000 Ns):"
    grep -o "N" ${FIXED_N_ASSEMBLY} | wc -l

    ############################################
    ## We have 49 gaps in the original assembly
    ## And 17 in the fixed one
    ############################################


    # count the number of contigs in the original assembly
    echo "Number of contigs in original assembly:"
    grep -c ">" ${GENOME_FASTA_IN}

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

close_gaps ${GENOME_FASTA_IN} reference


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

compare_assemblies ${GENOME_FASTA_IN} "reference" ${FIXED_N_ASSEMBLY} "fixed_n"




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

    echo "Distance deviation for ${NAME} is stored in file ${VIRT_PAIR_R_DIST}/${NAME}.distance_deviation"
}



virtual_paired_read_distance ${GENOME_FASTA_IN} "referece"


# check gap spanning
if [ ! -e ${VIRT_PAIR_R_DIST}/gap_spanning_reads ]; then
    python3 ${SCRIPTS_DIR}/spans_gap.py ${VIRT_PAIR_R_DIST}/referece.filtered.reads.sam ${VIRT_PAIR_R_DIST}/referece.filtered.mates.sam ${GENOME_FASTA_IN}.gaps.gff3 > ${VIRT_PAIR_R_DIST}/gap_spanning_reads
fi

virtual_paired_read_distance ${FIXED_N_ASSEMBLY} "fixed_n"


move_annotation(){
    ANNOTATION=$1

    conda deactivate
    conda activate ont_assembly_2
    
    module load ngs/bedtools2/2.26.0
    module load ncbi-blast/2.7.1+

    ## Get manual annotations from old genome version

    grep ${ANNOTATION} ${GFF_CENTRO_IN} > ${MOVE_ANNO_DIR}/${ANNOTATION}.filtered.gff
    # sed 's/core/A/g' ${MOVE_ANNO_DIR}/${ANNOTATION}.filtered.core.gff | sed 's/3A/A/g' - | sed 's/5A/A/g' - |  \
    #         grep -v "3B" | grep -v "5B" > ${MOVE_ANNO_DIR}/${ANNOTATION}.filtered.A.gff
    # sed 's/core/B/g' ${MOVE_ANNO_DIR}/${ANNOTATION}.filtered.core.gff | sed 's/3B/B/g' - | sed 's/5B/B/g' - |  \
    #         grep -v "3A" | grep -v "5A" > ${MOVE_ANNO_DIR}/${ANNOTATION}.filtered.B.gff
    # cat ${MOVE_ANNO_DIR}/${ANNOTATION}.filtered.A.gff ${MOVE_ANNO_DIR}/${ANNOTATION}.filtered.B.gff > ${MOVE_ANNO_DIR}/${ANNOTATION}.filtered.gff

    ## Crop coordinates on the entries based on chromosome length
    # Get chromosome lengths

    python3 $BIN_DIR/seq_length.py ${REF_CENTRO} > ${MOVE_ANNO_DIR}/genome.sizes

    # Crop coordinates based on chromosome length

    bedtools slop -i ${MOVE_ANNO_DIR}/${ANNOTATION}.filtered.gff -g ${MOVE_ANNO_DIR}/genome.sizes -b 0 > ${MOVE_ANNO_DIR}/${ANNOTATION}.cropped.gff

    ## Get fasta sequences for each entry and add it to the annotation file

    bedtools getfasta -tab -fi ${REF_CENTRO} -bed ${MOVE_ANNO_DIR}/${ANNOTATION}.cropped.gff > ${MOVE_ANNO_DIR}/${ANNOTATION}.sequences.fasta

    paste ${MOVE_ANNO_DIR}/${ANNOTATION}.cropped.gff ${MOVE_ANNO_DIR}/${ANNOTATION}.sequences.fasta > ${MOVE_ANNO_DIR}/${ANNOTATION}.sequence_annotation.out

    ## Run a modified version of Konrad's Förstner script to transfer annotation

    python3 $BIN_DIR/map_annotation_via_string_match.py ${FIXED_N_ASSEMBLY} ${MOVE_ANNO_DIR}/${ANNOTATION}.sequence_annotation.out ${MOVE_ANNO_DIR} ${MOVE_ANNO_DIR}/${ANNOTATION}.transfered.gff

    conda deactivate
    conda activate ont_assembly
}

# move_annotation "Centromere"


generate_overview_pic(){
    GFF_GAPLESS=${OVERVIEW_DIR}/reference.gapless.gff3
    GFF_FIXED_GAP=${OVERVIEW_DIR}/reference.gap_annotation_fixed.gff3


    if [ ! -e ${OVERVIEW_DIR}/annotation.gaps_closed.gff3 ]; then
        grep -v "gap" ${GFF_IN} > ${GFF_GAPLESS}
        cat ${GFF_GAPLESS} ${GENOME_FASTA_IN}.gaps.gff3 ${MOVE_ANNO_DIR}/Centromere.transfered.curated.gff > ${GFF_FIXED_GAP}

        python3 ${SCRIPTS_DIR}/close_gap_annotation_in_gff.py ${GFF_FIXED_GAP} ${FIXED_N_ASSEMBLY}.gaps.gff3 > ${OVERVIEW_DIR}/annotation.gaps_closed.gff3
    fi
    
    CONTIGS_WITH_GAPS=$(grep "gap" ${GFF_FIXED_GAP} | awk '{print $1} END { print "Chr2_B_Tb427v10" }' | sort | uniq)

    if [ ! -e ${OVERVIEW_DIR}/genome_overview_gaps.svg ]; then
        conda deactivate
        conda activate GENEastics_env
        python3 ${BIN_DIR}/geneastics.py \
            --replicons ${CONTIGS_WITH_GAPS} \
            --gff_file ${OVERVIEW_DIR}/annotation.gaps_closed.gff3 \
            --feature_types "gene" "Centromere" "gap" "closedgap" \
            --alpha 0.99 \
            --feature_color_mapping "Centromere=blue;gap=purple;gene=lightgrey;closedgap=purple" \
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
        if [ -e ${VIRT_PAIR_R_DIST}/closed_gaps_analysis.gff ]; then
            if [ ! -e ${OVERVIEW_DIR}/annotation.gaps_fixed.gff3 ]; then
                python3 ${SCRIPTS_DIR}/close_gap_annotation_in_gff.py ${GFF_FIXED_GAP} ${FIXED_N_ASSEMBLY}.gaps.gff3 ${VIRT_PAIR_R_DIST}/closed_gaps_analysis.gff > ${OVERVIEW_DIR}/annotation.gaps_fixed.gff3
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






## hmm this does not work as i want
# extend_repeats_long_stitch(){
#     conda deactivate
#     conda activate longstitch

#     TMP_READS=${LONG_STITCH_OUT}/reads
#     GENOME_SIZES=${LONG_STITCH_OUT}/genome.sizes
#     TMP_GNEOME=${LONG_STITCH_OUT}/genome

#     if [ ! -e ${TMP_READS}.fa.gz ]; then
#         zcat ${ONT_READS_IN} | sed -n '1~4s/^@/>/p;2~4p' | gzip > ${TMP_READS}.fa.gz
#     fi

#     if [ ! -e ${GENOME_SIZES} ]; then
#         faidx ${FIXED_N_ASSEMBLY} -i chromsizes > ${GENOME_SIZES}
#     fi

#     if [ ! -e ${TMP_GNEOME}.fa ]; then
#         cp ${FIXED_N_ASSEMBLY} ${TMP_GNEOME}.fa
#     fi

#     GENOME_SIZE=$(awk '{sum+=$2;} END{print sum;}' ${GENOME_SIZES})
#     echo "Genome size of fixed_n is ${GENOME_SIZE}"

#     cd ${LONG_STITCH_OUT}
#     longstitch run draft=genome reads=reads span=10 k_ntLink=24 w=100

#     REPEAT_EXTENDED_GENOME=$(realpath ./genome.k24.w100.tigmint-ntLink.longstitch-scaffolds.fa)

#     conda deactivate
#     conda activate ont_assembly
# }

#extend_repeats_long_stitch






