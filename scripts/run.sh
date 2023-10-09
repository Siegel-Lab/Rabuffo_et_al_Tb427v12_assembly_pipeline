#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH --mem=250G
#SBATCH --time=7-00:00:00
#SBATCH --partition=slim18
#SBATCH --job-name=samba
#SBATCH -o ../data/logfiles/slurm-out-samba-%j.out

source /home/mschmidt/.miniconda3/etc/profile.d/conda.sh

conda activate ont_assembly

set -e
## STEP 1: configure


GENOME_FASTA_IN=$(realpath ../data/genome_in/HGAP3_Tb427v10/HGAP3_Tb427v10.fasta)
GFF_IN=$(realpath ../data/genome_in/HGAP3_Tb427v10/HGAP3_Tb427v10_manual.gff3)
ONT_READS_IN=$(realpath ../data/ont_reads_in/duplex_reads.fastq.gz)
MASUCRA_BIN=$(realpath ../bin/MaSuRCA-4.1.0/bin)
MINIM_OUT_FILE=$(realpath ../data/comparison/genome_alignment)
DIFFERENCES_OUT_FILE=$(realpath ../data/comparison/genome_differences)
VIRT_PAIR_R_DIST=$(realpath ../data/virtual_paired_read_dist)

BIN_DIR=$(realpath ../bin/)
OVERVIEW_DIR=$(realpath ../data/overview)
DATA_DIR=$(realpath ../data)
SCRIPTS_DIR=$(realpath .)

## STEP 2: run samba to fill 'N' gaps in the assembly


fixup_assembly()
{
    cd ${DATA_DIR}/samba_out_1

    ${MASUCRA_BIN}/close_scaffold_gaps.sh -r ${GENOME_FASTA_IN} -q ${ONT_READS_IN} -t 18 -m 2000 -d ont

    # .valid_join_pairs.txt contains how scaffolds have been split into contigs
    # .fasta.split.joined.fa is the main outputfile

    # .patches.coords seems interesting

    FILLED_N_ASSEMBLY=${DATA_DIR}/samba_out_1/HGAP3_Tb427v10.fasta.split.joined.fa
    FIXED_N_ASSEMBLY=${DATA_DIR}/samba_out_1/HGAP3_Tb427v10.fixed_n.fasta

    if [ ! -e ${FIXED_N_ASSEMBLY} ];then
        python3 ${SCRIPTS_DIR}/fixup_number_of_n.py ${FILLED_N_ASSEMBLY} 100 1000 > ${FIXED_N_ASSEMBLY}
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

fixup_assembly




module load ngs/minimap2/2.10
module load ngs/samtools/1.9
virtual_paired_read_distance(){
    GENOME=$1
    NAME=$2
    # create virtual illumina reads from the ont reads
    # align them to both assemblies and check their expected to actual distance
    # this is a measure of how well the assembly is doing

    # lets create the virtual illumina reads

    if [ ! -e ${VIRT_PAIR_R_DIST}/${NAME}.reads.fasta ];then
        zcat ${ONT_READS_IN} | python3 ${SCRIPTS_DIR}/illumina_from_ont.py - ${VIRT_PAIR_R_DIST}/${NAME}.reads.fasta ${VIRT_PAIR_R_DIST}/${NAME}.mates.fasta ${VIRT_PAIR_R_DIST}/${NAME}.expected_distances 10000 2000
    fi 

    
    if [ ! -e ${VIRT_PAIR_R_DIST}/${NAME}.reads.sam ];then
        minimap2 -ax map-ont ${GENOME} ${VIRT_PAIR_R_DIST}/${NAME}.reads.fasta > ${VIRT_PAIR_R_DIST}/${NAME}.reads.sam
        minimap2 -ax map-ont ${GENOME} ${VIRT_PAIR_R_DIST}/${NAME}.mates.fasta > ${VIRT_PAIR_R_DIST}/${NAME}.mates.sam
    fi 

    # filter out low mapping quality reads
    if [ ! -e ${VIRT_PAIR_R_DIST}/${NAME}.filtered.reads.sam ];then
        samtools view -Sq 30 ${VIRT_PAIR_R_DIST}/${NAME}.reads.sam > ${VIRT_PAIR_R_DIST}/${NAME}.filtered.reads.sam 
        samtools view -Sq 30 ${VIRT_PAIR_R_DIST}/${NAME}.mates.sam > ${VIRT_PAIR_R_DIST}/${NAME}.filtered.mates.sam 
    fi 

    # compute distances
    if [ ! -e ${VIRT_PAIR_R_DIST}/${NAME}.distance_deviation ];then
        python3 ${SCRIPTS_DIR}/get_average_distance_deviation.py ${VIRT_PAIR_R_DIST}/${NAME}.filtered.reads.sam ${VIRT_PAIR_R_DIST}/${NAME}.filtered.mates.sam ${VIRT_PAIR_R_DIST}/${NAME}.expected_distances > ${VIRT_PAIR_R_DIST}/${NAME}.distance_deviation
    fi

    echo "Average distance deviation for ${NAME}:"
    cat ${VIRT_PAIR_R_DIST}/${NAME}.distance_deviation

}



# virtual_paired_read_distance ${GENOME_FASTA_IN} "referece"
# virtual_paired_read_distance ${FIXED_N_ASSEMBLY} "fixed_n"




generate_overview_pic(){    
    # "BES10_Tb427v10" "BES11_Tb427v10" "BES12_Tb427v10" "BES13_Tb427v10" "BES14_Tb427v10" "BES15_Tb427v10" "BES17_Tb427v10" "BES1_Tb427v10" "BES2_Tb427v10" "BES3_Tb427v10" "BES4_Tb427v10" "BES5_Tb427v10" "BES7_Tb427v10" "Chr10_3A_Tb427v10" "Chr10_3B_Tb427v10" "Chr10_5A_Tb427v10" "Chr10_5B_Tb427v10" "Chr10_core_Tb427v10" "Chr11_3A_Tb427v10" "Chr11_3B_Tb427v10" "Chr11_5A_Tb427v10" "Chr11_5B_Tb427v10" "Chr11_core_Tb427v10" "Chr1_3A_Tb427v10" "Chr1_3B_Tb427v10" "Chr1_5A_Tb427v10" "Chr1_5B_Tb427v10" "Chr1_core_Tb427v10" "Chr2_5A_Tb427v10" "Chr2_core_Tb427v10" "Chr3_3A_Tb427v10" "Chr3_5A_Tb427v10" "Chr3_core_Tb427v10" "Chr4_3A_Tb427v10" "Chr4_3B_Tb427v10" "Chr4_5A_Tb427v10" "Chr4_5B_Tb427v10" "Chr4_core_Tb427v10" "Chr5_3A_Tb427v10" "Chr5_3B_Tb427v10" "Chr5_core_Tb427v10" "Chr6_3A_Tb427v10" "Chr6_3B_Tb427v10" "Chr6_core_Tb427v10" "Chr7_5A_Tb427v10" "Chr7_core_Tb427v10" "Chr8_3A_Tb427v10" "Chr8_3B_Tb427v10" "Chr8_5A_Tb427v10" "Chr8_5B_Tb427v10" "Chr8_core_Tb427v10" "Chr9_3A_Tb427v10" "Chr9_3B_Tb427v10" "Chr9_5A_Tb427v10" "Chr9_5B_Tb427v10" "Chr9_core_Tb427v10" "Chr3_5B_Tb427v10" \

    # "BES17_Tb427v10" "BES2_Tb427v10" "Chr1_3A_Tb427v10" "Chr1_3B_Tb427v10" "Chr1_core_Tb427v10" "Chr3_5A_Tb427v10" "Chr3_core_Tb427v10" "Chr4_core_Tb427v10" "Chr5_3A_Tb427v10" "Chr5_3B_Tb427v10" "Chr5_core_Tb427v10" "Chr6_3A_Tb427v10" "Chr6_3B_Tb427v10" "Chr7_core_Tb427v10" "Chr8_3A_Tb427v10" "Chr8_5A_Tb427v10" "Chr8_5B_Tb427v10" "Chr8_core_Tb427v10" "Chr9_3A_Tb427v10" "Chr9_3B_Tb427v10" "Chr9_5A_Tb427v10" "Chr9_core_Tb427v10" "Chr10_3A_Tb427v10" "Chr11_3A_Tb427v10" "Chr11_3B_Tb427v10" "Chr11_core_Tb427v10"

    CONTIGS_WITH_GAPS=$(grep "gap" ${GFF_IN} | awk '{print $1}' | sort | uniq)

    python3 ${SCRIPTS_DIR}/close_gap_annotation_in_gff.py ${GFF_IN} ${FIXED_N_ASSEMBLY}.gaps.gff3 > ${OVERVIEW_DIR}/annotation.gaps_closed.gff3

    conda activate GENEastics_env
    python3 ${BIN_DIR}/geneastics.py \
        --replicons ${CONTIGS_WITH_GAPS} \
        --gff_file ${OVERVIEW_DIR}/annotation.gaps_closed.gff3 \
        --feature_types "pseudogene" "Centromere" "gap" "closedgap" \
        --alpha 0.99 \
        --feature_color_mapping "Centromere=blue;gap=red;pseudogene=lightgrey;closedgap=green" \
        --attribute_color_mapping 'signature_desc|Trypanosomal VSG domain|Grey|||Name|Similar to Tb427VSG|Grey|||product|Trypanosomal VSG|Grey' \
        --x_tick_distance 500000 \
        --font_size 1 \
        --output_file ${OVERVIEW_DIR}/genome_overview.svg
    conda deactivate
}




generate_overview_pic



# STEP ?: compare the assemblies

exit # for now we skip this step

if [ ! -e ${DATA_DIR}/comparison/HGAP3_Tb427v10.sorted.fasta ];then
    fassort ${GENOME_FASTA_IN} > ${DATA_DIR}/comparison/HGAP3_Tb427v10.sorted.fasta
    fassort ${OUTPUT_ASSEMBLY} > ${DATA_DIR}/comparison/new.sorted.fasta
fi 


if [ ! -e ${MINIM_OUT_FILE}.paf ];then
    minimap2 -x asm5 ${DATA_DIR}/comparison/new.sorted.fasta ${DATA_DIR}/comparison/HGAP3_Tb427v10.sorted.fasta > ${MINIM_OUT_FILE}.paf
fi 

if [ ! -e ${DIFFERENCES_OUT_FILE}.paf ];then
    python3 ../../scripts/differences_from_paf.py ${MINIM_OUT_FILE}.paf > ${DIFFERENCES_OUT_FILE}.paf
fi 

