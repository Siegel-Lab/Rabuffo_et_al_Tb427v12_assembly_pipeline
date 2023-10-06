#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH --mem=250G
#SBATCH --time=7-00:00:00
#SBATCH --partition=slim18
#SBATCH --job-name=samba
#SBATCH -o ../data/logfiles/slurm-out-samba-%j.out

source /home/mschmidt/.miniconda3/etc/profile.d/conda.sh

conda activate ont_assembly


## STEP 1: configure


GENOME_FASTA_IN=$(realpath ../data/genome_in/HGAP3_Tb427v10/HGAP3_Tb427v10.fasta)
GFF_IN=$(realpath ../data/genome_in/HGAP3_Tb427v10/HGAP3_Tb427v10.gff)
ONT_READS_IN=$(realpath ../data/ont_reads_in/duplex_reads.fastq.gz)
SAMBA_OUT=$(realpath ../data/samba_out)
MASUCRA_BIN=$(realpath ../bin/MaSuRCA-4.1.0/bin)
BLAST_OUT_FILE=$(realpath ../data/blast_out/blast_out.crunch)

DATA_DIR=$(realpath ../data)

## STEP 2: run samba

cd ${SAMBA_OUT}

${MASUCRA_BIN}/close_scaffold_gaps.sh -r ${GENOME_FASTA_IN} -q ${ONT_READS_IN} -t 18 -m 2000 -d ont

# .valid_join_pairs.txt contains how scaffolds have been split into contigs
# .fasta.split.joined.fa is the main outputfile

# .patches.coords seems interesting

OUTPUT_ASSEMBLY=HGAP3_Tb427v10.fasta.split.joined.fa

# count the number of Ns in the input genome
echo "N's in input genome (in this assembly a gap consists of 1,000 Ns):"
grep -o "N" ${GENOME_FASTA_IN} | wc -l

# count the number of Ns in the output genome
echo "N's in output genome (here, a gap consists of 100 Ns):"
grep -o "N" ${OUTPUT_ASSEMBLY} | wc -l

############################################
## We have 49 gaps in the original assembly
## And 17 in the fixed one
############################################


# count the number of contigs in the original assembly
echo "Number of contigs in original assembly:"
grep -c ">" ${GENOME_FASTA_IN}

# count the number of contigs in the output assembly
echo "Number of contigs in output assembly:"
grep -c ">" ${OUTPUT_ASSEMBLY}

############################################
## We have 317 contigs in the original assembly
## And 308 in the fixed one
## 
## In total 32 gaps have been closed?, 9 using a unitig present in the original assembly
##
############################################

# STEP 3: compare the assemblies
#pyfastx extract ${GENOME_FASTA_IN} Chr6_3A_Tb427v10 > ${DATA_DIR}/Chr6_3A_Tb427v10.input.fasta
#pyfastx extract ${OUTPUT_ASSEMBLY} Chr6_3A_Tb427v10 > ${DATA_DIR}/Chr6_3A_Tb427v10.output.fasta

cd ..


#blastn -outfmt 6 -num_threads=18 -query=Chr6_3A_Tb427v10.output.fasta -db=Chr6_3A_Tb427v10.input.fasta -evalue 1 -task megablast -out=${BLAST_OUT_FILE}

#liftoff -g ${GFF_IN} 