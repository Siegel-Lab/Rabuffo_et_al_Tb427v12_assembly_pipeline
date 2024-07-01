#!/bin/bash
#SBATCH --partition=slim16 		# (queue, see sinfo) 
#SBATCH --ntasks=16			# number of threads
#SBATCH --mem 128G 				#memory pool for all cores 
#SBATCH -o slurm.%u.%j.out 		#STDOUT 
#SBATCH -e slurm.%u.%j.err 		#STDERR

module load ngs/bwa/0.7.16 
module load ngs/samtools/1.8

# Check the file ending to see which decompression tool to use (currently either bunzip2 or zcat)

if [ ${INPUT_FILE_READ_1: -4} == ".bz2" ] # checks the four last characters of the variables to much with ".bz2"
then
	decompresser="bunzip2 -c"
else 
	decompresser="zcat"
fi

# Map reads to the genome index
bwa mem -t 16 \
	${GENOME_INDEX} \
	<(${decompresser} ${INPUT_FILE_READ_1}) \
	<(${decompresser} ${INPUT_FILE_READ_2}) \
	| samtools view -bh -F 3332 > ${OUTPUT_FILE}
wait

# sort and index the bam files
SORTED_BAM_FILE=$(echo ${OUTPUT_FILE} | sed "s/.bam$/.sorted.bam/")
samtools sort \
	-@ 16 \
	-o $SORTED_BAM_FILE ${OUTPUT_FILE} && \
samtools index $SORTED_BAM_FILE 
wait

# Remove duplicates with picard MarkDuplicates

OUTPUT_FILE_REMOVE_DUP_METRICS=${SORTED_BAM_FILE%sorted.bam}rmdup.metrics.txt
OUTPUT_FILE_REMOVE_DUP_METRICS=${OUTPUT_FILE_REMOVE_DUP_METRICS##*/}
java -jar $PICARD MarkDuplicates \
        I=$SORTED_BAM_FILE \
        O=${SORTED_BAM_FILE%sorted.bam}rmdup.sorted.bam \
        M=$REMOVE_DUP_METRICS_FOLDER/${OUTPUT_FILE_REMOVE_DUP_METRICS} \
        REMOVE_DUPLICATES=true \
        ASSUME_SORT_ORDER=coordinate \
        QUIET=true \
        VERBOSITY=WARNING && \
samtools index ${SORTED_BAM_FILE%sorted.bam}rmdup.sorted.bam
wait

# clean up
rm ${OUTPUT_FILE}

