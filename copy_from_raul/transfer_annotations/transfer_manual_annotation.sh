#!/bin/bash
#
#

module load ngs/bedtools2/2.26.0
module load ncbi-blast/2.7.1+

## Create temp folder to store temporary files

TEMP_FOLDER=${OUTPUT_FOLDER}/temp

mkdir -p ${TEMP_FOLDER}

## Get manual annotations from old genome version

grep 'Siegel_Group' $OLD_GENOME_ANNOTATION_PATH > $TEMP_FOLDER/Siegel_Group.txt

## Crop coordinates on the entries based on chromosome length
# Get chromosome lengths

python3.6 $SCRIPTS_FOLDER/seq_length.py ${OLD_GENOME_FILE_PATH} > $TEMP_FOLDER/Old_genome_length.txt

# Crop coordinates based on chromosome length

bedtools slop -i $TEMP_FOLDER/Siegel_Group.txt -g $TEMP_FOLDER/Old_genome_length.txt -b 0 > $TEMP_FOLDER/Siegel_Group_cropped.txt

## Get fasta sequences for each entry and add it to the annotation file

bedtools getfasta -tab -fi ${OLD_GENOME_FILE_PATH} -bed $TEMP_FOLDER/Siegel_Group_cropped.txt > $TEMP_FOLDER/Siegel_Group.fasta

paste $TEMP_FOLDER/Siegel_Group_cropped.txt $TEMP_FOLDER/Siegel_Group.fasta > $TEMP_FOLDER/Siegel_Group.out

## Run a modified version of Konrad's Förstner script to transfer annotation

python3.6 $SCRIPTS_FOLDER/map_annotation_via_string_match.py $GENOME_FILE_PATH $TEMP_FOLDER/Siegel_Group.out $TEMP_FOLDER $TEMP_FOLDER/Siegel_Group.gff


# Add this entries to the annotation file
cat $GFF_PATH $TEMP_FOLDER/Siegel_Group.gff > $GFF_NEW_PATH

# Remove temp folder

#rm -R ${TEMP_FOLDER}

