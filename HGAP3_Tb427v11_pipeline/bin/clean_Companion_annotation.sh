#!/bin/bash
#
#

module load ngs/bedtools2/2.26.0

## Create temp folder to store temporary files

TEMP_FOLDER=${OUTPUT_FOLDER}/temp

mkdir -p ${TEMP_FOLDER}

## Generate list of rRNA entries

grep -P '\trRNA\t' ${INPUT} > ${TEMP_FOLDER}/rRNA_gene_list.txt

## Generate list of pseudogene (protein-coding genes that are pseudogenes)

grep -P '\tpseudogene\t' ${INPUT} > ${TEMP_FOLDER}/pseudogene_list.txt

## Intersect the lists and get a list with only the pseudogenes overlapping with rRNA genes

bedtools intersect -wa -nonamecheck -a ${TEMP_FOLDER}/pseudogene_list.txt -b ${TEMP_FOLDER}/rRNA_gene_list.txt > ${TEMP_FOLDER}/overlapping.txt

## Filter the entries to get the list of IDs to remove

grep -oP "ID=.*?:" ${TEMP_FOLDER}/overlapping.txt | sed s/:$// > ${TEMP_FOLDER}/GeneIDs_to_remove.txt

## Remove entries from the Companion annotation file based on the list of IDs generated

grep -vf ${TEMP_FOLDER}/GeneIDs_to_remove.txt ${INPUT} > ${TEMP_FOLDER}/temp.gff3


# Generate list of entries from small protein coding genes

awk -v MIN_LENGTH="${MIN_LENGTH}" '{if (($5 -$4)< MIN_LENGTH && ($3 == "polypeptide")) print $0}' ${TEMP_FOLDER}/temp.gff3 > $TEMP_FOLDER/Small_protein_entries.txt

# Get the gene IDs

grep -oP "ID=.*?:" ${TEMP_FOLDER}/Small_protein_entries.txt | sed s/\.1:$// > ${TEMP_FOLDER}/Small_proteinIDs_to_remove.txt

# Remove the entries with that gene ID

grep -vf ${TEMP_FOLDER}/Small_proteinIDs_to_remove.txt ${TEMP_FOLDER}/temp.gff3 > ${OUTPUT}


## Remove temp folder

#rm -R ${TEMP_FOLDER}

