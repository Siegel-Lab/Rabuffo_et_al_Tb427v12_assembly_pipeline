#!/bin/bash
#  grep -E '\tfive_prime_UTR\t' | cat gff_without_five_prime_UTR.gff3 - | sort -k1,1 -k4,4n


cat <(grep '^#' ${INPUT_GFF_FILE}) <(cat <(grep -v '^#' ${INPUT_GFF_FILE}) <(grep -P '\tfive_prime_UTR\t' ${FIVE_ANNOTATION_FILE}) <(grep -P '\tthree_prime_UTR\t' ${THREE_ANNOTATION_FILE}) | sort -k1,1 -k4,4n) > ${OUTPUT_GFF_FILE}

#<(grep -E '\three_prime_UTR\t' ${THREE_ANNOTATION_FILE}) 

#cat ${INPUT_GFF_FILE} <(grep -E '\tfive_prime_UTR\t' ${FIVE_ANNOTATION_FILE} | cat <(grep -v '^#' ${INPUT_GFF_FILE}) - | sort -k1,1 -k4,4n ) > ${OUTPUT_GFF_FILE} 




#grep -E '\tfive_prime_UTR\t' ${GFF_FOLDER}/$UTRME_FOLDER/GFF/${UTRME_5_ANNOTATION_FILE%.gff}.gff | cat $GFF_FOLDER/${GFF_A_OUT} - > ${GFF_FOLDER}/${GFF_A_OUT%.gff3}_UTRann.gff3
#
#
#paste <(grep -P '\t(pseudo)?gene\t|\tmRNA\t|\tpseudogenic_transcript\t' ${GFF_FOLDER}/temp_Complete_${GENOME_11_FULLY_PHASED_FILE%.fasta}.gff3) <(grep -P '\t(pseudo)?gene\t|\tmRNA\t|\tpseudogenic_transcript\t' ${GFF_FOLDER}/temp_Complete_${GENOME_11_FULLY_PHASED_FILE%.fasta}.gff3 | cut -f 9 | sed 's/pseudogene;.*/pseudogene/' | sed 's/ID=.*Parent=/ID=/' | sed 's/;.*//' | sed 's/ID=//') | sort -t$'\t' -k10,10 -k3,3 > ${GFF_FOLDER}/${GENE_MRNA_ORDERED_FILE}
#
##paste ${GFF_FOLDER}/${GFF_IN} <(cut -f 9 ${GFF_FOLDER}/${GFF_IN} | sed 's/Name=Similar to //' | sed 's/;Note=.*//') <(paste <(cut -f 9 ${GFF_FOLDER}/${GFF_IN} | sed 's/Name=.*The alignment length is //' | sed 's/ nt, while.*//') <(cut -f 9 ${GFF_FOLDER}/${GFF_IN} | sed 's/ percent identity.*//' | sed 's/Name=.*, //') | awk '{ print $1 * $2 }') | sort -t$'\t' -k10,10 -k11r,11 > ${GFF_FOLDER}/${RANK_VSG_FILE}
