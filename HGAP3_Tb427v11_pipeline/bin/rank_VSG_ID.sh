#!/bin/bash


paste ${GFF_FOLDER}/${GFF_IN} <(cut -f 9 ${GFF_FOLDER}/${GFF_IN} | sed 's/Name=Similar to //' | sed 's/;Note=.*//') <(paste <(cut -f 9 ${GFF_FOLDER}/${GFF_IN} | sed 's/Name=.*The alignment length is //' | sed 's/ nt, while.*//') <(cut -f 9 ${GFF_FOLDER}/${GFF_IN} | sed 's/ percent identity.*//' | sed 's/Name=.*, //') | awk '{ print $1 * $2 }') | sort -t$'\t' -k10,10 -k11r,11 > ${GFF_FOLDER}/${RANK_VSG_FILE}
