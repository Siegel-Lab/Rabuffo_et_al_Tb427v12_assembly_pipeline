#!/bin/bash

awk 'BEGIN{OFS="\t"}{if ($7 == "+" && $3 == "three_prime_UTR") {$4+=1}; print $0}' ${INPUT_FOLDER}/${GFF_IN} > ${OUTPUT_FOLDER}/${GFF_OUT}

