#!/bin/bash
#SBATCH --partition=slim16              # (queue, see sinfo) 
#SBATCH --ntasks=16                     # number of threads
#SBATCH --mem 128G                      #memory pool for all cores 
#SBATCH -o slurm.%u.%j.out              #STDOUT 
#SBATCH -e slurm.%u.%j.err              #STDERR

#
# This module aims to generate an updated version of the allele-specific genome assembly of T. brucei Lister 427 to upload to EupathDB. The version will be 'Tb427v11'
#
#to run ./run.sh
#
# conda env create -f 2022-10-11-Generate-updated-haplotype-Lister427-genomes-for-EupathDB_env.yaml

# conda env update --file 2022-10-11-Generate-updated-haplotype-Lister427-genomes-for-EupathDB_env.yaml --prune

source /work/project/ladsie_003/miniconda3/etc/profile.d/conda.sh
conda activate 2022-10-11-Generate-updated-haplotype-Lister427-genomes-for-EupathDB_env # to use fastaq, blast, samtools, biopython, bioperl, seqtk

#module load R/3.5.0

#module load Anaconda3/2019.03
#source activate 2022-10-11-Generate-updated-haplotype-Lister427-genomes-for-EupathDB_env # to use fastaq, blast, samtools, biopython, bioperl

main(){
	set_variables
	create_folders
	link_input_genome_files

#	generate_samtools_ix
#	delete_wrong_region_from_BES2
#	adapt_BES2_gff
#	correct_intron_annotations
#	generate_blast_db
#	blast_VSGs_to_genome
#	blast2gff
#	clean_Companion_gff
#	merge_VSG_Companion_gff

#	link_reads
#	edit_configuration_file

#	conda deactivate
#	conda activate /home/rcosentino/.conda/envs/UTRme
#	run_UTRme
#	conda deactivate
#	conda activate R_env
#	add_UTR_annotation

#	conda deactivate
#	conda activate 2022-10-11-Generate-updated-haplotype-Lister427-genomes-for-EupathDB_env

#	scaffold_fasta_hap_chroms
#	transfer_manual_annotation
#	scaffold_gffs_hap_chroms

#	bwa_index_genome

#	map_RNAseq_to_genome



}

###############################################################################
# Setting of global variables
###############################################################################

set_variables(){
	INPUT_FOLDER=input
	OUTPUT_FOLDER=output
	SCRIPTS_FOLDER=bin

	#input files:
	GENOME_10_FULLY_PHASED_SOURCE=/ladsie/project/ladsie_003/TbAllelePhasingPaperPipelines/06_Fully_phased_scaffold/input
	GENOME_10_FULLY_PHASED_FILE=HGAP3_Tb427v10_phased_diploid_core.fasta
	ANNOTATION_10_FULLY_PHASED_FILE=Companion_927v41_HGAP3_Tb427v10_phased_diploid_core.gff3

	#to get coordinates for the transfering of manual annotations:
	GENOME_10_FILE=HGAP3_Tb427v10.fasta
	GENOME_10_SOURCE=/ladsie/project/ladsie_003/TbAllelePhasingPaperPipelines/04_Genome_annotation_conversion_v9.9_to_v10/output
	ANNOTATION_10_MANUAL_FILE=HGAP3_Tb427v10_manual.gff3
	ANNOTATION_10_MANUAL_SOURCE=/ladsie/project/ladsie_003/TbAllelePhasingPaperPipelines/04_Genome_annotation_conversion_v9.9_to_v10/output
	ANNOTATION_10_PTU_FILE=PTUs_HGAP3_Tb427v10_sorted_with_IDs.gff3
	ANNOTATION_10_PTU_SOURCE=/ladsie/project/ladsie_003/2019-projects/2019-03-11-PTU-annotation/output
	#output files:
	GENOME_11_FULLY_PHASED_FILE=Tb427v11_phased_diploid_core.fasta
	GENOME_11_PHASE_A_FILE=Tb427v11_hapA_scaffolded.fasta # file containing the scaffolded hapA of the eleven chromosomes
	GENOME_11_PHASE_B_FILE=Tb427v11_hapB_scaffolded.fasta # file containing the scaffolded hapB of the eleven chromosomes
	GENOME_11_PHASE_A_FULL_FILE=Tb427v11_hapA_full_not_scaffolded.fasta # file containing the hapA cores, one copy of each subtelomeric contig (including those present in both alleles), all the BESs, the minichromosomes, the unassigned contigs and the maxicircle.

	ANNOTATION_11_FULLY_PHASED_FILE=Tb427v11_phased_diploid_core.gff3
	ANNOTATION_11_PHASE_A_FILE=Tb427v11_hapA_scaffolded.gff3
	ANNOTATION_11_PHASE_B_FILE=Tb427v11_hapB_scaffolded.gff3
	ANNOTATION_11_PHASE_FULL_FILE=Tb427v11_hapA_full_not_scaffolded.gff3

	GENOME_FOLDER=${OUTPUT_FOLDER}/genomes
	BLAST_FOLDER=${OUTPUT_FOLDER}/blast
	GFF_FOLDER=${OUTPUT_FOLDER}/annotation

	VSG_FASTA_FILE=VSGs_CDS_Cross_list.fasta # it is also located in $GENOME_10_FULLY_PHASED_SOURCE
	GUIDE_SCAFFOLDING_FILE=guide_for_scaffolding_fully_phased_genome_v11.txt	
	#create read folder
	READS_SOURCE=/ladsie/project/ladsie_003/data/Kraus_et_al_RNAseq_data
	JEACOCK_READ_SOURCE=/ladsie/project/ladsie_003/data/Jeacock_et_al_2018_seq_data
	READS_FOLDER=${INPUT_FOLDER}/reads
	READ_ONE_FOLDER=${READS_FOLDER}/read_one
	READ_TWO_FOLDER=${READS_FOLDER}/read_two


}

###############################################################################
# Prelude
###############################################################################

create_folders(){
	mkdir -p \
	$SCRIPTS_FOLDER \
	$INPUT_FOLDER \
	$OUTPUT_FOLDER \
	$GENOME_FOLDER \
	$BLAST_FOLDER \
	$GFF_FOLDER \
	$READS_FOLDER \
	$READ_ONE_FOLDER \
	$READ_TWO_FOLDER 

}

link_input_genome_files(){

	ln -sf $GENOME_10_FULLY_PHASED_SOURCE/$GENOME_10_FULLY_PHASED_FILE $INPUT_FOLDER
	ln -sf $GENOME_10_FULLY_PHASED_SOURCE/$ANNOTATION_10_FULLY_PHASED_FILE $INPUT_FOLDER

	ln -sf $GENOME_10_SOURCE/$GENOME_10_FILE $INPUT_FOLDER
	ln -sf $ANNOTATION_10_MANUAL_SOURCE/$ANNOTATION_10_MANUAL_FILE $INPUT_FOLDER
	ln -sf $ANNOTATION_10_PTU_SOURCE/$ANNOTATION_10_PTU_FILE $INPUT_FOLDER
}


generate_samtools_ix(){

jid2=$(sbatch --job-name=samix --wrap="samtools faidx $INPUT_FOLDER/$GENOME_10_FULLY_PHASED_FILE")
jid2=$(echo $jid2 | sed "s/Submitted batch job //")


}

delete_wrong_region_from_BES2(){

#BES2
	#remove 50094-63987
#generate fasta files from those regions using samtools


jid3=$(sbatch --job-name=BES2_1st --wrap="samtools faidx $INPUT_FOLDER/$GENOME_10_FULLY_PHASED_FILE BES2_Tb427v10:1-50093 > ${GENOME_FOLDER}/temp_1st_part_BES2_v10.fasta")
jid3=$(echo $jid3 | sed "s/Submitted batch job //")

jid4=$(sbatch --job-name=BES2_2nd --wrap="samtools faidx $INPUT_FOLDER/$GENOME_10_FULLY_PHASED_FILE BES2_Tb427v10:63988- > ${GENOME_FOLDER}/temp_2nd_part_BES2_v10.fasta")
jid4=$(echo $jid4 | sed "s/Submitted batch job //")

# concat pieces into a file (adding a fasta sequence with 1000 Ns in between)

jid5=$(sbatch --job-name=cat_BES2 --dependency=afterok:${jid3}:${jid4} --wrap="cat ${GENOME_FOLDER}/temp_1st_part_BES2_v10.fasta $INPUT_FOLDER/1000_Ns.fasta ${GENOME_FOLDER}/temp_2nd_part_BES2_v10.fasta > ${GENOME_FOLDER}/temp_parts_BES2_v10.fasta")
jid5=$(echo $jid5 | sed "s/Submitted batch job //")

#merge into one sequence

jid6=$(sbatch --job-name=join-parts --dependency=afterok:${jid5} --wrap="fastaq merge --name BES2_Tb427v11 ${GENOME_FOLDER}/temp_parts_BES2_v10.fasta ${GENOME_FOLDER}/temp_fixed_BES2_v11.fasta")
jid6=$(echo $jid6 | sed "s/Submitted batch job //")

#replace old BES sequence with the new one and rename each chromosome to v11 and also the full genome

jid7=$(sbatch --job-name=rm-old-chr --wrap="fastaq filter --regex BES2_Tb427v10 -v ${INPUT_FOLDER}/${GENOME_10_FULLY_PHASED_FILE} ${GENOME_FOLDER}/temp_${GENOME_10_FULLY_PHASED_FILE%.fasta}_minus_BES2.fasta")
jid7=$(echo $jid7 | sed "s/Submitted batch job //")


jid8=$(sbatch --job-name=join-new --dependency=afterok:${jid6}:${jid7} --wrap="cat ${GENOME_FOLDER}/temp_fixed_BES2_v11.fasta ${GENOME_FOLDER}/temp_${GENOME_10_FULLY_PHASED_FILE%.fasta}_minus_BES2.fasta > ${GENOME_FOLDER}/temp_${GENOME_10_FULLY_PHASED_FILE%.fasta}_fixed_BES2.fasta")
jid8=$(echo $jid8 | sed "s/Submitted batch job //")

#shorten the contigs name 's/containing//g' and 's/restricted//g'
#replace chrom name : "Tb427VSG19_containing_unitig_Tb427v11" for "Tb427VSG-19_unitig_Tb427v11"
jid9=$(sbatch --job-name=version-update --dependency=afterok:${jid8} --wrap="sed 's/containing_//g' ${GENOME_FOLDER}/temp_${GENOME_10_FULLY_PHASED_FILE%.fasta}_fixed_BES2.fasta | sed 's/restricted_//g' | sed 's/Tb427VSG19/Tb427VSG-19/g' | sed 's/Tb427v10/Tb427v11/g' > ${GENOME_FOLDER}/${GENOME_11_FULLY_PHASED_FILE}")
jid9=$(echo $jid9 | sed "s/Submitted batch job //")


}

# Modify BES2 gff based on the sequence remove --> shift coordinates of genes downstream of the cut region, and delete entries within the cut region

adapt_BES2_gff(){
# Maybe i could do it with a perl script..
# go through the annotation file
# print every line with another chromosome
# print every line with of BES2 with a start position before the cut
# don't print lines with a BES2 start bigger than X and smaller than Y (that would be enough because no entries overlap the start and end position)
# for lines in BES2 bigger than Y, shift start and end coordinates by (Y - X) - 1000 --> double check if at the end there is no 1nt shift between annotation and sequence

# usage = perl adapt_gff.pl input.gff3 chrom_name cut_start cut_end new_gap_length > adapted.gff3
jid7=$(sbatch --job-name=adapt-gff  --wrap="perl ${SCRIPTS_FOLDER}/adapt_gff_after_deleting_chrom_region.pl ${INPUT_FOLDER}/${ANNOTATION_10_FULLY_PHASED_FILE} BES2_Tb427v10 50094 63987 1000 > ${GFF_FOLDER}/temp_${ANNOTATION_10_FULLY_PHASED_FILE}_fixed_BES2.gff3")
jid7=$(echo $jid7 | sed "s/Submitted batch job //")

NEW_TEMP_GFF_FILE=$(echo temp_${ANNOTATION_10_FULLY_PHASED_FILE}_fixed_BES2.gff3 | sed 's/Tb427v10/Tb427v11/g')


# add some sed to modify the naming of the chromosomes matching the fasta file: sed 's/containing_//g' ${GFF_FOLDER}/temp_${ANNOTATION_10_FULLY_PHASED_FILE}_fixed_BES2.gff3 | sed 's/restricted_//g' | sed 's/Tb427VSG19/Tb427VSG-19/g' | sed 's/Tb427v10/Tb427v11/g' > ${GFF_FOLDER}/${NEW_TEMP_GFF_FILE}
jid8=$(sbatch --job-name=adapt-gff --dependency=afterok:${jid7} --wrap="sed 's/containing_//g' ${GFF_FOLDER}/temp_${ANNOTATION_10_FULLY_PHASED_FILE}_fixed_BES2.gff3 | sed 's/restricted_//g' | sed 's/Tb427VSG19/Tb427VSG-19/g' | sed 's/Tb427v10/Tb427v11/g' > ${GFF_FOLDER}/${NEW_TEMP_GFF_FILE}")
jid8=$(echo $jid8 | sed "s/Submitted batch job //")




}

#       (gff file): Fix intron problems based on previous script

correct_intron_annotations(){

# Generate a list of the genes that are in this situation using the presence of "CDS:2" on the last column of the gff as evidence of gene that has an intron annotated (more than one CDS)

ID_LIST=ID_list_temp.txt

GFF_IN_FILE=$(echo temp_${ANNOTATION_10_FULLY_PHASED_FILE}_fixed_BES2.gff3 | sed 's/Tb427v10/Tb427v11/g')

jid1=$(sbatch --job-name=ID_list --wrap="grep 'CDS:2' $GFF_FOLDER/$GFF_IN_FILE | cut -f 9 | cut -f 1 -d '.' > $GFF_FOLDER/${ID_LIST}")
jid1=$(echo $jid1 | sed "s/Submitted batch job //")

GFF_OUT_FILE=$(echo $GFF_IN_FILE | sed 's/\.gff3$/\.1\.gff3/')

jid2=$(sbatch --job-name=correct_introns --dependency=afterok:${jid1} --wrap="perl $SCRIPTS_FOLDER/correct_Companion_intron_annotation_v2.pl $GFF_FOLDER/$GFF_IN_FILE $GFF_FOLDER/${ID_LIST} > $GFF_FOLDER/$GFF_OUT_FILE")
jid2=$(echo $jid2 | sed "s/Submitted batch job //")

#clean-up
jid3=$(sbatch --job-name=correct_introns --dependency=singleton --wrap="rm $GFF_FOLDER/ID_list_temp.txt")
jid3=$(echo $jid3 | sed "s/Submitted batch job //")


}



generate_blast_db(){
jid1=$(sbatch --job-name=blastdb --wrap="makeblastdb -in ${GENOME_FOLDER}/${GENOME_11_FULLY_PHASED_FILE} -dbtype nucl -parse_seqids")
jid1=$(echo $jid1 | sed "s/Submitted batch job //")
}

blast_VSGs_to_genome(){

OUTPUT_FILE="${VSG_FASTA_FILE%.fasta}_to_${GENOME_11_FULLY_PHASED_FILE%.fasta}.out"

jid2=$(sbatch --job-name=runBlast --wrap="blastn -db $GENOME_FOLDER/$GENOME_11_FULLY_PHASED_FILE -query $INPUT_FOLDER/$VSG_FASTA_FILE -outfmt 7 -out $BLAST_FOLDER/${OUTPUT_FILE}")
jid2=$(echo $jid2 | sed "s/Submitted batch job //")

}


blast2gff(){

BLAST_OUTPUT_FILE=$BLAST_FOLDER/"${VSG_FASTA_FILE%.fasta}_to_${GENOME_11_FULLY_PHASED_FILE%.fasta}.out"

jid3=$(sbatch --job-name=blast2gff --wrap="perl $SCRIPTS_FOLDER/Blast_output_to_gff3_v2.pl $BLAST_OUTPUT_FILE $INPUT_FOLDER/$VSG_FASTA_FILE > $GFF_FOLDER/VSGs_in_${GENOME_11_FULLY_PHASED_FILE%.fasta}.gff3")
jid3=$(echo $jid3 | sed "s/Submitted batch job //")

jid4=$(sbatch --job-name=sort_gff --dependency=afterok:${jid3} --wrap="sort -k1,1 -k4,4n -o $GFF_FOLDER/VSGs_in_${GENOME_11_FULLY_PHASED_FILE%.fasta}.gff3 $GFF_FOLDER/VSGs_in_${GENOME_11_FULLY_PHASED_FILE%.fasta}.gff3")
jid4=$(echo $jid4 | sed "s/Submitted batch job //")

}


clean_Companion_gff(){

#Remove pseudogene entries overlapping rRNA entries and small protein-coding gene entries
GFF_IN_FILE=${GFF_FOLDER}/$(echo temp_${ANNOTATION_10_FULLY_PHASED_FILE}_fixed_BES2.gff3 | sed 's/Tb427v10/Tb427v11/g' | sed 's/\.gff3$/\.1\.gff3/')
GFF_OUT_FILE=$(echo $GFF_IN_FILE | sed 's/\.1\.gff3$/\.2\.gff3/')
MIN_LENGTH=60 # Nucleotides


jid5=$(sbatch --job-name=rm_pseudo_rRNA_overlapping --export=INPUT=${GFF_IN_FILE},MIN_LENGTH=${MIN_LENGTH},OUTPUT_FOLDER=${GFF_FOLDER},OUTPUT=${GFF_OUT_FILE} ${SCRIPTS_FOLDER}/clean_Companion_annotation.sh)
jid5=$(echo $jid5 | sed "s/Submitted batch job //")

}


merge_VSG_Companion_gff(){

VSG_GFF=VSGs_in_${GENOME_11_FULLY_PHASED_FILE%.fasta}.gff3
COMPANION_GFF=$(echo temp_${ANNOTATION_10_FULLY_PHASED_FILE}_fixed_BES2.gff3 | sed 's/Tb427v10/Tb427v11/g' | sed 's/\.gff3$/\.2\.gff3/')
OUT_GFF=$(echo temp_${ANNOTATION_10_FULLY_PHASED_FILE}_fixed_BES2.gff3 | sed 's/Tb427v10/Tb427v11/g' | sed 's/\.gff3$/\.3\.gff3/')

#merge VSG-overlapping entries in the BLAST result

jid6=$(sbatch --job-name=join_overlapping --wrap="perl $SCRIPTS_FOLDER/join_overlapping_entries_VSGs_gff3_v3.pl $GFF_FOLDER/${VSG_GFF} > ${GFF_FOLDER}/${VSG_GFF%.gff3}_overlapping_joined.gff3")
jid6=$(echo $jid6 | sed "s/Submitted batch job //")

# modify the VSG_ID in case of repetition, so in case of multiple give Tb427VSG-XX to the one with highest similarity, and then Tb427VSG-XX.1 ; Tb427VSG-XX.2 to the rest..

#generate a table with a column with the VSG_ID and the Aln_length * perc_similarity , which will be used to rank the genes with the same VSG_ID
RANK_VSG_FILE=temp_VSG_rank_name.txt
GFF_IN=${VSG_GFF%.gff3}_overlapping_joined.gff3
jid7=$(sbatch --job-name=VSG_ID --dependency=afterok:${jid6} --export=GFF_FOLDER=${GFF_FOLDER},GFF_IN=${GFF_IN},RANK_VSG_FILE=${RANK_VSG_FILE} ${SCRIPTS_FOLDER}/rank_VSG_ID.sh)
jid7=$(echo $jid7 | sed "s/Submitted batch job //")

#modify the VSG_IDs in the VSG rank file adding .1 .2 to the successive genes with the same name in the list, add "description=" in front and use that as the attributes column, and remove the previous attribute column and the score column, converting it to a gff file again

RANK_SUFFIX_VSG_GFF_FILE=temp_VSG_rank_name_suffix_added.gff3
jid8=$(sbatch --job-name=join_overlapping --dependency=afterok:${jid7} --wrap="perl $SCRIPTS_FOLDER/add_suffix_to_ranked_VSG_IDs.pl $GFF_FOLDER/${RANK_VSG_FILE} > ${GFF_FOLDER}/${RANK_SUFFIX_VSG_GFF_FILE}")
jid8=$(echo $jid8 | sed "s/Submitted batch job //")


#join VSG_gff with Companion_gff, remove comment lines and sort them

jid9=$(sbatch --job-name=join_VSG_and_Companion_gffs --dependency=afterok:${jid8} --wrap="cat ${GFF_FOLDER}/${RANK_SUFFIX_VSG_GFF_FILE} ${GFF_FOLDER}/${COMPANION_GFF} | grep -v '^#' | sort -k1,1 -k4,4n > ${GFF_FOLDER}/temp_joined_VSG_Companion_in_${GENOME_11_FULLY_PHASED_FILE%.fasta}.gff3")
jid9=$(echo $jid9 | sed "s/Submitted batch job //")


#merge VSG entries with Companion_gff entries

jid10=$(sbatch --job-name=merge_VSG_and_gene_Companion_entries --dependency=afterok:${jid9} --wrap="perl $SCRIPTS_FOLDER/merge_VSG_Companion_gff3.pl ${GFF_FOLDER}/temp_joined_VSG_Companion_in_${GENOME_11_FULLY_PHASED_FILE%.fasta}.gff3 > ${GFF_FOLDER}/temp_merged_VSG_Companion_in_${GENOME_11_FULLY_PHASED_FILE%.fasta}.gff3")
jid10=$(echo $jid10 | sed "s/Submitted batch job //")


#add geneID to VSG entries from the BLAST result that were not merged with Companion VSG results. should i add transcript entries? pseudogenic_transcript? --> Yes
#should i add a ":pseudogene" at the end? --> Yes
GENE_ID_PREFIX=Tb427_000
FIRST_NEW_GENE_ID=845300

jid11=$(sbatch --job-name=Add_GeneIDs_orphan_VSG_entries --dependency=afterok:${jid10} --wrap="perl $SCRIPTS_FOLDER/Add_IDs_to_VSGs_v2.pl ${GFF_FOLDER}/temp_merged_VSG_Companion_in_${GENOME_11_FULLY_PHASED_FILE%.fasta}.gff3 $GENE_ID_PREFIX $FIRST_NEW_GENE_ID > ${GFF_FOLDER}/temp_Complete_${GENOME_11_FULLY_PHASED_FILE%.fasta}.gff3 ")
jid11=$(echo $jid11 | sed "s/Submitted batch job //")


#change pseudogenic_transcript mRNA coordinates to match the ones of the gene
	# generate a file containing the coordinates and IDs of all genes and pseudogenes
GENE_MRNA_ORDERED_FILE=temp_gene_mRNA_ordered.txt
jid12=$(sbatch --job-name=gene_mRNA --dependency=afterok:${jid11} --export=GFF_FOLDER=${GFF_FOLDER},GENOME_11_FULLY_PHASED_FILE=${GENOME_11_FULLY_PHASED_FILE},GENE_MRNA_ORDERED_FILE=${GENE_MRNA_ORDERED_FILE} ${SCRIPTS_FOLDER}/gene_mRNA_ordered.sh)
jid12=$(echo $jid12 | sed "s/Submitted batch job //")


	# modify the pseudogenic_transcript or mRNA coordinates of those that are not the same than their corresponding genes/pseudogenes
jid13=$(sbatch --job-name=adjust_transc_coord --dependency=afterok:${jid12} --wrap="perl $SCRIPTS_FOLDER/adjust_transc_coord.pl ${GFF_FOLDER}/${GENE_MRNA_ORDERED_FILE} > ${GFF_FOLDER}/temp_Complete_trans_coord_adjusted_${GENOME_11_FULLY_PHASED_FILE%.fasta}_only_genes_and_transcript.gff3")
jid13=$(echo $jid13 | sed "s/Submitted batch job //")

	# re-join the gene and mRNA entries (with coordinates corrected) with the rest of the gff entries and re-add the header
GFF_IN_ONE=temp_Complete_${GENOME_11_FULLY_PHASED_FILE%.fasta}.gff3
GFF_IN_TWO=temp_Complete_trans_coord_adjusted_${GENOME_11_FULLY_PHASED_FILE%.fasta}_only_genes_and_transcript.gff3
jid14=$(sbatch --job-name=re_join_adjusted --dependency=afterok:${jid13} --export=GFF_FOLDER=${GFF_FOLDER},COMPANION_GFF=${COMPANION_GFF},GFF_IN_ONE=${GFF_IN_ONE},GFF_IN_TWO=${GFF_IN_TWO},OUT_GFF=${OUT_GFF} ${SCRIPTS_FOLDER}/re_join_corrected_transcripts.sh)
jid14=$(echo $jid14 | sed "s/Submitted batch job //")


}

### ADD UTR ANNOTATION!! 
# link RNA-seq data to do the UTR annotation
link_reads(){

# link reads from Kraus_et_al_2018
for READ_ONE in $(ls $READS_SOURCE | grep '3wash\|ssRNAd_WT' | grep '_R1.fq.gz$')
do

	ln -sf $READS_SOURCE/$READ_ONE $READ_ONE_FOLDER
	
done

for READ_TWO in $(ls $READS_SOURCE | grep '3wash\|ssRNAd_WT' | grep '_R2.fq.gz$')
do

	ln -sf $READS_SOURCE/$READ_TWO $READ_TWO_FOLDER
	
done
# link reads from Jeacock_et_al_2018
for READ_ONE in $(ls $JEACOCK_READ_SOURCE | grep '^ERR' | grep '_1.fastq.gz$')
do
        ln -sf $JEACOCK_READ_SOURCE/$READ_ONE $READ_ONE_FOLDER/${READ_ONE%_1.fastq.gz}_R1.fq.gz
done

for READ_TWO in $(ls $JEACOCK_READ_SOURCE | grep '^ERR' | grep '_2.fastq.gz$')
do
        ln -sf $JEACOCK_READ_SOURCE/$READ_TWO $READ_TWO_FOLDER/${READ_TWO%_2.fastq.gz}_R2.fq.gz
done


}

edit_configuration_file(){
GFF_IN=$(echo temp_${ANNOTATION_10_FULLY_PHASED_FILE}_fixed_BES2.gff3 | sed 's/Tb427v10/Tb427v11/g' | sed 's/\.gff3$/\.3\.gff3/')

# Generate conf files removing cores of one haplotype (improves the annotation results)
GFF_IN_NO_CORE_B=${GFF_IN%.gff3}_no_core_B.gff3
GFF_IN_NO_CORE_A=${GFF_IN%.gff3}_no_core_A.gff3
grep -v '_core_Tb427v11_B' ${GFF_FOLDER}/${GFF_IN} > ${GFF_FOLDER}/${GFF_IN_NO_CORE_B}
grep -v '_core_Tb427v11_A' ${GFF_FOLDER}/${GFF_IN} > ${GFF_FOLDER}/${GFF_IN_NO_CORE_A}

jid1=$(sbatch --job-name= --wrap="fastaq filter --regex _core_Tb427v11_B -v ${GENOME_FOLDER}/$GENOME_11_FULLY_PHASED_FILE ${GENOME_FOLDER}/${GENOME_11_FULLY_PHASED_FILE%.fasta}_no_core_B.fasta")
jid1=$(echo $jid1 | sed "s/Submitted batch job //")

jid2=$(sbatch --job-name= --wrap="fastaq filter --regex _core_Tb427v11_A -v ${GENOME_FOLDER}/$GENOME_11_FULLY_PHASED_FILE ${GENOME_FOLDER}/${GENOME_11_FULLY_PHASED_FILE%.fasta}_no_core_A.fasta")
jid2=$(echo $jid2 | sed "s/Submitted batch job //")

# conf file removing core_B
cp /home/rcosentino/UTRme/Configuration_Files/Example_configuration_file.txt $INPUT_FOLDER/Configuration_no_core_B_file.txt

sed -i 's/T. cruzi/T. brucei/g' $INPUT_FOLDER/Configuration_no_core_B_file.txt
sed -i "s#/home/Maria/Reference/annotation.gff#${PWD}/${GFF_FOLDER}/${GFF_IN_NO_CORE_B}#g" $INPUT_FOLDER/Configuration_no_core_B_file.txt
sed -i "s#/home/Maria/Reference/genome.fasta#${PWD}/${GENOME_FOLDER}/${GENOME_11_FULLY_PHASED_FILE%.fasta}_no_core_B.fasta#g" $INPUT_FOLDER/Configuration_no_core_B_file.txt
sed -i "s#adapter_value:AGATCGGAAGAGC#adapter_value:NO#g" $INPUT_FOLDER/Configuration_no_core_B_file.txt
sed -i "s#feature_value:CDS#feature_value:mRNA pseudogenic_transcript#" $INPUT_FOLDER/Configuration_no_core_B_file.txt
sed -i "s#/home/Maria/Experiment/FASTQ2/#${PWD}/${READ_TWO_FOLDER}/#g" $INPUT_FOLDER/Configuration_no_core_B_file.txt
sed -i "s#/home/Maria/Experiment/FASTQ1/#${PWD}/${READ_ONE_FOLDER}/#g" $INPUT_FOLDER/Configuration_no_core_B_file.txt
echo "" >> $INPUT_FOLDER/Configuration_no_core_B_file.txt
echo "output_directory:${PWD}/${GFF_FOLDER}/UTRme_no_core_B_results" >> $INPUT_FOLDER/Configuration_no_core_B_file.txt

# conf file removing core_A

cp /home/rcosentino/UTRme/Configuration_Files/Example_configuration_file.txt $INPUT_FOLDER/Configuration_no_core_A_file.txt

sed -i 's/T. cruzi/T. brucei/g' $INPUT_FOLDER/Configuration_no_core_A_file.txt
sed -i "s#/home/Maria/Reference/annotation.gff#${PWD}/${GFF_FOLDER}/${GFF_IN_NO_CORE_A}#g" $INPUT_FOLDER/Configuration_no_core_A_file.txt
sed -i "s#/home/Maria/Reference/genome.fasta#${PWD}/${GENOME_FOLDER}/${GENOME_11_FULLY_PHASED_FILE%.fasta}_no_core_A.fasta#g" $INPUT_FOLDER/Configuration_no_core_A_file.txt
sed -i "s#adapter_value:AGATCGGAAGAGC#adapter_value:NO#g" $INPUT_FOLDER/Configuration_no_core_A_file.txt
sed -i "s#feature_value:CDS#feature_value:mRNA pseudogenic_transcript#" $INPUT_FOLDER/Configuration_no_core_A_file.txt
sed -i "s#/home/Maria/Experiment/FASTQ2/#${PWD}/${READ_TWO_FOLDER}/#g" $INPUT_FOLDER/Configuration_no_core_A_file.txt
sed -i "s#/home/Maria/Experiment/FASTQ1/#${PWD}/${READ_ONE_FOLDER}/#g" $INPUT_FOLDER/Configuration_no_core_A_file.txt
echo "" >> $INPUT_FOLDER/Configuration_no_core_A_file.txt
echo "output_directory:${PWD}/${GFF_FOLDER}/UTRme_no_core_A_results" >> $INPUT_FOLDER/Configuration_no_core_A_file.txt

}

run_UTRme(){

	jid7=$(sbatch --job-name=clean --wrap="rm -rf ${PWD}/${GFF_FOLDER}/UTRme_no_core_B_results")
	jid7=$(echo $jid7 | sed "s/Submitted batch job //")


	jid8=$(sbatch --job-name=UTRme --mem=60G --dependency=afterok:${jid7} --wrap="python ${HOME}/UTRme/utrme.py $INPUT_FOLDER/Configuration_no_core_B_file.txt")
	jid8=$(echo $jid8 | sed "s/Submitted batch job //")

	jid9=$(sbatch --job-name=clean --wrap="rm -rf ${PWD}/${GFF_FOLDER}/UTRme_no_core_A_results")
	jid9=$(echo $jid9 | sed "s/Submitted batch job //")


	jid10=$(sbatch --job-name=UTRme --mem=60G --dependency=afterok:${jid9} --wrap="python ${HOME}/UTRme/utrme.py $INPUT_FOLDER/Configuration_no_core_A_file.txt")
	jid10=$(echo $jid10 | sed "s/Submitted batch job //")

}

add_UTR_annotation(){

GFF_IN=$(echo temp_${ANNOTATION_10_FULLY_PHASED_FILE}_fixed_BES2.gff3 | sed 's/Tb427v10/Tb427v11/g' | sed 's/\.gff3$/\.3\.gff3/')
GFF_IN_NO_CORE_B=${GFF_IN%.gff3}_no_core_B.gff3
GFF_IN_NO_CORE_A=${GFF_IN%.gff3}_no_core_A.gff3

UTRME_5_ANNOTATION_FILE=UTRme-Run-5-best-score.gff
UTRME_3_ANNOTATION_FILE=UTRme-Run-3-best-score.gff


# Using the Rscript developed by Benedikt Brink to integrate UTRme-derived annotation modyfing gene and mRNA entries
# for no_core_B:
UTRME_FOLDER=UTRme_no_core_B_results
GFF_IN=${GFF_IN_NO_CORE_B}
SUFFIX=$(echo $UTRME_FOLDER | sed 's/_no_core_B_results//')
GFF_A_OUT=${GFF_IN%.gff3}_${SUFFIX}.gff3

jid1=$(sbatch --job-name=int-UTRs --wrap="Rscript $SCRIPTS_FOLDER/GFF_fixes_simplified.R ${GFF_FOLDER}/$UTRME_FOLDER/GFF/${UTRME_5_ANNOTATION_FILE%.gff}.gff ${GFF_FOLDER}/$UTRME_FOLDER/GFF/${UTRME_3_ANNOTATION_FILE%.gff}.gff $GFF_FOLDER/$GFF_IN $GFF_FOLDER/${GFF_A_OUT}")
jid1=$(echo $jid1 | sed "s/Submitted batch job //")


jid2=$(sbatch --job-name=add-UTR --dependency=afterok:${jid1} --export=INPUT_GFF_FILE=${GFF_FOLDER}/${GFF_A_OUT},FIVE_ANNOTATION_FILE=${GFF_FOLDER}/${UTRME_FOLDER}/GFF/${UTRME_5_ANNOTATION_FILE},THREE_ANNOTATION_FILE=${GFF_FOLDER}/${UTRME_FOLDER}/GFF/${UTRME_3_ANNOTATION_FILE},OUTPUT_GFF_FILE=${GFF_FOLDER}/${GFF_A_OUT%.gff3}_UTRann.gff3 ${SCRIPTS_FOLDER}/add_UTR_annotations.sh)
jid2=$(echo $jid2 | sed "s/Submitted batch job //")



# for no_core_A:

UTRME_FOLDER=UTRme_no_core_A_results
GFF_IN=${GFF_IN_NO_CORE_A}
SUFFIX=$(echo $UTRME_FOLDER | sed 's/_no_core_A_results//')
GFF_B_OUT=${GFF_IN%.gff3}_${SUFFIX}.gff3

jid3=$(sbatch --job-name=int-UTRs --wrap="Rscript $SCRIPTS_FOLDER/GFF_fixes_simplified.R ${GFF_FOLDER}/$UTRME_FOLDER/GFF/${UTRME_5_ANNOTATION_FILE%.gff}.gff ${GFF_FOLDER}/$UTRME_FOLDER/GFF/${UTRME_3_ANNOTATION_FILE%.gff}.gff $GFF_FOLDER/$GFF_IN $GFF_FOLDER/${GFF_B_OUT}")
jid3=$(echo $jid3 | sed "s/Submitted batch job //")

jid4=$(sbatch --job-name=add-UTR --dependency=afterok:${jid3} --export=INPUT_GFF_FILE=${GFF_FOLDER}/${GFF_B_OUT},FIVE_ANNOTATION_FILE=${GFF_FOLDER}/${UTRME_FOLDER}/GFF/${UTRME_5_ANNOTATION_FILE},THREE_ANNOTATION_FILE=${GFF_FOLDER}/${UTRME_FOLDER}/GFF/${UTRME_3_ANNOTATION_FILE},OUTPUT_GFF_FILE=${GFF_FOLDER}/${GFF_B_OUT%.gff3}_UTRann.gff3 ${SCRIPTS_FOLDER}/add_UTR_annotations.sh)
jid4=$(echo $jid4 | sed "s/Submitted batch job //")




##fix the problem with start 1nt to early for "three_prime_UTR" annotations of "+" strand for both annotation files
#
#jid3=$(sbatch --job-name=fix3UTR --dependency=afterok:${jid1} --export=GFF_IN=${GFF_A_OUT},INPUT_FOLDER=${GFF_FOLDER},OUTPUT_FOLDER=${GFF_FOLDER},GFF_OUT=${GFF_A_OUT%.gff3}_adj3UTR.gff3 ${SCRIPTS_FOLDER}/adjUTRme_three_prime_UTR.sh)
#jid3=$(echo $jid3 | sed "s/Submitted batch job //")
#
#jid4=$(sbatch --job-name=fix3UTR --dependency=afterok:${jid2} --export=GFF_IN=${GFF_B_OUT},INPUT_FOLDER=${GFF_FOLDER},OUTPUT_FOLDER=${GFF_FOLDER},GFF_OUT=${GFF_B_OUT%.gff3}_adj3UTR.gff3 ${SCRIPTS_FOLDER}/adjUTRme_three_prime_UTR.sh)
#jid4=$(echo $jid4 | sed "s/Submitted batch job //")
#
#
#change ID=utr_Tb427 to ID=5utr_Tb427 OR ID=3utr_Tb427 in the gff for every genome version

jid5=$(sbatch --job-name=modUTRID --dependency=afterok:${jid2} --wrap="sed -i '/five_prime_UTR/{s!ID=utr_Tb427_!ID=5utr_Tb427_!g}' ${GFF_FOLDER}/${GFF_A_OUT%.gff3}_UTRann.gff3")
jid5=$(echo $jid5 | sed "s/Submitted batch job //")

jid6=$(sbatch --job-name=modUTRID --dependency=afterok:${jid4} --wrap="sed -i '/five_prime_UTR/{s!ID=utr_Tb427_!ID=5utr_Tb427_!g}' ${GFF_FOLDER}/${GFF_B_OUT%.gff3}_UTRann.gff3")
jid6=$(echo $jid6 | sed "s/Submitted batch job //")

jid7=$(sbatch --job-name=modUTRID --dependency=afterok:${jid5} --wrap="sed -i '/three_prime_UTR/{s!ID=utr_Tb427_!ID=3utr_Tb427_!g}' ${GFF_FOLDER}/${GFF_A_OUT%.gff3}_UTRann.gff3")
jid7=$(echo $jid7 | sed "s/Submitted batch job //")

jid8=$(sbatch --job-name=modUTRID --dependency=afterok:${jid6} --wrap="sed -i '/three_prime_UTR/{s!ID=utr_Tb427_!ID=3utr_Tb427_!g}' ${GFF_FOLDER}/${GFF_B_OUT%.gff3}_UTRann.gff3")
jid8=$(echo $jid8 | sed "s/Submitted batch job //")

## temporary fix to remove the additional ";" at the end of the annotations that is appearing after runinng GFF_fixes.R

jid9=$(sbatch --job-name=rm_semicolon --dependency=afterok:${jid7} --wrap="sed -i 's/;$//g' ${GFF_FOLDER}/${GFF_A_OUT%.gff3}_UTRann.gff3")
jid9=$(echo $jid9 | sed "s/Submitted batch job //")

jid10=$(sbatch --job-name=rm_semicolon --dependency=afterok:${jid8} --wrap="sed -i 's/;$//g' ${GFF_FOLDER}/${GFF_B_OUT%.gff3}_UTRann.gff3")
jid10=$(echo $jid10 | sed "s/Submitted batch job //")

#join together core_B annotations with full_core_A
GFF_OUT=temp_Tb427v11_all.gff3
jid11=$(sbatch --job-name=join --dependency=afterok:${jid9}:${jid10} --wrap="grep '_Tb427v11_B' ${GFF_FOLDER}/${GFF_B_OUT%.gff3}_UTRann.gff3 | cat ${GFF_FOLDER}/${GFF_A_OUT%.gff3}_UTRann.gff3 - > ${GFF_FOLDER}/${GFF_OUT}")
jid11=$(echo $jid11 | sed "s/Submitted batch job //")



}
# once the addition of the UTR annotation is done for each chromosome, now we need to keep the no_core_B.gff as the one it will go ahead as standalone, then rejoin the UTR-annotated core_B part from no_core_A.gff into no_core_B.gff, generating full_UTR_annotated.gff . Now i need to generate the scaffolds of the 11 chromosomes for each haplotype on the fasta and of the gff following a table

# Generation of haplotype_A fasta: scaffold haplotype A and seperate from unassembled contigs
# Generation of haplotype_A gff:
# Generation of haplotype_B fasta: scaffold haplotype B and separate from unassembled contigs
# Generation of haplotype_B gff:
# Generatio of core_A genome fasta: remove core_B chromosomes
# Generation of core_A genome gff:

scaffold_fasta_hap_chroms(){

#scaffolding of the two haplotypes and splitting into three files : haplotype_A
#perl Extract_and_scaffold_sequences_v4.pl FASTA_IN GUIDE > FASTA_OUT
INPUT_GENOME=${GENOME_FOLDER}/${GENOME_11_FULLY_PHASED_FILE}

jid1=$(sbatch --job-name=scaffold --wrap="perl ${SCRIPTS_FOLDER}/Extract_and_scaffold_sequences_v6.pl $INPUT_GENOME $INPUT_FOLDER/${GUIDE_SCAFFOLDING_FILE} > $GENOME_FOLDER/scaffold.fasta")
jid1=$(echo $jid1 | sed "s/Submitted batch job //")


jid2=$(sbatch --job-name=wrap-scaffold --dependency=afterok:${jid1} --wrap="fastaq to_fasta ${GENOME_FOLDER}/scaffold.fasta ${GENOME_FOLDER}/${GENOME_11_FULLY_PHASED_FILE%.fasta}_scaffolded.fasta") # by default it generates 60 nt wrapped sequences
jid2=$(echo $jid2 | sed "s/Submitted batch job //")


#separate the two haplotypes fasta

OUTPUT_GENOME=${GENOME_FOLDER}/Tb427v11_hapA.fasta
jid3=$(sbatch --job-name=hapA --dependency=afterok:${jid2} --wrap="fastaq filter --regex hapA_Tb427v11 ${GENOME_FOLDER}/${GENOME_11_FULLY_PHASED_FILE%.fasta}_scaffolded.fasta ${OUTPUT_GENOME}")
jid3=$(echo $jid3 | sed "s/Submitted batch job //")

OUTPUT_GENOME=${GENOME_FOLDER}/Tb427v11_hapB.fasta
jid4=$(sbatch --job-name=hapB --dependency=afterok:${jid2} --wrap="fastaq filter --regex hapB_Tb427v11 ${GENOME_FOLDER}/${GENOME_11_FULLY_PHASED_FILE%.fasta}_scaffolded.fasta ${OUTPUT_GENOME}")
jid4=$(echo $jid4 | sed "s/Submitted batch job //")

# generate core_A final genome fasta by removing core_B chromosomes and removing the "_A" from the chromosome cores name

OUTPUT_GENOME=${GENOME_FOLDER}/Tb427v11.fasta
jid5=$(sbatch --job-name=coreA  --wrap="fastaq filter --regex _Tb427v11_B -v ${INPUT_GENOME} ${OUTPUT_GENOME%.fasta}_temp.fasta")
jid5=$(echo $jid5 | sed "s/Submitted batch job //")

jid6=$(sbatch --job-name=rm_A --dependency=afterok:${jid5} --wrap="sed -i 's/_core_Tb427v11_A/_core_Tb427v11/g' ${OUTPUT_GENOME%.fasta}_temp.fasta")
jid6=$(echo $jid6 | sed "s/Submitted batch job //")

jid7=$(sbatch --job-name=sortfasta --dependency=afterok:${jid6} --wrap="fastaq sort_by_name ${OUTPUT_GENOME%.fasta}_temp.fasta ${OUTPUT_GENOME}")
jid7=$(echo $jid7 | sed "s/Submitted batch job //")

}

# transfer manual annotation to Tb427v11.fasta

# --> TRANSFER MANUAL ANNOTATION.. ONLY FOR THE UNSCAFFOLD MAIN GENOME?
transfer_manual_annotation(){
#/ladsie/project/ladsie_003/TbAllelePhasingPaperPipelines/04_Genome_annotation_conversion_v9.9_to_v10/HGAP3_Tb427v10_manual.gff3
#/ladsie/project/ladsie_003/TbAllelePhasingPaperPipelines/04_Genome_annotation_conversion_v9.9_to_v10/HGAP3_Tb427v10.fasta
OLD_GENOME_ANNOTATION=${ANNOTATION_10_MANUAL_FILE%.gff3}_PTU.gff3 #/ladsie/project/ladsie_003/TbAllelePhasingPaperPipelines/04_Genome_annotation_conversion_v9.9_to_v10/HGAP3_Tb427v10_manual.gff3
OLD_GENOME_FILE=${GENOME_10_FILE} #/ladsie/project/ladsie_003/TbAllelePhasingPaperPipelines/04_Genome_annotation_conversion_v9.9_to_v10/HGAP3_Tb427v10.fasta
GENOME_FILE=Tb427v11.fasta
#ANNOTATION_IN="Full_${GENOME_FILE%.fasta}.gff3"
OLD_PTU_ANNOTATION_FILE=${ANNOTATION_10_PTU_FILE} #/ladsie/project/ladsie_003/2019-projects/2019-03-11-PTU-annotation/output/PTUs_HGAP3_Tb427v10_sorted_with_IDs.gff3
#ID=Tb427v10_37
# add PTU annotation to the "manual.gff3
jid1=$(sbatch --job-name=addPTUs --wrap="grep -v '^#' ${INPUT_FOLDER}/${ANNOTATION_10_PTU_FILE} | sed 's/ID=Tb427v10/ID=PTU/' | cat ${INPUT_FOLDER}/${ANNOTATION_10_MANUAL_FILE} - > ${GFF_FOLDER}/${OLD_GENOME_ANNOTATION}")
jid1=$(echo $jid1 | sed "s/Submitted batch job //")

OLD_GENOME_FILE_PATH=${INPUT_FOLDER}/${GENOME_10_FILE}
OLD_GENOME_ANNOTATION_PATH=${GFF_FOLDER}/${OLD_GENOME_ANNOTATION}
GENOME_FILE_PATH=${GENOME_FOLDER}/${GENOME_FILE}
#GFF_PATH=${GFF_FOLDER}/temp_Companion_927v41_HGAP3_Tb427v11_phased_diploid_core.gff3_fixed_BES2.3_no_core_B_UTRme.gff3
GFF_PATH=${GFF_FOLDER}/temp_Companion_927v41_HGAP3_Tb427v11_phased_diploid_core.gff3_fixed_BES2.3_no_core_B_UTRme_UTRann.gff3
GFF_NEW_PATH=${GFF_FOLDER}/temp_Companion_927v41_HGAP3_Tb427v11_phased_diploid_core.gff3_fixed_BES2.3_no_core_B_UTRme_UTRann_manual.gff3
#GFF_NEW_PATH=${GFF_FOLDER}/temp_Companion_927v41_HGAP3_Tb427v11_phased_diploid_core.gff3_fixed_BES2.3_no_core_B_UTRme_manual.gff3

jid2=$(sbatch --job-name=transfer-manual-annotation --dependency=afterok:${jid1} --export=OUTPUT_FOLDER=${GFF_FOLDER},OLD_GENOME_ANNOTATION_PATH=${OLD_GENOME_ANNOTATION_PATH},OLD_GENOME_FILE_PATH=${OLD_GENOME_FILE_PATH},SCRIPTS_FOLDER=${SCRIPTS_FOLDER},GENOME_FILE_PATH=${GENOME_FILE_PATH},GFF_PATH=${GFF_PATH},GFF_NEW_PATH=${GFF_NEW_PATH} ${SCRIPTS_FOLDER}/transfer_manual_annotation.sh)
jid2=$(echo $jid2 | sed "s/Submitted batch job //")


}




scaffold_gffs_hap_chroms(){
#echo "under construction\n"
#GFF_IN=temp_Companion_927v41_HGAP3_Tb427v11_phased_diploid_core.gff3_fixed_BES2.3_no_core_B_UTRme.gff3
GFF_IN=temp_Companion_927v41_HGAP3_Tb427v11_phased_diploid_core.gff3_fixed_BES2.3_no_core_B_UTRme_UTRann_manual.gff3
HEADER_IN=temp_Companion_927v41_HGAP3_Tb427v11_phased_diploid_core.gff3_fixed_BES2.3_no_core_B.gff3
GFF_OUT=Tb427v11.gff3
# to get the core_A final genome gff, just remove comment lines of the no_core_B_UTRme.gff3 and add the comment lines from a previous version that had the sizes and so on

jid1=$(sbatch --job-name=gff_header --wrap="grep '^#' ${GFF_FOLDER}/${HEADER_IN} > ${GFF_FOLDER}/header.gff3")
jid1=$(echo $jid1 | sed "s/Submitted batch job //")


jid2=$(sbatch --job-name=final_genome_gff --dependency=afterok:${jid1} --wrap="grep -v '^#' ${GFF_FOLDER}/${GFF_IN} | cat ${GFF_FOLDER}/header.gff3 - > ${GFF_FOLDER}/${GFF_OUT}")
jid2=$(echo $jid2 | sed "s/Submitted batch job //")

# the chromosomes end up with a name "Chr10_core_Tb427v11_A" (for the non-manual annotation), i have to change that to "Chr10_core_Tb427v11"

jid3=$(sbatch --job-name=adj-name-gff --dependency=afterok:${jid2} --wrap="sed -i 's/_core_Tb427v11_A/_core_Tb427v11/g' ${GFF_FOLDER}/${GFF_OUT}")
jid3=$(echo $jid3 | sed "s/Submitted batch job //")



# to get the scaffolds gff 

FASTA_IN=${GENOME_11_FULLY_PHASED_FILE}
GFF_IN=temp_Tb427v11_all.gff3

GFF_HAPS_OUT=Tb427v11_haps_without_header.gff3

jid4=$(sbatch --job-name=scaf_gff  --wrap="perl ${SCRIPTS_FOLDER}/scaffold_gff_v2.pl ${GENOME_FOLDER}/${FASTA_IN} ${GFF_FOLDER}/${GFF_IN} $INPUT_FOLDER/${GUIDE_SCAFFOLDING_FILE} > ${GFF_FOLDER}/${GFF_HAPS_OUT}")
jid4=$(echo $jid4 | sed "s/Submitted batch job //")

#scaffold A
HEADER_IN=header_hapA.gff3
GFF_OUT=Tb427v11_hapA.gff3

#remove hap_B chromosomes from "haps_without_header.gff3" file

jid5=$(sbatch --job-name=hapAgff3 --dependency=afterok:${jid4} --wrap="grep -v 'hapB_Tb427v11' ${GFF_FOLDER}/${GFF_HAPS_OUT} > ${GFF_FOLDER}/${GFF_OUT%.gff3}_without_header.gff3")
jid5=$(echo $jid5 | sed "s/Submitted batch job //")


#Generate a header needed to read chromosome sizes by GENEastics
jid6=$(sbatch --job-name=add_header --wrap="python ${SCRIPTS_FOLDER}/gff_header.py ${GENOME_FOLDER}/Tb427v11_hapA.fasta > ${GFF_FOLDER}/${HEADER_IN}")
jid6=$(echo $jid6 | sed "s/Submitted batch job //")

jid7=$(sbatch --job-name=cat --dependency=afterok:${jid5}:${jid6} --wrap="echo '##gff-version 3' | cat - ${GFF_FOLDER}/${HEADER_IN} ${GFF_FOLDER}/${GFF_OUT%.gff3}_without_header.gff3 > ${GFF_FOLDER}/${GFF_OUT}")
jid7=$(echo $jid7 | sed "s/Submitted batch job //")

#Modify "ID" and "Parent in haplotype "A" by adding an "A"

jid8=$(sbatch --job-name=modID --dependency=afterok:${jid7} --wrap="sed -i '/^Chr/{s!ID=Tb427_!ID=Tb427_A!g}' ${GFF_FOLDER}/${GFF_OUT}")
jid8=$(echo $jid8 | sed "s/Submitted batch job //")

jid9=$(sbatch --job-name=modID --dependency=afterok:${jid8} --wrap="sed -i '/^Chr/{s!Parent=Tb427_!Parent=Tb427_A!g}' ${GFF_FOLDER}/${GFF_OUT}")
jid9=$(echo $jid9 | sed "s/Submitted batch job //")

jid10=$(sbatch --job-name=modID --dependency=afterok:${jid9} --wrap="sed -i '/^Chr/{s!ID=5utr_Tb427_!ID=5utr_Tb427_A!g}' ${GFF_FOLDER}/${GFF_OUT}")
jid10=$(echo $jid10 | sed "s/Submitted batch job //")

jid11=$(sbatch --job-name=modID --dependency=afterok:${jid10} --wrap="sed -i '/^Chr/{s!ID=3utr_Tb427_!ID=3utr_Tb427_A!g}' ${GFF_FOLDER}/${GFF_OUT}")
jid11=$(echo $jid11 | sed "s/Submitted batch job //")





#scaffold B
HEADER_IN=header_hapB.gff3
GFF_OUT=Tb427v11_hapB.gff3

#remove hap_A chromosomes from "haps_without_header.gff3" file

jid12=$(sbatch --job-name=hapAgff3 --dependency=afterok:${jid4} --wrap="grep -v 'hapA_Tb427v11' ${GFF_FOLDER}/${GFF_HAPS_OUT} > ${GFF_FOLDER}/${GFF_OUT%.gff3}_without_header.gff3")
jid12=$(echo $jid12 | sed "s/Submitted batch job //")

#Generate a header needed to read chromosome sizes by GENEastics
jid13=$(sbatch --job-name=add_header --wrap="python ${SCRIPTS_FOLDER}/gff_header.py ${GENOME_FOLDER}/Tb427v11_hapB.fasta > ${GFF_FOLDER}/${HEADER_IN}")
jid13=$(echo $jid13 | sed "s/Submitted batch job //")

jid14=$(sbatch --job-name=cat --dependency=afterok:${jid12}:${jid13} --wrap="echo '##gff-version 3' | cat - ${GFF_FOLDER}/${HEADER_IN} ${GFF_FOLDER}/${GFF_OUT%.gff3}_without_header.gff3 > ${GFF_FOLDER}/${GFF_OUT}")
jid14=$(echo $jid14 | sed "s/Submitted batch job //")

#Modify "ID" and "Parent in haplotype "B" by adding a "B"

jid15=$(sbatch --job-name=modID --dependency=afterok:${jid14} --wrap="sed -i '/^Chr/{s!ID=Tb427_!ID=Tb427_B!g}' ${GFF_FOLDER}/${GFF_OUT}")
jid15=$(echo $jid15 | sed "s/Submitted batch job //")


jid16=$(sbatch --job-name=modID --dependency=afterok:${jid15} --wrap="sed -i '/^Chr/{s!Parent=Tb427_!Parent=Tb427_B!g}' ${GFF_FOLDER}/${GFF_OUT}")
jid16=$(echo $jid16 | sed "s/Submitted batch job //")

jid17=$(sbatch --job-name=modID --dependency=afterok:${jid16} --wrap="sed -i '/^Chr/{s!ID=5utr_Tb427_!ID=5utr_Tb427_B!g}' ${GFF_FOLDER}/${GFF_OUT}")
jid17=$(echo $jid17 | sed "s/Submitted batch job //")

jid18=$(sbatch --job-name=modID --dependency=afterok:${jid17} --wrap="sed -i '/^Chr/{s!ID=3utr_Tb427_!ID=3utr_Tb427_B!g}' ${GFF_FOLDER}/${GFF_OUT}")
jid18=$(echo $jid18 | sed "s/Submitted batch job //")

}

#map RNA-seq data genome to sannity check the UTR annotations

bwa_index_genome(){

module load ngs/bwa/0.7.16

BWA_INDEX_FOLDER=${OUTPUT_FOLDER}/bwa_index
mkdir -p ${BWA_INDEX_FOLDER}
GENOME=Tb427v11.fasta
GENOME_FOLDER=${GENOME_FOLDER}
jid1=$(sbatch --job-name=index --wrap="bwa index -p $BWA_INDEX_FOLDER/$GENOME $GENOME_FOLDER/$GENOME")
jid1=$(echo $jid1 | sed "s/Submitted batch job //")

}




map_RNAseq_to_genome(){

module load ngs/bwa/0.7.16
module load ngs/samtools/1.8

PICARD=~/picard.jar
REMOVE_DUP_METRICS_FOLDER=${OUTPUT_FOLDER}/rmdup_metrics
BWA_ALIGNMENT_FOLDER=${OUTPUT_FOLDER}/bwa_alignments
mkdir -p ${BWA_ALIGNMENT_FOLDER} \
	${REMOVE_DUP_METRICS_FOLDER}

GENOME_FILE=Tb427v11.fasta
BWA_INDEX_FOLDER=${OUTPUT_FOLDER}/bwa_index

for INPUT_FILE_READ_1 in $(ls $READ_ONE_FOLDER | grep "fq.gz$\|fq.bz2$\|fastq.gz$" | grep _R1 )
        do
        PREFIX=${INPUT_FILE_READ_1%_R1.*}
        OUTPUT_BAM=${BWA_ALIGNMENT_FOLDER}/${PREFIX}"_to_"${GENOME_FILE%.fasta}"_bwa.bam"
	INPUT_FILE_READ_2=$(echo $INPUT_FILE_READ_1 | sed "s/_R1/_R2/")
        INPUT_FILE_READ_1=${READ_ONE_FOLDER}/${INPUT_FILE_READ_1}
        INPUT_FILE_READ_2=${READ_TWO_FOLDER}/${INPUT_FILE_READ_2}
        GENOME_INDEX=${BWA_INDEX_FOLDER}/${GENOME_FILE}

        jid2=$(sbatch --partition=slim18 --job-name=alignment --export=GENOME_INDEX=${GENOME_INDEX},INPUT_FILE_READ_1=${INPUT_FILE_READ_1},INPUT_FILE_READ_2=${INPUT_FILE_READ_2},OUTPUT_FILE=${OUTPUT_BAM},PICARD=${PICARD},REMOVE_DUP_METRICS_FOLDER=${REMOVE_DUP_METRICS_FOLDER} ${SCRIPTS_FOLDER}/run_paired_end_alignment.sh)
        jid2=$(echo $jid2 | sed "s/Submitted batch job //")

        done

jid3=$(sbatch --job-name=alignment --dependency=singleton --wrap="echo Alignment step finished.")
jid3=$(echo $jid3 | sed "s/Submitted batch job //")


}


main

