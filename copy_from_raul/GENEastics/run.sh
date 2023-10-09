#!/bin/bash
#SBATCH --partition=slim18              # (queue, see sinfo) 
#SBATCH --ntasks=16                     # number of threads
#SBATCH --mem 128G                      #memory pool for all cores 
#SBATCH -o slurm.%u.%j.out              #STDOUT 
#SBATCH -e slurm.%u.%j.err              #STDERR


# This module aims to generate a GENEastics plot of the fully phased genome with WT RNA-seq data on top

#source /work/project/ladsie_003/GENEastics_env/bin/activate
module load Anaconda3/2019.03

source activate GENEastics_env

main(){
	set_variables
	create_folders
	link_genomes
	link_wigs
	plot_RNAseq_full_diploid_scaffold

}

set_variables(){
        INPUT_FOLDER=input
        OUTPUT_FOLDER=output
        SCRIPTS_FOLDER=bin
	GENOME_FOLDER=${INPUT_FOLDER}/genomes
	GENOME_PATH=/work/project/ladsie_001/Tbrucei_Genomes/HGAP3_Tb427v10_diploid
	DIPLOID_SCAFFOLDED_ANNOTATION_FILE=HGAP3_Tb427v10_diploid_scaffolded.gff3
	FULL_PLOT_FOLDER=${OUTPUT_FOLDER}/Full_scaffold_plots

}

create_folders(){
        mkdir -p \
        $SCRIPTS_FOLDER \
        $INPUT_FOLDER \
	$GENOME_FOLDER \
        $OUTPUT_FOLDER \
	$FULL_PLOT_FOLDER

}



link_genomes(){

	ln -sf $GENOME_PATH/$DIPLOID_SCAFFOLDED_ANNOTATION_FILE $GENOME_FOLDER

}

link_wigs(){

WIG_PATH=/work/project/ladsie_003/2020-06-09-Mapping-datasets-on-fully-phased-scaffold/output/wig_folder/ws10001_ss1000

	ln -sf $WIG_PATH/WT_merged_over_WT_merged_ws10001ss1000_denominator_WT_merged_to_HGAP3_Tb427v10_diploid_scaffolded_bwa.rmdup.sorted.wig $INPUT_FOLDER/WT_RNAseq_on_HGAP3_Tb427v10_diploid_scaffolded.wig

}


plot_RNAseq_full_diploid_scaffold(){

#light gray = #DADEDF
#light yellow = #FFFFB7
#product=term%3DTrypanosomal VSG domain containing protein
#Name=Similar to Tb427VSG
#'product|VSG domain containing protein|#B833FF'
#signature_desc=Trypanosomal VSG domain
#'product|VSG domain|#FF0000|||ID|Boing|Blue'
#Name=Similar to Tb427VSG
#product=term%3DTrypanosomal VSG domain containing protein

python ${SCRIPTS_FOLDER}/geneastics.py \
	--gff_file ${GENOME_FOLDER}/${DIPLOID_SCAFFOLDED_ANNOTATION_FILE} \
	--replicons "Chr11_A_Tb427v10" "Chr11_B_Tb427v10" "Chr10_A_Tb427v10" "Chr10_B_Tb427v10" "Chr9_A_Tb427v10" "Chr9_B_Tb427v10" "Chr8_A_Tb427v10" "Chr8_B_Tb427v10" "Chr7_A_Tb427v10" "Chr7_B_Tb427v10" "Chr6_A_Tb427v10" "Chr6_B_Tb427v10" "Chr5_A_Tb427v10" "Chr5_B_Tb427v10" "Chr4_A_Tb427v10" "Chr4_B_Tb427v10" "Chr3_A_Tb427v10" "Chr3_B_Tb427v10" "Chr2_A_Tb427v10" "Chr2_B_Tb427v10" "Chr1_A_Tb427v10" "Chr1_B_Tb427v10" \
	--feature_types "polypeptide" "mRNA" "pseudogene" "rRNA" "protein_match" \
	--alpha 0.99 \
	--feature_color_mapping "polypeptide=None;mRNA=#DADEDF;rRNA=green;protein_match=None" \
	--attribute_color_mapping 'signature_desc|Trypanosomal VSG domain|Red|||Name|Similar to Tb427VSG|Red|||product|Trypanosomal VSG|Red' \
	--x_tick_distance 500000 \
	--output_file ${FULL_PLOT_FOLDER}/Full_scaffold_RNAseq_GENEastics.svg \
	--wiggle_files ${INPUT_FOLDER}/WT_RNAseq_on_HGAP3_Tb427v10_diploid_scaffolded.wig \
	--wiggle_y_max 50

}


main

