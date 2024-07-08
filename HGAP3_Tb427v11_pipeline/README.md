# Correcting some errors in the HGAP3_Tb427v10 genome assembly

## Introduction

This pipeline removes some erroneous sequence in the BES2 contig of the HGAP3_Tb427v10 assembly and correct some other minor mistakes. 
These corrections were made prior to obtaining the ONt reads for the gap closing pipeline.

## Usage

This pipeline is run using the `run.sh` file.

First, you will have to install the 2022-10-11-Generate-updated-haplotype-Lister427-genomes-for-EupathDB_env.yaml (the file is provided). 
Next, you will have to download and configure all necessary input data (paths can be configured in the `run.sh` file, set_variables function). 
- The HGAP3_Tb427v10 input assembly can be found here: https://tritrypdb.org/common/downloads/Current_Release/TbruceiLister427_2018/fasta/data/
- The RNAseq data from [(Krauss et al., 2020), GSE145812](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145812).
- The RNAseq data from [(Jeacock et al., 2018), PRJEB22797](https://www.ebi.ac.uk/ena/browser/view/PRJEB22797).


Finally, you can run the pipeline by executing the `run.sh` file.