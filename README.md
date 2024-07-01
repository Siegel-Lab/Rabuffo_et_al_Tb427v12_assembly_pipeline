# Closing gaps in the T. brucei Lister 427 genome assembly using ONT reads

## Introduction

This repository contains the scripts used to improve the T. brucei Lister 427 genome assembly using ONT long reads. Roughly, the pipeline does the following: It runs SAMBA multiple times to close gaps in the assembly. Next, we also expand collapsed repeats by turning them into gaps and closing those gaps using SAMBA. We annotate the newly created sequences and do quality control on them.

## Usage

This pipeline is run using the `run.sh` file.
First, you will have to install the ont_assembly & ont_assembly_2 conda environments (the .yaml files are provided). Then you will have to install MaSucRCA (version 4.1.0, https://github.com/alekseyzimin/masurca) and seqtk (version 1.4, https://github.com/lh3/seqtk) from their githubs into the bin folder. Next, you will have to download and configure all necessary input data (paths can be configured at the top of the `run.sh` file, see list below). Finally, you can run the pipeline by executing the `run.sh` file.


You will need the following data:
- ONT long reads: @todo link
- the HGAP3_Tb427v11 genome assembly & annotation: @todo link, this version of the assembly is generated from the HGAP3_Tb427v10 assembly using the pipeline in the HGAP3_Tb427v11_pipeleine folder


## Citation

If you use the HGAP3_Tb427v12 genome assembly or the scripts in this repository in your research, please cite the following publication:

@todo fill in once manuscript is public