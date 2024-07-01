# Closing gaps in the T. brucei Lister 427 genome assembly

## Introduction

## Usage

This pipeline is run using the `run.sh` file.
First, you will have to install the ont_assembly & ont_assembly_2 conda environments (the .yaml files are provided). Then you will have to install MaSucRCA (version 4.1.0, https://github.com/alekseyzimin/masurca) and seqtk (version 1.4, https://github.com/lh3/seqtk) from their githubs into the bin folder. Next, you will have to download and configure all necessary input data (paths can be configured at the top of the `run.sh` file, see list below). Finally, you can run the pipeline by executing the `run.sh` file.


You will need the following data:
- ONT long reads: @todo add link
- the HGAP3_Tb427v11 genome assembly & annotation: https://tritrypdb.org/common/downloads/Current_Release/TbruceiLister427_2018/fasta/data/ @todo this is version 10, but we need version 11


## Citation

If you use the HGAP3_Tb427v12 genome assembly or the scripts in this repository in your research, please cite the following publication:

@todo fill in once manuscript is public