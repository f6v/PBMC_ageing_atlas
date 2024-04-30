#### PBMC Ageing Atlas

Source code for the PBMC Ageing Atlas analysis that includes over a million cells from young and old donors and seven studies.

<img width="800" alt="image" src="https://github.com/f6v/PBMC_ageing_atlas/assets/13019221/8f03e688-d20f-4b68-b09c-5b6086b7257f">


##### Data

The datasets can be found under the following accessions on NCBI, GSA, and Synapse: GSE157007, GSE213516, GSE214546, HRA000203, HRA000624, HRA003766, syn22255433.

##### Directory structure

`GSE157007`, `GSE213516`, `GSE214546`, `HRA000203`, `HRA000624`, `HRA003766`, `syn22255433` contain scripts to combine the sample data for each study and perform the quality control.

`atlas` contains code for integrating the seven datasets and performing the downstream analyses.

`T_markers` - visualization of the T cell reference for validating marker genes.

`utils_py` - shared Python utilities.