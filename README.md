#### PBMC Ageing Atlas

Source code for the PBMC Ageing Atlas analysis that includes over a million cells from young and old donors and seven studies.

<img width="750" alt="image" src="https://github.com/f6v/PBMC_ageing_atlas/assets/13019221/55eb7ce5-1772-4b24-8bf7-56e3e977d10a">


### Data acquisition

The datasets can be found under the following accessions on NCBI, GSA, and Synapse: GSE157007, GSE213516, GSE214546, HRA000203, HRA000624, HRA003766, syn22255433.

### Resource requirements

Please note that `scVI` scripts run much faster on GPU. 

### Directory structure

#### Dataset-specific folders

`GSE157007`, `GSE213516`, `GSE214546`, `HRA000203`, `HRA000624`, `HRA003766`, `syn22255433` contain scripts to combine the sample data for each study and perform the quality control. The scripts in each folder should be run in the following order:
1. `create_adatas.py` to import the data into AnnData format.
2. `doublets.R` to perform the doublet calls for each sample.
3. `combine.py` to create a single AnnData file with all the samples and doublet calls for each dataset.
4. `qc.py` to peform the quality control for each dataset.

#### Combined atlas analysis

`atlas` contains code for integrating the seven datasets and performing the downstream analyses. The scripts in each folder should be run in the following order:
1. `combine_datasets.py` to create a single AnnData object for seven datasets.
2. `prepare_combined.py` to peform the clean up, such as removing V(D)J genes.
3. `integrate.py` to run scVI integration (preferably on a GPU).
4. `viz_integrated.py` to perform additional QC (removing doublets and RBC contamination) and PBMC annotation.

`T_integrate.py`, `B_integrate.py`, `MALAT1_integrate.py` perform the T, B, and MALAT1+ cell re-analysis, respectively.

`celltypist_run.py` and `celltypist_viz.py` perform CellTypist classification for T cells.

`scanpy_DE.py` performs DE test between young and old individuals.

`cell_type_props_plots.R` plots of cell type proportions in the datasets.

`harmony_integrate.py`, `harmony_cluster.py` - alternative integration with Harmony.

`abundance_heatmap.R` - cell type proportions heatmap in each sample and dataset.


#### T cell marker gene analysis

`T_markers` - visualization of the T cell reference for validating marker genes.

#### Utility functions

`utils_py` - shared Python utilities.

### Acknowledgements

This study was supported by the funding from the European Union's Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No.: 955321, as well by Estonian Research Council grant PRG2011. This publication is based upon work partially supported by the Google Cloud Research Credits program award No.: GCP19980904.
