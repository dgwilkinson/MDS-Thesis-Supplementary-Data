# MDS-Thesis-Supplementary-Data

Supplementary Data for my Master of Data Science Research Project entitled *"Implementing a local analysis pipeline for the integration of single-cell and spatial transcriptomics data"*

### Raw Visualisations + Tables:
- `Alternative_snRNA-seq_Integrations` - A collection of plots depicting various procedures for the integration of snRNA-seq data incl. Subsetting samples (`P14` suffix), RPCA-based on downsampled population, Referenced CCA-based (`ref` suffix)
- `Final_snRNA-seq_Integration` - Plots from downstream analysis on the CCA-based Downsampled snRNA-seq Integration
- `Region_Death_Signature_and_Loss` - Plots depicting trends in spot loss and death signature between major regions (IZ, FZ, BZ, RZ, and Ctrl)
- `Sample QC` - Quality control plots and tables across snRNA-seq and spatial-seq samples before and after processing 
- `Spatial-seq_Sample_Visualisations` - Visualisations produced throughout the pipeline for individual spatial-seq samples 
- `snRNA-seq_Sample_Visualisations` - Visualisations produced throughout the pipeline for individual snRNA-seq samples
- `test_figures` - Visualisations generated in preliminary execution of pipeline steps

### Analysis Pipeline Source Code:
- `jobs` - Bash scripts used to run scripts on the SLURM workload array of Hamilton 8
- `log` - Log files output from the execution of scripts on the SLURM workload array of Hamilton 8
- `snRNA_analysis` - Scripts coordinating processing, integration, and annotation of snRNA-seq samples
- `spatialRNA_analysis` - Scripts coordinating processing, death signature detection, and cell type deconvolution on spatial-seq samples
- `ArchR_test.R` - Very cursory investigation of ArchR for snATAC-seq processing (deemed out of scope)
- `download_dataset.sh` - Bash script for downloading the entire dataset from a `.txt` file of links
- `get_download_links.py` - Retrieves download links for files in a Zenodo repository
- `optimise_clustering.R` - Contains a cluster optimisation function used by both snRNA-seq and spatial-seq analysis

***
**Note**: The implementation of this pipeline was inspired by the publication ["*Spatial multi-omic map of human myocardial infarction*"](https://www.nature.com/articles/s41586-022-05060-x) and greatly facilitated by guidance from code available in the [saezlab/visium_heart](https://github.com/saezlab/visium_heart) GitHub repository. 
***
