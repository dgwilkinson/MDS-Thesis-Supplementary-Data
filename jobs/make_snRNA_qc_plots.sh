#!/bin/bash
#SBATCH -p shared
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH --gres=tmp:30G
#SBATCH -t 01:00:00
#SBATCH --job-name=plot_snRNA_qcs
#SBATCH --output=/home/%u/research_project/src/log/%x_%A_%a.out
#SBATCH --error=/home/%u/research_project/src/log/%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.g.wilkinson@dur.ac.uk


CUR_DIR=$PWD

cd "/home/${USER}/research_project/"

Rscript ./src/snRNA_analysis/plot_snRNA_qc.R \
	--data_path "/home/dwzj28/research_project/dataset/snRNA" \
	--fig_path "/home/dwzj28/research_project/media/"

# Processes data_path and fig_path within the Rscript
Rscript ./src/snRNA_analysis/plot_snRNA_qc_filtered.R

cd $CUR_DIR
