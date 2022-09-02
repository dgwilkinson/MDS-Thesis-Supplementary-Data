#!/bin/bash
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --gres=tmp:1G
#SBATCH -t 00:10:00
#SBATCH --job-name=regional_death_sig
#SBATCH --output=/home/%u/research_project/src/log/%x_%A_%a.out
#SBATCH --error=/home/%u/research_project/src/log/%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.g.wilkinson@dur.ac.uk

CUR_DIR=$PWD
cd "/home/${USER}/research_project/dataset/spatialRNA/"
SRC_DIR=../../src/spatialRNA_analysis
MEDIA_DIR=../../media

Rscript "${SRC_DIR}/plot_DS_distribution.R" \
	--sample_ids "sample_names.txt" \
        --fig_path "${MEDIA_DIR}/" \
	
cd $CUR_DIR
