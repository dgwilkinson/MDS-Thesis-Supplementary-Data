#!/bin/bash
#SBATCH -p shared
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH --gres=tmp:30G
#SBATCH -t 01:00:00
#SBATCH --job-name=spatial_qc
#SBATCH --output=/home/%u/research_project/src/log/%x_%A_%a.out
#SBATCH --error=/home/%u/research_project/src/log/%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.g.wilkinson@dur.ac.uk

CUR_DIR=$PWD
cd "/home/${USER}/research_project/dataset/spatialRNA"
SRC_DIR=../../src/spatialRNA_analysis
MEDIA_DIR=../../media

#Rscript "${SRC_DIR}/initial_spatial_qc.R" \
	--sample_dir "./" \
	--fig_path "${MEDIA_DIR}/" \

Rscript "${SRC_DIR}/plot_spatial_qc_filtered.R"

cd $CUR_DIR
