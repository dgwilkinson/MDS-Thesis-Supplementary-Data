#!/bin/bash
#SBATCH -p shared
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH --gres=tmp:30G
#SBATCH -t 00:10:00
#SBATCH --job-name=generate_integrated_annotations
#SBATCH --output=/home/%u/research_project/src/log/%x_%A_%a.out
#SBATCH --error=/home/%u/research_project/src/log/%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.g.wilkinson@dur.ac.uk

CUR_DIR=$PWD
cd "/home/${USER}/research_project/markers"
SRC_DIR=../src/snRNA_analysis
MEDIA_DIR=../media

Rscript "${SRC_DIR}/generate_annotations.R" \
        --sample_marker_path "./integrated_DS_cluster_markers.csv" \
        --sample_id "integrated_DS" \
        --fig_path "${MEDIA_DIR}/" \
	--anno_marker_path "./hca_ct_markers.csv" \
	
cd $CUR_DIR
