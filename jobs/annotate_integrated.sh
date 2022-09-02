#!/bin/bash
#SBATCH -p shared
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH --gres=tmp:30G
#SBATCH -t 01:00:00
#SBATCH --job-name=annotate_integrated
#SBATCH --output=/home/%u/research_project/src/log/%x_%A_%a.out
#SBATCH --error=/home/%u/research_project/src/log/%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.g.wilkinson@dur.ac.uk

CUR_DIR=$PWD
cd "/home/${USER}/research_project/dataset/snRNA/"
SRC_DIR=../../src/snRNA_analysis
MEDIA_DIR=../../media

Rscript "${SRC_DIR}/annotate_obj.R" \
        --sample_path "./integrated_object/cca_integrated_obj_DS.rds" \
        --sample_id "integrated_DS" \
        --fig_path "${MEDIA_DIR}/" \
	--anno_path "${MEDIA_DIR}/integrated_DS_cluster_annotations.csv"
	
cd $CUR_DIR
