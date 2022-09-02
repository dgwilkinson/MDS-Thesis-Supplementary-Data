#!/bin/bash
#SBATCH -p shared
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --gres=tmp:30G
#SBATCH -t 04:00:00
#SBATCH --job-name=find_integrated_markers
#SBATCH --output=/home/%u/research_project/src/log/%x_%A_%a.out
#SBATCH --error=/home/%u/research_project/src/log/%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.g.wilkinson@dur.ac.uk

CUR_DIR=$PWD
cd "/home/${USER}/research_project/dataset/snRNA/"
SRC_DIR=../../src/snRNA_analysis
MEDIA_DIR=../../media

Rscript "${SRC_DIR}/integrated_dea.R" \
        --sample_path "./integrated_object/cca_integrated_obj_DS.rds" \
        --sample_id "integrated_DS" \
        --marker_path "~/research_project/markers/" \
	
cd $CUR_DIR
