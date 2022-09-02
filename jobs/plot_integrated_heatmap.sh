#!/bin/bash
#SBATCH -p bigmem
#SBATCH -c 32
#SBATCH --mem=356G
#SBATCH --gres=tmp:96G
#SBATCH -t 01:00:00
#SBATCH --job-name=make_integrated_heatmap
#SBATCH --output=/home/%u/research_project/src/log/%x_%A_%a.out
#SBATCH --error=/home/%u/research_project/src/log/%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.g.wilkinson@dur.ac.uk

CUR_DIR=$PWD
cd "/home/${USER}/research_project/dataset/snRNA/"
SRC_DIR=../../src/snRNA_analysis
MEDIA_DIR=../../media

Rscript "${SRC_DIR}/marker_heatmap.R" \
        --sample_path "./integrated_object/integrated_snRNA.rds" \
        --sample_id "integrated" \
        --fig_path "${MEDIA_DIR}/" \
        --marker_path "/home/dwzj28/research_project/markers/paper_cell_type_markers.rds" \
        --marker_id "paper"

cd $CUR_DIR
