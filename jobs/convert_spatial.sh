#!/bin/bash
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --gres=tmp:1G
#SBATCH -t 01:00:00
#SBATCH --job-name=convert_spatial_H5AD
#SBATCH --output=/home/%u/research_project/src/log/%x_%A_%a.out
#SBATCH --error=/home/%u/research_project/src/log/%x_%A_%a.err
#SBATCH --array=0-31
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.g.wilkinson@dur.ac.uk

CUR_DIR=$PWD
cd "/home/${USER}/research_project/dataset/spatialRNA"
SRC_DIR=../../src/snRNA_analysis

# Create sample names array 
SAMPLE_PATHS=($(cat ./sample_names.txt))

# Set sample using SLURM_ARRAY_TASK_ID
SAMPLE=${SAMPLE_PATHS[SLURM_ARRAY_TASK_ID]}

Rscript "${SRC_DIR}/convert_H5AD.R" \
	--seurat_path "./processed_objects/" \
	--seurat_name "${SAMPLE}.rds" \

cd $CUR_DIR
