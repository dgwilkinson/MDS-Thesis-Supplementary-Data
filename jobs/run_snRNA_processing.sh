#!/bin/bash
#SBATCH -p shared
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH --gres=tmp:30G
#SBATCH -t 01:00:00
#SBATCH --job-name=run_snRNA_processing
#SBATCH --output=/home/%u/research_project/src/log/%x_%A_%a.out
#SBATCH --error=/home/%u/research_project/src/log/%x_%A_%a.err
#SBATCH --array=0-31
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.g.wilkinson@dur.ac.uk

CUR_DIR=$PWD
cd "/home/${USER}/research_project/dataset/snRNA"
SRC_DIR=../../src/snRNA_analysis
MEDIA_DIR=../../media

SAMPLE_PATHS=($(ls -d [A-Z][A-Z][0-9][0-9][0-9]))

# Check Number of Samples 
# echo ${#SAMPLE_PATHS[@]}
# Print a Specific Sample
# echo ${SAMPLE_PATHS[31]}

SAMPLE=${SAMPLE_PATHS[SLURM_ARRAY_TASK_ID]}

if [[ ! -d "${MEDIA_DIR}/${SAMPLE}" ]]
then
	mkdir "${MEDIA_DIR}/${SAMPLE}"
fi

# Run processing script
Rscript "${SRC_DIR}/process_single_sample.R" \
	--sample_path "${SAMPLE}/outs/filtered_feature_bc_matrix/" \
	--sample_id $SAMPLE \
	--fig_path "${MEDIA_DIR}/${SAMPLE}/" \
	--obj_path "processed_objects/" \
	--dis_path "/home/dwzj28/research_project/markers/dissociation_markers.csv"

cd $CUR_DIR
