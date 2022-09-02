#!/bin/bash
#SBATCH -p shared
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH --gres=tmp:30G
#SBATCH -t 01:00:00
#SBATCH --job-name=get_death_sig
#SBATCH --output=/home/%u/research_project/src/log/%x_%A_%a.out
#SBATCH --error=/home/%u/research_project/src/log/%x_%A_%a.err
#SBATCH --array=0-27
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.g.wilkinson@dur.ac.uk

CUR_DIR=$PWD
cd "/home/${USER}/research_project/dataset/spatialRNA"
SRC_DIR=../../src/spatialRNA_analysis
MEDIA_DIR=../../media

# Create sample names array 
SAMPLE_PATHS=($(cat ./sample_names.txt))

# Set sample using SLURM_ARRAY_TASK_ID
SAMPLE=${SAMPLE_PATHS[SLURM_ARRAY_TASK_ID]}

Rscript "${SRC_DIR}/death_sig.R" \
	--sample_path "${SAMPLE}/outs/" \
	--sample_id $SAMPLE \
	--fig_path "${MEDIA_DIR}/${SAMPLE}/" \
	--marker_path "../../markers/death_signature/" 

cd $CUR_DIR
