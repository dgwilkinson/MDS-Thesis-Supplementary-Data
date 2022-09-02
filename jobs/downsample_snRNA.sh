#!/bin/bash
#SBATCH -p shared
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH --gres=tmp:30G
#SBATCH -t 01:00:00
#SBATCH --job-name=downsample_snRNA
#SBATCH --output=/home/%u/research_project/src/log/%x_%A_%a.out
#SBATCH --error=/home/%u/research_project/src/log/%x_%A_%a.err
#SBATCH --array=0-31
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.g.wilkinson@dur.ac.uk

CUR_DIR=$PWD
cd "/home/${USER}/research_project/dataset/snRNA/"
SRC_DIR=../../src/snRNA_analysis

SAMPLE_PATHS=($(ls -d [A-Z][A-Z][0-9][0-9][0-9]))

SAMPLE=${SAMPLE_PATHS[SLURM_ARRAY_TASK_ID]}

Rscript "${SRC_DIR}/downsample_snRNA.R" \
        --sample_path "./processed_objects/${SAMPLE}.rds" \

cd $CUR_DIR
