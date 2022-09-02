#!/bin/bash
#SBATCH -p shared
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH --gres=tmp:48G
#SBATCH -t 01:00:00
#SBATCH --job-name=integrate_snRNA_DS
#SBATCH --output=/home/%u/research_project/src/log/%x_%A_%a.out
#SBATCH --error=/home/%u/research_project/src/log/%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.g.wilkinson@dur.ac.uk

CUR_DIR=$PWD
cd "/home/${USER}/research_project/dataset/snRNA"
SRC_DIR=../../src/snRNA_analysis

export R_PARALLELLY_AVAILABLECORES_METHODS="system,mc.cores,Slurm,fallback,custom"

Rscript "${SRC_DIR}/integrate_snRNA_DS.R" \
        --objects_dir "./processed_objects" \
        --out_path "./integrated_object" \
        --fig_path "/home/dwzj28/research_project/media" \
        --metadata_path "metadata-RNA-seq.xlsx"

cd $CUR_DIR
