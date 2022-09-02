#!/bin/bash
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --gres=tmp:1G
#SBATCH -t 00:10:00
#SBATCH --job-name=extract_inf_aver
#SBATCH --output=/home/%u/research_project/src/log/%x_%A_%a.out
#SBATCH --error=/home/%u/research_project/src/log/%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.g.wilkinson@dur.ac.uk

CUR_DIR=$PWD
cd "/home/${USER}/research_project/dataset/spatialRNA"
SRC_DIR=../../src/spatialRNA_analysis

# Set miniconda environment 
source activate /home/dwzj28/miniconda3/envs/cell2loc_env
export PATH=/home/dwzj28/miniconda3/envs/cell2loc_env/bin:${PATH}

python "${SRC_DIR}/extract_inf_aver.py" \
        "deconv_results/" \
        "integrated_DS" \
        "deconv_results/" \
	
cd $CUR_DIR
