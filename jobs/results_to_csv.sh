#!/bin/bash
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --gres=tmp:3G
#SBATCH -t 00:10:00
#SBATCH --job-name=res_to_csv
#SBATCH --output=/home/%u/research_project/src/log/%x_%A_%a.out
#SBATCH --error=/home/%u/research_project/src/log/%x_%A_%a.err
#SBATCH --array=0-31
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.g.wilkinson@dur.ac.uk

CUR_DIR=$PWD
cd "/home/${USER}/research_project/dataset/spatialRNA"
SRC_DIR=../../src/spatialRNA_analysis
MEDIA_DIR=../../media

# Set miniconda environment
source activate /home/dwzj28/miniconda3/envs/cell2loc_env
export PATH=/home/dwzj28/miniconda3/envs/cell2loc_env/bin:${PATH}

# Create sample names array 
SAMPLE_PATHS=($(cat ./sample_names.txt))

# Set sample using SLURM_ARRAY_TASK_ID
SAMPLE=${SAMPLE_PATHS[SLURM_ARRAY_TASK_ID]}

python "${SRC_DIR}/results_to_csv.py" \
	"deconv_results/${SAMPLE}_abundances.h5ad" \
	"deconv_results/abd_csvs/" \
	${SAMPLE}

cd $CUR_DIR
