#!/bin/bash
#SBATCH -p shared
#SBATCH -c 32
#SBATCH -N 1
#SBATCH --mem=24G
#SBATCH --gres=tmp:30G
#SBATCH -t 05:00:00
#SBATCH --job-name=run_nb_model
#SBATCH --output=/home/%u/research_project/src/log/%x_%A_%a.out
#SBATCH --error=/home/%u/research_project/src/log/%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.g.wilkinson@dur.ac.uk

CUR_DIR=$PWD
cd "/home/${USER}/research_project/dataset/"
SRC_DIR=../src/spatialRNA_analysis

# Set miniconda environment 
source activate /home/dwzj28/miniconda3/envs/cell2loc_env
export PATH=/home/dwzj28/miniconda3/envs/cell2loc_env/bin:${PATH}

python "${SRC_DIR}/run_nb_model.py" \
        "snRNA/integrated_object/cca_integrated_obj_DS.h5ad" \
        "spatialRNA/deconv_results/" \
        "integrated_DS" \
	
cd $CUR_DIR
