#!/bin/bash
#SBATCH -p shared
#SBATCH -n 10
#SBATCH --mem=20G
#SBATCH --gres=tmp:30G
#SBATCH -t 03:00:00
#SBATCH --job-name=test_cores
#SBATCH --output=/home/%u/research_project/src/log/%x_%A_%a.out
#SBATCH --error=/home/%u/research_project/src/log/%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.g.wilkinson@dur.ac.uk

export R_PARALLELLY_AVAILABLECORES_METHODS="system,mc.cores,Slurm,fallback,custom"

echo $OMP_NUM_THREADS

Rscript /home/dwzj28/research_project/src/snRNA_analysis/test_available_cores.R

