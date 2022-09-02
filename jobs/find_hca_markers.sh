#!/bin/bash
#SBATCH -p shared
#SBATCH -n 16
#SBATCH --mem=128G
#SBATCH --gres=tmp:30G
#SBATCH -t 05:00:00
#SBATCH --job-name=find_hca_markers
#SBATCH --output=/home/%u/research_project/src/log/%x_%A_%a.out
#SBATCH --error=/home/%u/research_project/src/log/%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.g.wilkinson@dur.ac.uk

OMP_NUM_THREADS=${SLURM_CPUS_ON_NODE}

Rscript /home/dwzj28/research_project/src/snRNA_analysis/find_hca_markers.R

