#!/bin/bash
#SBATCH -p shared
#SBATCH -c 10
#SBATCH --mem=6G
#SBATCH --gres=tmp:1G
#SBATCH -t 00:01:00
#SBATCH --job-name=check_threads
#SBATCH --output=/home/%u/research_project/src/log/%x_%A_%a.out
#SBATCH --error=/home/%u/research_project/src/log/%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.g.wilkinson@dur.ac.uk

CUR_DIR=$PWD
cd "/home/${USER}/research_project/dataset/spatialRNA"
SRC_DIR=../../src/spatialRNA_analysis
MEDIA_DIR=../../media

python "${SRC_DIR}/check_threads.py" \

cd $CUR_DIR
