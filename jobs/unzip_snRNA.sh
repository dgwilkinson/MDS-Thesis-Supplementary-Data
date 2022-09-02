#!/bin/bash

CUR_DIR=$PWD
cd "/home/${USER}/research_project"

SAMPLES=$(ls ./dataset/snRNA/*.tar.gz)

for sample in $SAMPLES
do
	echo "Unzipping:" $sample
	tar -xf $sample -C ./dataset/snRNA/
done

ls ./dataset/snRNA/!(*.tar.gz|sample_data|*.xlsx)

cd $CUR_DIR
