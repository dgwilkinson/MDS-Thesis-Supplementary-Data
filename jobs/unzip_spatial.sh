#!/bin/bash

CUR_DIR=$PWD
cd "/home/${USER}/research_project"

SAMPLES=$(ls ./dataset/spatialRNA/*.tar.gz)

for sample in $SAMPLES
do
        echo "Unzipping:" $sample
        tar -xf $sample -C ./dataset/spatialRNA/
done

ls ./dataset/spatialRNA/!(*.tar.gz|sample_data|*.csv)

cd $CUR_DIR
