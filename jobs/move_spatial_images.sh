#!/bin/bash

CUR_DIR=$PWD
cd "/home/${USER}/research_project/dataset/spatialRNA"

SAMPLES=$(ls -d !(*tar.gz|*csv|sample_data))

for sample in $SAMPLES
do
        #Unzip the spatial filtered feature data
        echo "Processing: $sample"
        cd ${sample}/outs
	if [[ ! ( -d spatial || -d Volumes ) ]]
        then
                unzip "spatial_${sample}_spatial.zip"
        fi

        #Move the filtered_feature_bc_matrix to a more sensible filepath
        if [[ ! -d spatial ]]
        then
                #Will throw an error claiming file cannot be found because it was just moved
		#mindepth necessary as two dirs are named spatial (we want last one)
                find Volumes -mindepth 7 -type d -name '*spatial*' -exec mv {}/ . \;
                rm -rf Volumes
        fi
	
	# Remove unnecessary files
	rm -rf filtered_feature_bc_matrix

	cd ../../
done
cd $CUR_DIR
