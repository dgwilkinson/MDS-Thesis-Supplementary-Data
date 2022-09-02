#!/bin/bash

CUR_DIR=$PWD
cd "/home/${USER}/research_project/dataset/snRNA"

SAMPLES=$(ls -d !(*tar.gz|*xlsx|sample_data))

for sample in $SAMPLES
do
	#Unzip the snRNA filtered feature data
	echo "Processing: $sample"
	cd ${sample}/outs
	if [[ ! ( -d Volumes || -d filtered_feature_bc_matrix ) ]]
	then
		unzip "snRNA_${sample}_filtered_feature_bc_matrix.zip"
	fi
	
	#Move the filtered_feature_bc_matrix to a more sensible filepath
	if [[ ! -d filtered_feature_bc_matrix ]]
	then 
		#Will throw an error claiming file cannot be found because it was just moved
		find Volumes -type d -name '*filtered*' -exec mv {}/ . \;
		rm -rf Volumes
	fi

	cd ../../
done

cd $CUR_DIR
