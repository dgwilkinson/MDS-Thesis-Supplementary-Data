#!/bin/bash

data_dir=/nobackup/dwzj28/research_project/dataset

all_links=$(cat "${data_dir}/data_links/"*)

#Subset links for which download failed
sub_links=$(sed -n '55,$p' <<< $all_links)

while read link
do
	filename=${link#* }
	url=${link% *}
	echo "Performing download on..."
	echo "URL: $url"
	echo "FILENAME: $filename"
	curl -o "${data_dir}/${filename}" $url
	echo "OUTPUT: ${data_dir}/${filename}"
done <<< $sub_links

