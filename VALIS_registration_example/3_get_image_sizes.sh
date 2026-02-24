#!/bin/bash
#This script gets the image sizes of the given directory, easier to use if you pipe it with > yourfile.log

#load this in session to run these faster (if checking image sizes of every smaple
module load ImageMagick/7.1.1-34-GCCcore-13.2.0

#sample_dir="$1"
image_loc="/scratch/svc_td_cri/projects/multiplex/senesence"
#image_loc=$sample_dir

min_height=1000000000000000
min_width=1000000000000000

#for image in $(find $image_loc -type f); do
for image in $image_loc/*/x20_images/*; do
	#identify $image

	out=$(identify "$image")
	res=$(echo $out | cut -d' ' -f3)

	height=$(echo "$res" | cut -d'x' -f2)
    	width=$(echo "$res" | cut -d'x' -f1)

    	if (( height < min_height )); then
        	min_height=$height
    	fi

    	if (( width < min_width )); then
        	min_width=$width
    	fi
	echo $out
done

echo "Smallest resolution: ${min_width}x$min_height"
