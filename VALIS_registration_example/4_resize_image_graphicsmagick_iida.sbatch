#!/bin/bash
#
#SBATCH -J resize_image
#
#SBATCH --output=logs/resize_%j.txt
#SBATCH --error=logs/resize_%j.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#
# time format hours:minutes:seconds or days-hours
#SBATCH --time=3:59:00
#SBATCH --mem=135G
#SBATCH --partition=small

module load ImageMagick/7.1.1-34-GCCcore-13.2.0

#This script gets the image sizes of the given directory and crops them to be the same size, 
#fills the edges with black if the smallest image is smaller than the cropped size

image_loc="/scratch/svc_td_cri/projects/multiplex/hypoxia/IDHwt5/x20_images"
output_directory="/scratch/svc_td_cri/projects/multiplex/hypoxia/IDHwt5/raw_images_resized/"

for image in $(find $image_loc -type f); do
	image_output_dir="${image}"
	image_output_dir="$output_directory/$(basename $image_output_dir)"
	echo $image
	magick $image -compress lzw -colorspace GRAY -gravity Center -crop 55552x26880+0+0 +repage $image_output_dir
done

