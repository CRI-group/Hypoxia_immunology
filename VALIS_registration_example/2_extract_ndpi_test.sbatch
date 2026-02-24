#!/bin/bash
#
#SBATCH -J extract_ndpi
#
#SBATCH --output=logs/extract_ndpi_%j.txt
#SBATCH --error=logs/extract_ndpi_%j.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#
# time format hours:minutes:seconds or days-hours
#SBATCH --time=3:59:00
#SBATCH --mem=50G
#SBATCH --partition=small

# (no updated Tcl for 13.3.0)
module load ImageMagick/7.1.1-34-GCCcore-13.2.0

#This scripts calls the tool ndpisplit to extract x20, x10, x5 or x2.5 .tif images from raw Hamamatsu .ndpi files

#set the raw image location
raw_images='/scratch/svc_td_cri/projects/multiplex/hypoxia/IDHwt5/ndpi_images'

#set the output location (untested, can also be put in command)
output_loc='/scratch/svc_td_cri/projects/multiplex/hypoxia/IDHwt5/x20_images'

for image in $(find $raw_images -type f); do
	echo $image 
	# -cn NO COMPRESSION -cl LZW compression
	/scratch/svc_td_cri/projects/multiplex/packages/ndpitools/ndpisplit -cl -x20 -O $output_loc  "$image"
done