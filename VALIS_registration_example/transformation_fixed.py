# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 13:22:23 2023

@author: Iida
"""
import numpy as np
import os
from valis import registration, slide_io

registrar = registration.load_registrar("/mnt/project/VALIS/valis_files/micro/DAPIs/data/DAPIs_registrar.pickle")

dapi_folder="/mnt/project/VALIS/registered_DAPIs_micro/"
raw_image_dir="/mnt/project/raw_images_resized/"
dest_dir="/mnt/project/VALIS/registered_images_micro/"

for dapi in os.listdir(dapi_folder):
	base = dapi.split('_DAPI_')[0]

	cy5 = f"{raw_image_dir}{base}_Cy5_x20_z0.tif"
	cy7 = f"{raw_image_dir}{base}_Cy7_x20_z0.tif"
	dapi = f"{dapi_folder}{dapi}"
	
	cy5_dest = f"{dest_dir}{base}_Cy5_x20_z0.ome.tiff"
	cy7_dest = f"{dest_dir}{base}_Cy7_x20_z0.ome.tiff"

	if os.path.exists(cy5_dest):
        	print(f"Skipping {cy5_dest}, it already exists.")
        	continue

	print(f"Getting dapi slide: {dapi}")
	dapi_slide = registrar.get_slide(dapi)

	if os.path.exists(cy5):
		print(f"Processing image: {cy5}")
		dapi_slide.warp_and_save_slide(src_f=cy5, dst_f=cy5_dest, non_rigid=False, crop=True, channel_names=None, colormap=None, interp_method="bicubic", tile_wh=None, compression="lzw")

	if os.path.exists(cy7):
		print(f"Processing image: {cy7}")
		dapi_slide.warp_and_save_slide(src_f=cy7, dst_f=cy7_dest, non_rigid=False, crop=True, channel_names=None, colormap=None, interp_method="bicubic", tile_wh=None, compression="lzw")



