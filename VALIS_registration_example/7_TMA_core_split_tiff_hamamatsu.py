import os
import tifffile
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np
from skimage import io, img_as_uint
from argparse import ArgumentParser

Image.MAX_IMAGE_PIXELS = None


def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument("core_name", type=str)
    parser.add_argument("x", type=int)
    parser.add_argument("y", type=int)
    parser.add_argument("side", type=int)
    parser.add_argument("full_res_image", type=str)
    parser.add_argument("subsampled_x", type=int)
    parser.add_argument("subsampled_y", type=int)
    parser.add_argument("output_path", type=str)
    parser.add_argument("image_name", type=str)
    parser.add_argument("--mask", default=False, action="store_true")
    args = vars(parser.parse_args())
    return args

def fromarray_large(array, mode=None):
    """ This replaces PIL's Image.fromarray() method for images with more than 2^28 pixels.
    Image is split into N pieces horizontally and the pieces are then pasted ( Image.paste() )
    into PIL image to circumvent using fromarray() fo large images.

    Args:
        array: input image as numpy ndarray
        mode: PIL image mode for the output image, if none mode is deduced automatically

    Returns:
        Pil image that has the same size and content as the input array
    """

    # Use fromarray if image is small enough
    if array.shape[0] * array.shape[1] < 2 ** 28:
        print("fromarray_large(): Using fromarray")
        return Image.fromarray(array, mode=mode)

    else:
        print("fromarray_large(): Using fromarray large")
        if mode is None:
            # get mode by creating one pixel image
            mode = Image.fromarray(array[0:1, 0:1]).mode
            print("fromarray_large(): Determined output image mode: {}".format(mode))

        # divide image into N equal sized pieces.
        # Here the N is selected to be as large as possible so that
        # size of every piece is below 2^28.
        N = int(np.ceil((array.shape[0] * array.shape[1]) / (2 ** 28)))
        print("fromarray_large(): Splitting image into N={} pieces".format(N))

        # divide array's x-axis into N slices
        piece_inds = np.linspace(0, array.shape[1], N + 1).astype(int)

        # create empty image having the same size as input array
        # PIL uses coordinates in order x, y whereas numpy uses y, x
        pil_im = Image.new(mode, array.shape[:2][::-1])  # 8-bits / channel (3x8)

        # iterate over all pieces and paste them into the output image
        y_min = 0
        y_max = array.shape[0]
        for i in range(N):
            x_min = piece_inds[i]
            x_max = piece_inds[i + 1]
            # left, upper, right, lower coordinates of the pasted image in the output coordinates
            box = (x_min, y_min, x_max, y_max)
            print(
                "fromarray_large(): Pasting piece {} with box coordinates (l,t,r,b) {} "
                "into image with size {}".format(i, box, pil_im.size)
            )
            piece = Image.fromarray(array[:, x_min:x_max])
            pil_im.paste(piece, box)

        return pil_im

def read_image_region(filepath, x1, y1, x2, y2):
    file_ext = os.path.splitext(filepath)[1].lower()
    print(file_ext)
    if file_ext in [".tiff", ".tif", ".ome.tiff", ".ome.tif"]:
        print("reading tiff with tifffile")
        with tifffile.TiffFile(filepath) as tif:
            if hasattr(tif, 'series') and len(tif.series) > 0:
                data = tif.series[0].asarray()
            else:
                data = tif.asarray()
            im = data[y1:y2, x1:x2]
            im = im.astype(np.uint8)
    else:
        with Image.open(args['full_res_image']) as pil_image:
            im = np.array(pil_image.crop([x1, y1, x2, y2]), dtype='uint8')
    return im

def get_image_dimensions(filepath):
    file_ext = os.path.splitext(filepath)[1].lower()
    
    if file_ext in ['.tiff', '.tif', '.ome.tiff', '.ome.tif']:
        with tifffile.TiffFile(filepath) as tif:
            if hasattr(tif, 'series') and len(tif.series) > 0:
                shape = tif.series[0].shape
            else:
                shape = tif.asarray().shape
            
            if len(shape) == 2:
                h, w = shape
            else:
                raise ValueError(f"Unexpected image shape: {shape}")
            
            return w, h
    else:
        with Image.open(filepath) as pim:
            return pim.size



args = parse_arguments()
os.makedirs(os.path.join(args['output_path'], args['core_name']), exist_ok=True)
table_start = False

# Get ratio of downsampled and upsampled images
# downsampled image dimensions
h1 = args["subsampled_y"]
w1 = args["subsampled_x"]

# full_res_image_dimensions
w2, h2 = get_image_dimensions(args['full_res_image'])
rh = h2/h1
rw = w2/w1

# Process one line from TMA grid file (one line is passed as an argument)
core = args['core_name']
print(core)
x = args['x']
y = args['y']
x1 = x*rw
y1 = y*rh

# Use set side length (2500 - 4000 ok for analysis)
side = 4000

# Use core length from qptma
#side = args['side']*rh

x2 = x1+side
y2 = y1+side

x1 = np.round(x1).astype(np.int32)
y1 = np.round(y1).astype(np.int32)
x2 = np.round(x2).astype(np.int32)
y2 = np.round(y2).astype(np.int32)

print(f"Cropping region (x1, y1, x2, y2): {[x1, y1, x2, y2]}")

# Read and crop the image region
im = read_image_region(args['full_res_image'], x1, y1, x2, y2)

im_name = args['image_name'] + '_{}_{}.png'.format(x1, y1)
path_out = os.path.join(args['output_path'], args['core_name'], im_name)
print(path_out)

print(f"Core piece maximum: {np.amax(im)}")
print(f"Core piece minimum: {np.amin(im)}")
print(f"Image dtype: {im.dtype}")
print(f"Image shape: {im.shape}")

if args["mask"]:
    im = fromarray_large(im, 'L')
else:
    im = fromarray_large(im)

im.save(path_out)



