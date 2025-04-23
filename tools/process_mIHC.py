import openslide
import cv2
from PIL import Image
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from skimage import transform
from openslide import OpenSlide
import tifffile
import os
import shutil
from skimage.color import rgb2hsv
from skimage.filters import threshold_otsu
from scipy.ndimage import binary_dilation, binary_erosion

import numpy as np
import tifffile
import matplotlib.pyplot as plt
import os
from skimage.measure import block_reduce
from skimage.exposure import rescale_intensity
import cv2


def adjust_brightness(input_path, target_brightness=150):

    savefig_path = input_path.replace(".tiff", "_adjust_brightness.tiff")    
    if os.path.exists(savefig_path):
        create_thumbnail_mIHC(savefig_path)
        print(f"Thumbnail exists - skipping: {os.path.basename(input_path)}")
        
        return savefig_path       
    
    with tifffile.TiffReader(input_path) as tif:
        image = tif.asarray()

        local_weight = cv2.blur((image > 0).astype(float), (51, 51))
        local_mean = cv2.blur(image.astype(float), (51, 51)) / (local_weight + 1e-10)
        image = (image * (target_brightness / local_mean).clip(0, 3)).clip(0, 255).astype(np.uint8)

    # save image
    with tifffile.TiffWriter(savefig_path) as tif:
        tif.write(image, compression=None, photometric='minisblack', tile=(256, 256))

    # save thumbnail of mIHC
    create_thumbnail_mIHC(savefig_path)      

    return savefig_path

def process_perc_clipped_image(input_path):

    savefig_path = input_path.replace(".tiff", "_cliped.tiff")
    if os.path.exists(savefig_path):
        create_thumbnail_mIHC(savefig_path)
        print(f"Processed image exists - skipping: {os.path.basename(input_path)}")
        return savefig_path

    with tifffile.TiffReader(input_path) as tif:
        image = tif.asarray()

        # 展平图像并过滤0值
        nonzero_values = image.flatten()[image.flatten() != 0]  # 一维数组
        # 展平图像为1D数组
        pixel_values = nonzero_values.flatten()

        q01, q99 = np.percentile(pixel_values, [1, 99]) # q01作为背景干扰

        cliped_image = (image - q01).clip(0) # remove bg

        # set cut-off using mean + 3 * std
        mean = np.mean(cliped_image)
        std = np.std(cliped_image)

        cliped_image = rescale_intensity(cliped_image, in_range=(0, q99), out_range=(0, 255))    

    # save image
    with tifffile.TiffWriter(savefig_path) as tif:
        tif.write(cliped_image, compression=None, photometric='minisblack', tile=(256, 256))

    # save thumbnail of mIHC
    create_thumbnail_mIHC(savefig_path)  

    return(savefig_path)

def downsample_manual(image, factor=2):
    """使用块平均降采样"""
    return block_reduce(image, block_size=(factor, factor), func=np.mean)


# 调整生成预览图的函数。
def create_thumbnail_mIHC(input_path, level=0, scale_factor=0.1):

    thumbnail_path = input_path.replace(".tiff", "_thumbnail.png")
    if os.path.exists(thumbnail_path):
        print(f"Thumbnail exists - skipping: {os.path.basename(input_path)}")
        return    

    with tifffile.TiffReader(input_path) as tif:
        image = tif.asarray()     
        image = downsample_manual(image, factor=10)  # 下采样10倍        
        
        plt.figure(figsize=(8,8))
        plt.imshow(image, cmap="gray")
        plt.axis("off")
        plt.title(f"{os.path.basename(thumbnail_path)}")        
        plt.savefig(thumbnail_path)     
        print(f"Saved thumbnail: {os.path.basename(thumbnail_path)}")

def scaled_image_mIHC(input_path, input_per_um_pixel=0.163343757390976, target_per_um_pixel=0.5, level=0):
    # save the scaled image
    savefig_path = input_path.replace(".tiff", "_scaled.tiff") # return path of scaled tiff
    if os.path.exists(savefig_path):
        print(f"Scaled image exists - skipping: {os.path.basename(input_path)}")
        return savefig_path

    scaled_factor = input_per_um_pixel / target_per_um_pixel
    with OpenSlide(input_path) as slide:
        width, height = slide.level_dimensions[level]
        img_scaled = np.array(
                slide.read_region((0, 0), level, (width, height))
                .convert('RGB')
                .resize((int(width*scaled_factor), int(height*scaled_factor)))
            )
        

        img_scaled = cv2.cvtColor(img_scaled, cv2.COLOR_BGR2GRAY)  # 注意OpenCV是BGR顺序


    with tifffile.TiffWriter(savefig_path) as tif:
        tif.write(img_scaled, compression=None, photometric='minisblack', tile=(256, 256))

    # save the thumbnail of scaled image
    create_thumbnail_mIHC(input_path=savefig_path)

    return savefig_path

# ---------------------------------------------------------
# Main processing
project_dir = "/aaa/zihanwu/yyyli2/project2_Spatial_CLIP"
data_dir = f"{project_dir}/data"
slides_dir = f"{project_dir}/slides"

# Create directories if they don't exist
os.makedirs(slides_dir, exist_ok=True)

samples_df_path = f"{project_dir}/data/Mouse_KPC_sample.csv"
samples_df = pd.read_csv(samples_df_path)
samples_df = samples_df[samples_df['finish'] == 1]
channels = ["DAPI", "520", "570", "700"]

# ---------------------------------------------------------------
print(f"\nStarting processing for {samples_df.shape[0]} samples")
for i in range(samples_df.shape[0]):
    dataset_id = samples_df["Dataset_ID"].iloc[i]    
    mIHC_slide_id = samples_df["mIHC_ID"].iloc[i]

    print(f"\nProcessing dataset: {dataset_id}")
    # Set up paths
    mIHC_slides_dir = f"{slides_dir}/{dataset_id}/mIHC"    
        
    os.makedirs(mIHC_slides_dir, exist_ok=True)

    # Copy files if needed
    source_mIHC_folder = f"{data_dir}/{dataset_id}/mIHC"
        
    if os.path.exists(source_mIHC_folder):
        for file in os.listdir(source_mIHC_folder):
            if file.lower().endswith(('.tiff', '.tif')):
                src_file = os.path.join(source_mIHC_folder, file)
                dst_file = os.path.join(mIHC_slides_dir, file)
                if not os.path.exists(dst_file):
                    shutil.copy2(src_file, dst_file)
        
    # Process mIHC channels
    for channel in channels:
        mIHC_channel_path = f"{mIHC_slides_dir}/{mIHC_slide_id}-{channel}.tiff"
        if not os.path.exists(mIHC_channel_path):
            print(f"Warning: mIHC channel {channel} not found at {mIHC_channel_path}")
            continue
        else:
            print(f"Processing mIHC channel {channel}: {mIHC_channel_path}")

        # downscaled to 0.5 um/pixel, save as scaled tiff and thumbnail tiff
        print(f"Scaled to 0.5 um/pixel...")    
        scaled_image_path = scaled_image_mIHC(input_path=mIHC_channel_path) # scaled tiff
        print(scaled_image_path)

        print(f"Cliped background of channel...")        
        clipped_image_path = process_perc_clipped_image(input_path=scaled_image_path) # remove background
        
        # add brightness
        print(f"Adjust brightness of channel...")        
        Adjusted_image_path = adjust_brightness(input_path=clipped_image_path)        


                    

print("\nProcessing complete.")        
