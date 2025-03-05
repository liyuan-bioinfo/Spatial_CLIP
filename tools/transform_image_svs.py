# transform image to svs

import pandas as pd
import numpy as np
import torch

import scanpy as sc
import pandas as pd
import numpy as np

import os
import matplotlib.pyplot as plt

import openslide
from PIL import Image

import tifffile

import cv2
image = cv2.imread("UKF_304_RT.tif")
output_path = "GBM_304.tiff"

# save tiff
with tifffile.TiffWriter(output_path) as tif:
    tif.save(image, photometric='rgb', compression='lzw', tile=(256, 256))
    
    for level in range(1, 4):  # 保存 3 层金字塔
        downsampled = image[::2**level, ::2**level]  # 缩小图像
        tif.save(downsampled, photometric='rgb', compression='lzw', tile=(256, 256))

# openslide
# test whethor could be open with thumbnail image
slide = openslide.OpenSlide("/aaa/zihanwu/yyyli2/10XVisium/UKF304_T_ST/HE/GBM_304.tiff")
target_width = 2000

# 获取最高分辨率层级的尺寸（level=0）
base_width, base_height = slide.dimensions

# 计算缩放比例
scale_factor = target_width / base_width
scaled_size = (int(base_width * scale_factor), int(base_height * scale_factor))

thumbnail = slide.get_thumbnail(scaled_size)

# 显示或保存
thumbnail.save("preview_thumbnail_gbm_304.jpg")  # 保存为文件
# thumbnail.show()
