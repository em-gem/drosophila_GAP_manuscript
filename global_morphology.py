# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 18:32:58 2025

@author: Emily Gemmill
"""

import os
import numpy as np
import pandas as pd
from skimage.io import imread
from skimage.filters import threshold_li
from skimage.measure import label, regionprops
from skimage.morphology import closing, remove_small_objects
from skimage.io import imsave
from scipy.spatial import ConvexHull
from skimage.filters import try_all_threshold
import matplotlib.pyplot as plt

# Use the below lines to test thresholding methods if desired
# fig, ax = try_all_threshold(image, figsize=(10, 8), verbose=False)
# plt.show()


# Set path to processed tif images
root_dir = "PATH/TO/DIR"

# Iterate over genotype folders
for genotype in os.listdir(root_dir):
    genotype_path = os.path.join(root_dir, genotype)
    if not os.path.isdir(genotype_path):
        continue

    result_dict = {
        "Perimeter": {},
        "Circularity": {},
        "ConvexHullArea": {},
        "Solidity": {}
    }

    # Iterate over .tif files
    for BCellCluster in os.listdir(genotype_path):
        if not BCellCluster.endswith(".tif"):
            continue
        fpath = os.path.join(genotype_path, BCellCluster)

        # Load and binarize image
        image = imread(fpath)
        if image.ndim == 3:
            image = image[:, :, 0]  # convert RGB to grayscale
                
        # Threshold image           
        thresh = threshold_li(image)
        binary = closing(image > thresh)

        # Clean binary
        binary = remove_small_objects(binary, min_size=500)

        # Label and get largest region (assume 1 cluster per image)
        label_img = label(binary)
        props = regionprops(label_img)

        if not props:
            print(f"No object found in {BCellCluster}")
            continue

        # Use largest region
        region = max(props, key=lambda x: x.area)

        area = region.area
        perimeter = region.perimeter
        circularity = 4 * np.pi * area / (perimeter ** 2) if perimeter > 0 else np.nan

        coords = region.coords
        try:
            hull = ConvexHull(coords)
            convex_area = hull.volume  # for 2D, "volume" is area
        except:
            convex_area = np.nan

        solidity = area / convex_area if convex_area > 0 else np.nan

        # Save to dict
        result_dict["Perimeter"][BCellCluster] = perimeter
        result_dict["Circularity"][BCellCluster] = circularity
        result_dict["ConvexHullArea"][BCellCluster] = convex_area
        result_dict["Solidity"][BCellCluster] = solidity
        
        # Create a blank mask with the same shape
        binary_mask = np.zeros_like(label_img, dtype=np.uint8)
        binary_mask[label_img == region.label] = 255  # white mask for cluster

        # Save binary mask to same folder, with "_mask.tif" suffix
        mask_cluster = os.path.splitext(BCellCluster)[0] + "_mask.tif"
        mask_path = os.path.join(genotype_path, mask_cluster)

        imsave(mask_path, binary_mask.astype(np.uint8))

    # Save to CSV
    df = pd.DataFrame(result_dict)
    df = df.transpose()  # so rows = features
    output_path = os.path.join(root_dir, f"{genotype}_morphology_summary.csv")
    df.to_csv(output_path)

    print(f"Saved: {output_path}")
