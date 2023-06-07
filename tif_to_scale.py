#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:11:22 2023

@author: nikkivanhandel
"""

from glob import glob
import cv2 
from skimage.measure import label, regionprops
from csv import DictWriter
import numpy as np
from realign import realign, bbox2


path = '/Users/nikkivanhandel/Projects/Aperiodic Lattices/'
# Where to read files from
read_from = path + '/tifs/'
write_to = path + '/out/'
# Where to write files to
dest = path + ''

def rename(key):
    # This function is to fix naming convention. Rewrite as needed.
    # Go from V0.2_N0_s -> V.2_N0_S 
    split = key.split('0.2')
    if len(split) > 1:
        name = split[0]+'.2'+split[1]
    else:
        name = key
    renamed = name[:-1] + name[-1].capitalize()

    return renamed

def process_scan(path, dest, im_dest=False):
    '''
    Takes path of a .TIF blue light scan. Assumes two images - a 10mm by 10mm calibration square and a much larger honeycomb image.  Optionally saves the larger honeycomb to an image destination folder. 
    Parameters
    ----------
    path : str
        path of images to read in.
    dest : str
        DESCRIPTION.
    im_dest : str, optional
        destination to create images . The default is False.

    Returns
    -------
    None.

    '''
    tifs = glob(read_from + '*.TIF')
    scale = {'geom': 'conversion'} # Initialize titles for csv file
    for tif in tifs:
        
        key = tif[len(read_from):-4] 
        name = rename(key)
        #print(name)   
        img = cv2.imread(tif, cv2.IMREAD_GRAYSCALE) 
        thresh = img < 255 # Assuming bg is white 
        lab = label(thresh) # Connectivity map
        lat_props, box_props = regionprops(lab)
        box_x1, _, box_x2, _ = box_props['bbox']
        #box_props = regionprops_table(lab, properties=('bbox', 'area'))
        if im_dest != False: # if you want to save a new image
            rot_mat = realign(lab==1, lat_props)
            lattice = bbox2(rot_mat)
            cv2.imwrite(write_to+name+'.png', lattice*255)

        y = box_x2 - box_x1
        scale[name] = y/5 # px per mm
    return scale
    
scale = process_scan(path, dest, im_dest)

csv_file = dest + "tif_scale.csv"
with open(csv_file, 'w') as f:
    for key in scale.keys():
        f.write("%s,%s\n"%(key,scale[key]))