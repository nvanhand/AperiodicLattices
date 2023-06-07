#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 12:56:45 2023

@author: nikkivanhandel
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 13:58:19 2022

@author: nikkivanhandel

This script is mostly just functions that I use in the analysis script. 
Functions are used on a per image basis

"""

import cv2
import numpy as np
import pandas as pd
from skimage import measure
from scipy.ndimage import convolve
from scipy.ndimage import binary_dilation
from scipy.ndimage import distance_transform_edt as dt_edt
from skimage.morphology import skeletonize

def logical_minus(keep, *arrays):
    n = len(arrays)
    diff = keep.astype(int) * n
    for arr in arrays:
        diff -= arr.astype(int)
    
    return diff==n

def process(p): 
    '''
    Performs preprocessing for lattice images. Intended for Keyence/Blue light scans. 

    Parameters
    ----------
    p : str
        path to file to read and process.

    Returns
    -------
    image : bool array
        Binarized lattice image.

    '''
    gray = cv2.imread(p, cv2.IMREAD_GRAYSCALE)
    blur = cv2.GaussianBlur(gray, (13,13),0)
    _, image = cv2.threshold(blur, np.mean(blur), 1, cv2.THRESH_BINARY)
    if np.sum(image) > image.size/2:
        image = np.logical_not(image)
    return image

def segmentize(image):
    '''
    Accepts a 2D lattice and returns the nodes and segments of the lattice.  
    Parameters
    ----------
    image : bool type image array
        Preprocessed image of a lattice (from process fxn).

    Returns
    -------
    image : bool type array
        Standardizes the array s.t. the positive space is 1 and negative space is 0.
    nodes : bool type array
        image with all skeleton node points = 1.
    segments : bool type array
        skeleton minus the nodes - completely segmented struts of thickness 1.

    '''
    sk = skeletonize(image).astype(np.uint8)
    kernel = np.ones((3,3))
    kernel[1,1] = 9 # this means that the middle pixel is most weighted
    nodes = convolve(sk, kernel, mode='constant') >= 12 
    nodes = cv2.dilate(nodes.astype(np.uint8), np.ones((3,3)))
    segments = ((sk+1) - nodes)>1
    return nodes.astype(np.bool), segments #, sk.astype(np.bool)


def expand_nodes(bit_nodes, k=15):
    '''
    Isolates the full size struts from the binary lattice using the information from segmentize.  

    Parameters
    ----------
    cell_im : bool type image array
        Binary image of the lattice where the lattice is 1 and background is 0.
    bit_nodes : bool type image array
        Binary image of lattice nodes (from segmentize). Nodes are 1px in size
    k : int
        dilation size (based on thickness)

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    # Removes the nodes from the binary image to isolate struts
    kern = np.ones((k,k))
    nodes = cv2.dilate(bit_nodes.astype(np.uint8), kern)
    nodes = cv2.morphologyEx(nodes, cv2.MORPH_CLOSE, kern)
    ell_kern =  cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(k,k))
    dil_nodes = cv2.dilate(nodes.astype(np.uint8), ell_kern)
    # Find Corners using Harris 
    #dst = cv2.cornerHarris(cell_im.astype(np.float32),k,3,0.01)
    #dst = cv2.dilate(dst, np.ones((k,k)))
    #corners = dst > 0.001
    return dil_nodes

def connectivity(segments, nodes, dil_nodes):
    '''
    Determine how many struts meet at a given node

    Parameters
    ----------
    segments : TYPE
        bool array where all strut segments are .
    nodes : bool array
        Bool array where all nodes are marked as 1.
    dil_nodes : TYPE
        DESCRIPTION.

    Returns
    -------
    node_connectivity : TYPE
        DESCRIPTION.

    '''
    # How many branches per region? 
    #connect_map = dil_nodes.copy()   # Unsuppress for 
    node_labels = measure.label(dil_nodes)
    node_props = measure.regionprops(node_labels, cache=True)
    c = []
    for r in node_props:
        (x1, y1, x2, y2) = r.bbox
        sgno, _ = cv2.connectedComponents(segments[x1:x2, y1:y2].astype(np.uint8))
        nono, _ = cv2.connectedComponents(nodes[x1:x2, y1:y2].astype(np.uint8))
        total = sgno - nono + 1
        c.append(total)
        #connect_map[node_labels==r.label] = total
    node_connectivity = np.mean(c)
    return node_connectivity
        
    
def thickness(im, nodes, sk):
    dist_edt = dt_edt(im)
    strut_thick = 2*np.mean(dist_edt[sk])
    node_thick = np.mean(dist_edt[nodes])
    return strut_thick, node_thick


def color_by(labelled, reg, feat='area',  amin=0):
    '''
    Replaces labels with desired feature values

    Parameters
    ----------
    labelled : int32 array
        result of measure.label.
    regions : pandas DataFrame
        DataFrame from regionprops_table.
    feat : TYPE, optional
        DESCRIPTION. The default is 'area'.

    Returns
    -----
    colored : TYPE
        DESCRIPTION.

    '''
    colored = np.zeros_like(labelled)
    feat_list = []
    for r in reg:
        coords, f = r['coords'], r[feat]
        colored[tuple(zip(*coords))] = f
        feat_list.append(f)
    return colored, feat_list

def shape_analyze(im, bg=0, minval=0, crit='area', table=False, xprops=[]):
    '''
    processing for a labelled connected components image

    Parameters
    ----------
    label_im : int32 array
        result of measure.label.
    bg : bool, either 0 or 1
        What value to treat as the background. The default is 0.
    amin : int, optional
        Optionally filter the minimum size of components. The default is 0.
    xprops: list of str, optional
        Additional properties besides those assumed (see props)
    
    Returns
    -------
    label_im: int32 array
        label image of connected components.
    tab: pandas DataFrame
        region properties (adjust based on props below)

    '''
    props = ['label', 'area', 'orientation', 
             'minor_axis_length', 'major_axis_length'] + xprops

    label_im = measure.label(im, background=bg)
    if bg != 0: #if characterizing negative space 
        label_im[np.where(label_im>=1)] -= 1 # Remove bg
    regions = measure.regionprops(label_im, cache=False)
    
    if minval > 0:
        for p in regions:
            if p[crit] < minval:
                label_im[label_im==p.label] = 0
                label_im[label_im>p.label] -= 1
      
    if table:
        propout = measure.regionprops_table(label_im, properties=props)
                                          #properties = props)
    else:                                   
        propout = measure.regionprops(label_im, cache=False)
    #tab = pd.DataFrame(tab)
    return label_im, propout
    



'''
# Attempted to use watershed segmentation to isolate. 
# Didn't work, but know knows if i'll ever want that again. 
def segment_centroids(binary_image, segments):
    struts = measure.label(segments)
    tab = measure.regionprops_table(struts, 
                                    properties = ['centroid'])
    centroids = []
    markers = np.zeros_like(binary_image).astype(np.uint16)
    n=1
    for n in range(len(tab['centroid-0'])):
        c1, c2 = int(tab['centroid-0'][n]), int(tab['centroid-1'][n])
        markers[c1, c2] = n
        n+=1

    return markers 
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
from scipy.ndimage import distance.... transform somethign or other. 
dist = dt_edt(cell_im)
loc_max = peak_local_max(dist, min_distance=51, labels=cell_im)
markers = np.zeros_like(cell_im, dtype=np.uint8)
markers[loc_max[:, 0], loc_max[:, 1]] = 1

'''