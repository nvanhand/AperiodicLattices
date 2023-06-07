#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 10:29:04 2023

@author: nikkivanhandel
"""

import numpy as np
import cv2




def bbox2(img):
    rows = np.any(img, axis=1)
    cols = np.any(img, axis=0)
    ymin, ymax = np.where(rows)[0][[0, -1]]
    xmin, xmax = np.where(cols)[0][[0, -1]]
    return img[ymin:ymax+1, xmin:xmax+1]


def realign(img, region_props=False):
    '''
    Takes in a given grayscale image and the correpsonding labelled region of interest. Realigns image about the centroid of the region .

    Parameters
    ----------
    img : grayscale array
    region_props : region_props item

    Returns
    -------
    rot_img : rotated img.
    '''
    if region_props == False:
        img
    if type(region_props) == dict:
        points = region_props['coords'] 
        xc, yc = region_props['centroid']
    elif type(region_props) == list:
        points = region_props[0].coords
        xc, yc = region_props[0].centroid
    else:
        points = region_props.coords
        xc, yc = region_props.centroid
    _, (w,h), rot = cv2.minAreaRect(points)
    if rot < 90-rot:
        rot = -rot
    else:
        rot = 90-rot
    M = cv2.getRotationMatrix2D((xc, yc), rot ,1) 
    rot_img = cv2.warpAffine(img.astype('float32'),
                             M, img.shape[:2], 
                             borderMode = cv2.BORDER_CONSTANT,
                             )
    return rot_img


def crop_lattice(im, save=True):
    scale, delta, ddepth = 1, 0, cv2.CV_8UC1
    gray = cv2.cvtColor(im, cv2.COLOR_RGB2GRAY)
    opp = np.logical_not(gray>200)
    coords = np.argwhere(opp!=0)
    
    x_min, y_min = coords.min(axis=0)
    x_max, y_max = coords.max(axis=0)
    cropped = cv2.bitwise_not(gray)[x_min:x_max+1, y_min:y_max+1]   
    if save:
        cv2.imwrite(f, cropped)
