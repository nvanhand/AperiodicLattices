#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 08:03:07 2023

@author: nikkivanhandel
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 14:35:49 2023

@author: nikkivanhandel
"""

import matplotlib.pyplot  # as plt
import numpy as np
import cv2
from glob import glob
row_pts, col_pts = 9,10 # Size of checkerboard
path = '/Users/nikkivanhandel/Projects/Aperiodic Lattices/'

# Points prepare
objp = np.zeros((row_pts*col_pts, 3), np.float32)
objp[:, :2] = np.mgrid[0:col_pts, 0:row_pts].T.reshape(-1, 2)

objpoints = []  # 3d point in real world space
imgpoints = []  # 2d points in image plane.

img = cv2.imread(path+'chessboard_example.jpeg')
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
ret, thresh = cv2.threshold(gray, 127, 255, cv2.THRESH_BINARY)
msk = thresh.astype(np.uint8)

krn = cv2.getStructuringElement(cv2.MORPH_RECT, (30, 30))
dlt = cv2.dilate(msk, krn, iterations=5)
res = 255 - cv2.bitwise_and(dlt, msk)

h, w = img.shape[:2]

ret, corners = cv2.findChessboardCorners(res.astype(
    np.uint8), (row_pts, col_pts), flags=cv2.CALIB_CB_ADAPTIVE_THRESH + cv2.CALIB_CB_FAST_CHECK + cv2.CALIB_CB_NORMALIZE_IMAGE)

if ret == True:
    print('Success')
    objpoints.append(objp)
    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 30, 0.001)
    corners2 = cv2.cornerSubPix(gray, corners, (21, 21), (-1, -1), criteria)
    imgpoints.append(corners2)

    # Draw and copy corners
    #cv2.drawChessboardCorners(img, (row_pts, col_pts), corners, ret)
    # plt.imshow(img)
    #cv2.imwrite(f"{path}annotated_{row_pts}x{col_pts}.png", img)

    ret, mtx, dist, rvecs, tvecs = cv2.calibrateCamera(
        objpoints, imgpoints, gray.shape[::-1], None, None)

    newcameramtx, roi = cv2.getOptimalNewCameraMatrix(
        mtx, dist, (w, h), 1, (w, h))
    undistort_img = cv2.undistort(img, mtx, dist, None, newcameramtx)
    cv2.imwrite(f"{path}undistort_{row_pts}x{col_pts}.png", undistort_img)

'''
mapx, mapy = cv2.initUndistortRectifyMap(mtx, dist, None, newcameramtx, (w,h), 5)
dst = cv2.remap(img2, mapx, mapy, cv2.INTER_LINEAR)
# crop the image
x, y, w, h = roi
dst = dst[y:y+h, x:x+w]
#plt.imshow(dst)
#cv2.imwrite('calibresult.png', dst)

# Determine error 
mean_error = 0
for i in range(len(objpoints)):
    imgpoints2, _ = cv2.projectPoints(objpoints[i], rvecs[i], tvecs[i], mtx, dist)
    error = cv2.norm(imgpoints[i], imgpoints2, cv2.NORM_L2)/len(imgpoints2)
    mean_error += error
print( "total error: {}".format(mean_error/len(objpoints)) )

newcameramtx, roi = cv2.getOptimalNewCameraMatrix(mtx, dist, (w,h), 1, (w,h))
'''
