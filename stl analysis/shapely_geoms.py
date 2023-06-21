#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 07:22:59 2023

@author: nikkivanhandel
"""
import numpy as np
from math import sqrt, degrees, acos
import pandas as pd
import geopandas as gpd
import shapely 

def vertex_angles(points):
    angles = []
    points = points[:-1]
    num_points = len(points) 
    
    for i in range(num_points):
        p1 = points[i]
        p2 = points[(i + 1) % num_points]
        p3 = points[(i + 2) % num_points]
        
        v1 = p1 - p2
        v2 = p3 - p2
        
        norm_v1 = np.linalg.norm(v1)
        norm_v2 = np.linalg.norm(v2)
        
        if norm_v1 != 0 and norm_v2 != 0:
            dot_product = np.dot(v1, v2)
            norm_product = norm_v1 * norm_v2
            
            angle = np.arccos(dot_product / norm_product)
            angles.append(degrees(angle))
        else:
            print(p1, p2, p3)
    return angles

def min_maj(cell):
    # from https://stackoverflow.com/questions/13536209/efficient-way-to-measure-region-properties-using-shapely
    mrr = cell.minimum_rotated_rectangle.exterior.coords
    lens = []
    for n in range(4):
        x1, y1 = mrr[n]
        x2, y2 = mrr[n+1] 
        d = ((x2-x1)**2+(y2-y1)**2)**0.5
        lens.append(d)
    
    # get major/minor axis measurements
    minor_axis = min(lens)
    major_axis = max(lens)
    return minor_axis, major_axis

def get_angles(poly): 
    coord = np.array(poly.coords)
    return vertex_angles(coord)

def geosummary(shapely_lattice):
    glr = gpd.GeoSeries(shapely_lattice.interiors)
    gp = glr.map(shapely.geometry.Polygon)
    gdf = gpd.GeoDataFrame({'geometry': gp, 'rings': glr})
    return gdf

def measure_polygon(cell):
    # Input: individual polygon 
    # Returns: width, height, aspect, n_sides, 
    #          angles, perimeter, area
    geom = shapely.geometry.Polygon(cell)
    x1, y1, x2, y2 = cell.bounds
    w, h= x2-x1, y2-y1
    pts = cell.coords

    return w, h, w/h, len(pts)-1,  cell.length, geom.area,vertex_angles(pts)

def delaunay_stats(poly):
    # Input: Honeycomb shapely
    # Returns delta, avg edge, std of edges
    #https://stackoverflow.com/questions/69443921/how-can-i-count-the-length-of-the-edge-associated-with-each-point
    
    from scipy.spatial import Delaunay
    from scipy.spatial.distance import pdist
    
    points = np.array(centroids(poly.interiors))
    tri = Delaunay(points)
    edges = np.array([pdist(x) for x in points[tri.simplices]])
    
    return [np.amin(edges), np.mean(edges), np.std(edges)]

def centroids(poly_ring):
    return [x.centroid.coords[0] for x in poly_ring]

def exterior_stats(poly, code='h'):
    # Structural Area, Width, Height, 
    # NegativeArea, NumberofCells, Zhu's d0
    mult = {'h': 2/(3**0.5), 's': 1, 't': 1/(3**1.5)}

    bb = poly.exterior.bounds
    w, h = bb[2] - bb[0], bb[3]-bb[1]
    N = len(poly.interiors)
    d0 = (w*h/N*mult[code])**0.5
    
    nA = (bb[2]-bb[0])*(bb[3]-bb[1]) - poly.area
    return [poly.area, w,h, nA, N, d0]


def analyze_lattice(poly, code='h'):
    # input whole shapeply lattice
    ext_stats = dict(zip(['structural_area', 'lattice_width', 'lattice_height', 'negative_area', 'total_cells', 'd0'], exterior_stats(poly)))
    del_stats= dict(zip(['delta', 'mean_del', 'std_del'], delaunay_stats(poly)))
    df = []
    cell_stats = {}
    for cell in poly.interiors:
        df.append(measure_polygon(cell))
    df = pd.DataFrame(df, columns=['x', 'y', 'aspect','n_sides', 'perimeter', 'area', 'angles'])
    
    angles_summary = df['angles'].apply(lambda x: pd.Series([np.min(x), np.max(x), np.mean(x), np.std(x)]))
    angles_summary.columns = ['min_angle', 'max_angle', 'mean_angle', 'std_angle']
    all_data = pd.concat([df.drop('angles',axis=1),
                    angles_summary], axis=1)
    
    summary = all_data.describe()
    for col in all_data.columns:
        for stat in ['max', 'min', 'mean', 'std']:
            cell_stats[f"{stat}_{col}"] = summary.loc[stat, col]
    cell_stats.update(ext_stats)
    cell_stats.update(del_stats)
    return cell_stats

## Visualization

def plot_delaunay(simple, title='Delaunay.png'):
    # Plot the delaunay triangulation 
    import matplotlib.pyplot as plt
    from scipy.spatial import Delaunay, delaunay_plot_2d 
    
    
    c = centroids(simple.interiors)
    fig, ax = plt.subplots(1)
    ax.plot(*simple.exterior.xy, 'k')
    for poly in simple.interiors:
        plt.plot(*poly.xy, 'k')
    
    dtx = Delaunay(np.array(c))
    delaunay_plot_2d(dtx, ax)
    ax.axis('equal')
    plt.axis('off')
    
    plt.savefig(title, bbox_inches='tight', dpi=200)
    
    
