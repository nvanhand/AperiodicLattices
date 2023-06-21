#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 21:29:55 2023

@author: nikkivanhandel
Inspiration to start: https://github.com/mikedh/trimesh/issues/1492

"""
import trimesh
import geopandas as gpd
import matplotlib.pyplot as plt
import shapely 
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning) 

def compare(out, inn):
    outer = gpd.GeoSeries(out)
    inner = gpd.GeoSeries(inn)
    ax = gpd.GeoSeries.plot(outer, figsize=(12,12), edgecolor='b', linewidth=1)
    gpd.GeoSeries.plot(inner, ax=ax, edgecolor='r', linewidth=1)
    ax.axis('off')
    return ax

def tround(tuples_list, pr=3):
    return [tuple(round(value, pr) for value in tuple_item) for tuple_item in tuples_list]


def calculate_average_side_length(linear_rings, pr=3):
    total_length = 0
    total_sides = 0

    unique_segments = set()

    for ring in linear_rings:
        coords = tround(list(ring.coords), pr)
        for i in range(len(coords)-1):
            p1 = coords[i]
            p2 = coords[i + 1]

            sorted_segment = tuple(sorted([p1, p2]))

            if sorted_segment in unique_segments:
                continue
            length = ((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2) ** 0.5

            total_length += length

            unique_segments.add(sorted_segment)

            # Increment the count of total sides
            total_sides += 1
      
    average_length = total_length / total_sides if total_sides > 0 else 0

    return unique_segments


def honeycomb_union(stl_f):
    # Loads STL and returns Shapely multipolygon
    mesh = trimesh.load_mesh(stl_f, enable_post_processing=True, solid=True) # Import Objects
    
    # get a single cross section of the mesh
    section = mesh.section(plane_origin=[0,0,0], 
                         plane_normal=[0,0,1])
    slice_2D, _ = section.to_planar()
    
    # Multipolygon
    poly_union = shapely.geometry.MultiPolygon([poly for poly in slice_2D.polygons_full])
    return poly_union.geoms[0]

from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from shapely.geometry import LineString
from shapely.ops import unary_union

def plot_segments(line_segments):
    # Calculate the lengths of line segments
    lengths = [LineString(segment).length for segment in line_segments]
    
    # Normalize the lengths to a range from 0 to 9
    length_min = min(lengths)
    length_max = max(lengths)
    normalized_lengths = []
    denom = length_max - length_min
    if denom>0:
        normalized_lengths = [(length - length_min)/denom*9 for length in lengths]
    else:
        normalized_lengths = [length for length in lengths]
    
    # Create a colormap
    cmap = plt.cm.get_cmap('coolwarm')
    
    # Create a color map based on the normalized lengths
    norm = Normalize(vmin=0, vmax=8.5)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    colors = [sm.to_rgba(length) for length in normalized_lengths]
    
    # Create a figure and axis
    fig, ax = plt.subplots(facecolor='white')

    
    # Plot each line segment with its corresponding color
    for segment, color in zip(line_segments, colors):
        ax.plot(*zip(*segment), color=color)
    
    # Set the aspect ratio to equal
    ax.set_aspect('equal')
    ax.set_facecolor('black')
    # Show the colorbar
    cbar = plt.colorbar(sm, ticks=[0,2,4,6,8])
    ax.axis('off')
    # Show the plot
    return fig
