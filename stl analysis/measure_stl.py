#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 07:29:23 2023

@author: nikkivanhandel
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 21:29:55 2023

@author: nikkivanhandel
Inspiration to start: https://github.com/mikedh/trimesh/issues/1492

"""
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm

from shapely_geoms import *
from stl_to_shapely import *

aperiod = '/Users/nikkivanhandel/Projects/Aperiodic/'
stldir = aperiod + 'Imgs/STL Files/'

perturb = ['V0_N0_', 'V.2_N0_', 'V0_N.2_', ]
geom = ['H','S','T']


# Open STLs and turn into images
stl_data = pd.DataFrame()
w = 0.7
precision = 4

#
hist = False
vmin, vmax= 30, 180
m = 0

if hist:
    histcolor = ['red', 'blue', 'yellow']
    fig, axes = plt.subplots(3,1, figsize=(6,10))
else:
    fig, axes = plt.subplots(3,3, figsize=(12.5,10))

axeses = axes.flatten()
row0=True
for i in geom:
    
    colorpick = 0
    for j in perturb:
        name = j + i
        stl_f = f'{stldir}{name}.stl'

        poly = honeycomb_union(stl_f)
        simple = poly.simplify(.5) # Douglas-pecker
        
        gdf = geosummary(simple)
        gdf['angles'] = gdf.rings.map(get_angles)
        gdf['sides'] = gdf.geometry.map(shapely.get_num_coordinates)-1
        gdf['max_angle'] = gdf.angles.map(np.amax)
        
        if row0 and not hist: # Set Top Titles 
            axeses[m].set_title(j[:-1])
            
        if False: # Medial Axis info
            eroded = simple.buffer(-w/2+(10**-precision), join_style='mitre')
            avg = calculate_average_side_length(eroded.interiors, pr=precision-1)
            segments = list(avg)

        if hist:
            cnt, vals, obj = axeses[m].hist(gdf.max_angle, 
                           alpha=0.5,
                           color=histcolor[colorpick],
                           label=j[:-1], bins=np.arange(30, 180, 15))
            
        if  not hist and False:

            axeses[m] = gdf.plot(column="max_angle", ax=axeses[m], vmin=vmin, vmax=vmax, 
                                 cmap='coolwarm',
                                 edgecolor='white')
            
            axeses[m].plot(*simple.exterior.xy, 'k')
            axeses[m].set_yticks([])
            axeses[m].set_xticks([])
            for spine in  axeses[m].spines.values():
                spine.set_edgecolor('white')
        #plt.savefig(f'sides_{name}.png', bbox_inches='tight')
        if not hist and True:
            plt.clf()
            fig, ax = plt.subplots(figsize=(5,6))
            im = gdf.plot(column="sides", vmin=3, vmax=9, ax=ax,
                                 cmap='coolwarm',legend=False,
                                 edgecolor='white')
            y0 = 0
            handles =[]
            cw = cm.get_cmap('coolwarm')
            for x in range(3,9):
                handles.append(mpatches.Patch(color=cw(y0), label=str(x)))
                y0 += 1/len(range(3,9))
            plt.legend(handles=handles, loc='upper right')
            plt.plot(*simple.exterior.xy, 'k')
            plt.axis('off')
        
            plt.savefig(f'sides_{name}.png', bbox_inches='tight')
            
        colorpick += 1 
        m +=1
    row0=False   


if  hist:
    k = 5
    ytick = np.arange(0, 200, 25)
    axeses[0].legend()
    axeses[0].set_xticks([]), axeses[1].set_xticks([])
    axeses[2].set_xlabel('max_angle')
    for ax in range(3):
        axeses[ax].set_yticks(ytick)
        axeses[ax].set_ylabel(geom[ax])

    
    #plt.suptitle('Geometries by Y-Strut Length')
    plt.savefig('max_angle_distribution.png', dpi=500, bbox_inches='tight')

else:
    fig.subplots_adjust( wspace=0.025, hspace=0.025)

    fig = axeses[8].get_figure()
    cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    if False:
        cmap = plt.get_cmap('coolwarm', 7)
    sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    fig.colorbar(sm, cax=cax, ticks=np.arange(vmin, vmax, 15))

        
    for idx, a in enumerate([0,3,6]):
        axeses[a].set_ylabel(geom[idx])
    
    #plt.suptitle('Geometries by Y-Strut Length')
    plt.savefig('max_angle_grid.png', dpi=500, bbox_inches='tight')
    
    
    
    
    
    
    
