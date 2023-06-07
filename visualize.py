#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 13:47:34 2023

@author: nikkivanhandel
"""

import matplotlib.pyplot as plt
from math  import sin, cos
from preprocess import *
from matplotlib.colors import ListedColormap
from seaborn import heatmap

def color_struts(colored, title=None):
    fig, ax = plt.subplots()
    colored = cv2.dilate(colored.astype(np.float32), np.ones((15,15)))
    mask = colored != 0
    im1= ax.imshow(colored, cmap= 'coolwarm')
    fig.colorbar(im1)
    ax.imshow(np.ma.masked_array(colored, mask), cmap='binary_r', vmin=1, vmax=np.max(colored))
    ax.axis('off')
    if title is not None:
        plt.suptitle(title)
    return fig

def color_by_thickness(im, segments, minval=50):
    dt = dt_edt(im)
    label = measure.label(segments)
    props = measure.regionprops(label, cache=False)
    colored = np.zeros_like(im).astype(np.float32)
    for p in props:
        if p.area > minval:
            coords = p['coords']
            tr = tuple(zip(*coords))
            colored[tr] = np.mean(dt[tr])
    return colored

def pixel_represent(arr,annot=True, mask=None):
    fig, ax = plt.subplots()
    heat = heatmap(arr, cmap='viridis', ax=ax,
                   xticklabels=[], yticklabels=[], 
                   cbar=False, linewidths=1,
                   annot=annot)
    if mask is not None:
        for text, show_annot in zip(ax.texts, mask.ravel()):
            text.set_visible(show_annot)
            

def comp_sk(skim, node_im, dil, bounds=[0, 600, 0, 600]):
    fig, axes = plt.subplots(1,3)
    x1, x2, y1, y2 = bounds
    vers = [skim, node_im, logical_minus(skim, dil)]
    kern = np.ones((15, 15))

    for v in range(3):
        vers[v] = cv2.dilate(((vers[v])[x1:x2, y1:y2]).astype(np.uint8), kern)
    vers[1] = vers[0] - vers[1]
    til = ['skeleton', 'nodeless', 'segments']
    for i, ax in enumerate(axes.reshape(-1)):
        ax.axis('off') 
        #ax.set_title(til[i])
        ax.imshow(vers[i])
    
    plt.savefig('test.png', dpi=500, bbox_inches='tight')
        

def compare_thicknesses(cell_im, skim, node_im):
    thick1 = thickness(cell_im, skim)
    dist_edt = dt_edt(cell_im)
    dil_nodes = expand_nodes(node_im, k=round(thick1))
    thick2 = np.mean(dist_edt[segment_im])
    thick3 = np.mean(dist_edt[logical_minus(skim, dil_nodes)])
    return [thick1, thick2, thick3]

    df  = pd.DataFrame(thicknesses, columns=['specimen', 'skeleton', 'sk-nodeless', 'segments'] )
    df = 2*df.set_index('specimen')
    return df

def discrete_cbar(image, colormaps='viridis'):
    ticks = np.arange(np.amin(image), np.amax(image+1))

    # Create a colormap by combining colors from existing colormaps
    colors = []
    cmap = plt.get_cmap(colormaps)
    colors.extend(cmap(np.linspace(0, 1, len(ticks))))
    newcmap = ListedColormap(colors)
    return newcmap, ticks

def node_connectivity_map(cell_im, c):
    fig, ax = plt.subplots()
    im1 = ax.imshow(cell_im, cmap='binary_r')
    ax.axis('off')
    ax.patch.set_edgecolor('black')  
    ax.patch.set_linewidth('1')  

    cmask = np.ma.masked_where(c==0, c)
    custom, ticks = discrete_cbar(cmask)
    im2 = ax.imshow(cmask, cmap=custom, alpha=0.6)
    
    cbar = ax.figure.colorbar(im2, ax=ax, ticks=ticks)
    cbar.ax.set_yticklabels(ticks)
    return fig


def see_nodes(cells, dnodes, nodes, segments):
    white = logical_minus(cells, dnodes)
    red = logical_minus(cells, segments)
    yellow = logical_minus(cells, nodes)
    stack = np.dstack((white, red, yellow))
    return stack
    #cv2.imwrite(push+'nodes_'+name+'.png', 255*stack)

def draw_min_maj(image):
    # Visualize 
    label_img = measure.label(image)
    regions = measure.regionprops(label_img)
    # Remove BG if necessary
    if label_img[0,0] != 0:
        label_img[label_img==label_img[0,0]] = 0
    fig, ax = plt.subplots()
    ax.imshow(image, cmap=plt.cm.gray)
    
    for p in regions:
        if p.area > 1000:
            y0, x0 = p.centroid
            pmax, pmin =0.5* p.axis_major_length, 0.5*p.axis_minor_length
            orient = p.orientation
            x1 = x0 + cos(orient) * pmin
            y1 = y0 - sin(orient) * pmin
            x2 = x0 - sin(orient) * pmax
            y2 = y0 - cos(orient) * pmax
        
            ax.plot((x0, x1), (y0, y1), '-r', linewidth=1)
            ax.plot((x0, x2), (y0, y2), '-b', linewidth=1)
            ax.plot(x0, y0, '.g', markersize=1)
    plt.axis('off')
    return fig

def plot_hist(skel, seg_reg, space_reg, path, filename):
    segmap, spmap = np.zeros_like(skel*1), np.zeros_like(skel*1) 
    seglen = []
    for seg in seg_reg:
        coord, dia = seg.coords, seg.major_axis_length
        segmap[tuple(zip(*coord))] = dia
        seglen.append(dia)
    
    spsz=[]
    for s in space_reg:
        coord, dia = s.coords, s.area
        spmap[tuple(zip(*coord))] = dia
        spsz.append(dia)
    spsz = sorted(spsz)
    
    plt.figure()
    plt.imshow(spmap, cmap='hot', vmax=spsz[-2])
    plt.axis("off")
    
    plt.savefig(path + "cellA-"+ filename + ".png", dpi=1000, bbox_inches="tight")
    
    fig, ax = plt.subplots(1,2)
    spplt = ax[0].imshow(spmap, cmap='hot', vmax=spsz[-2])
    plt.colorbar(spplt,ax =  ax[0])
    
    segplt = ax[1].imshow(segmap, cmap='hot', vmin=1)
    plt.colorbar(segplt,ax = ax[1])
    
    ax[0].axis('off')
    ax[1].axis('off')
    
    ax[0].title.set_text('Cell Area')
    ax[1].title.set_text('Max Axis Length')
    
    plt.savefig(path + "map-"+ filename + ".png", dpi=1000, bbox_inches="tight")
    
    fig, ax = plt.subplots(1,2)
    
    ax[0].hist(seglen, range=(1, max(seglen)))
    ax[1].hist(spsz, range=(spsz[0], spsz[-2]))
    ax[0].title.set_text('Cell Area')
    ax[1].title.set_text('Max Axis Length')   
    
def quickplot(im, cax=None):
    if cax is None:
        fig, ax = plt.subplots()
    ax.imshow(im)
    ax.axis('off')
    