#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 16:15:33 2023

@author: nikkivanhandel
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 19:07:49 2022

@author: nikkivanhandel
https://docs.opencv.org/3.4/d2/dbd/tutorial_distance_transform.html

"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 16:13:17 2022

@author: nikkivanhandel
"""
from preprocess import *
from visualize import *

def y_comp(r):
    y1, x1, y2, x2 = r.bbox
    return y2-y1

def aniso(r):
    return r['major_axis_length']/r['minor_axis_length']

def orient(r):
    return np.degrees(r['orientation'])


path = '/Users/nikkivanhandel/Projects/Aperiodic/'
impath = path + '/Imgs/'

# Read scale values for tif images 
tif_scale = pd.read_csv(f'{path}/code/image/tif_scale.csv')
tif_scale.set_index('geom', inplace=True)

# Lazy way to iterate through files
perturb = ['V0_N0_', 'V.2_N0_', 'V0_N.2_', ]
geom = ['H','S','T']

# Compatability for both types of images 
def choose(name, which):
    tifc = float(tif_scale.loc[name, 'conversion'])
    ch = {'tif': [tifc, 'BlueLight/',  '.TIF'],
          'png': [21.3, 'Keyence/', '.png'],
          'stlim': [1, 'STL Images/', '.png']}
    return ch[which]

ims = ['bin_', 'skel_', 'node_', 'seg_']


# (12.5, 10) for 9x9 with colorbar

fig, axes = plt.subplots(3,1, figsize=(6,10))

m = 0
axeses = axes.flatten()

histcolor = ['red', 'blue', 'yellow']
cnt = []
cntmin, cntmax = 100, 0
valmin, valmax = 100, 0

row0 = True
traits = []

for i in geom:
    colorpick = 0
    for j in perturb:
        name = j+i
        print(name)
        c, pull, end = choose(name, 'tif')
        
        if True: # Load from saved 
            bin_im, skel_im, node_im, seg_im = [cv2.imread(f'{impath}{pull}{code}{name}.png', cv2.IMREAD_GRAYSCALE).astype(bool) for code in ims]
            thick_im = 2*dt_edt(bin_im)*seg_im/c
        if False: # Resegment images
            bin_im = process(f'{impath}{pull}{name}{end}')
            node_im, seg_im, skel_im = segmentize(bin_im)
            # save images
    
        bt, nt = thickness(bin_im, node_im, seg_im)
        dil_nodes = expand_nodes(node_im, int(nt))
        short_seg = logical_minus(skel_im, dil_nodes)

        if False: # Node Connectivity 
    
            bt, nt = thickness(bin_im, node_im, seg_im)
            bt, nt = round(bt/c, 6), round(nt, 2)
            
            connections = connectivity(seg_im, node_im, dil_nodes,
                                          return_map=False, return_all=True)
            
            cnt, vals, obj = axeses[m].hist(connections, alpha=0.5,
                                       color=histcolor[colorpick],
                                       label=j[:-1], bins=np.arange(3,10), )
        
        # Length by Width Distribution
        if True: # Characterize Struts
            len_im, seg_table =  shape_analyze(seg_im, minval=1,table=True,
                                               xprops=['bbox'])
                                               
            seg_table = pd.DataFrame(seg_table)
            seg_table['y_comp'] = seg_table['bbox-2'] - seg_table['bbox-0']
            #traits.append([name, np.mean(seg_table.major_axis_length)/c, np.mean(seg_table.y_comp)/c])
            cnt, vals, obj = axeses[m].hist((seg_table['bbox-2']- seg_table['bbox-0'])/c, alpha=0.5,
                                           color=histcolor[colorpick],
                                           label=j[:-1], bins=np.arange(0,10,0.25))
            #avg_length = np.mean(seg_table['major_axis_length'])/c

            #plt.savefig('maxaxis_'+name+'.png', dpi=300, bbox_inches='tight')
        if False:
            # Shape Anisotropy                                                                                                     
            cell_im, cell_table =  shape_analyze(bin_im, bg=1,
                                                 minval=20, table=True)
            cell_table['aniso'] = cell_table['major_axis_length'] / cell_table['minor_axis_length']

            cnt, vals, obj = axeses[m].hist(cell_table['aniso'], alpha=0.5,
                                       color=histcolor[colorpick],
                                       label=j[:-1], bins=np.arange(1,6,0.25))
        if False:
            lthick = np.array(length_thick(seg_im, dt_edt(bin_im), listed=True))/2
            cnt, vals, obj = axeses[m].hist(lthick, alpha=0.5,
                                       color=histcolor[colorpick],
                                       label=j[:-1], bins=np.arange(0,14,1))
        if False:
            newmaxcnt, newmincnt = np.amax(cnt), np.amin(cnt)
            if newmaxcnt > cntmax:
                cntmax = newmaxcnt
            if newmincnt < cntmin:
                cntmin = newmincnt
        colorpick += 1 
        
    m +=1
    row0=False


if True:
    k = 5
    ytick = np.arange(0, 250, 50)
    axeses[0].legend()
    axeses[0].set_xticks([]), axeses[1].set_xticks([])
    axeses[2].set_xlabel('strut max axis length')
    for ax in range(3):
        axeses[ax].set_yticks(ytick)
        axeses[ax].set_ylabel(geom[ax])
    
    #plt.suptitle('Geometries by Y-Strut Length')
    plt.savefig('max_axis_distribution.png', dpi=500, bbox_inches='tight')