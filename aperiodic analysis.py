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

path = '/Users/nikkivanhandel/Projects/Aperiodic Lattices/imgs/'


# Read scale values for tif images 
tif_scale = pd.read_csv('tif_scale.csv')
px_per_cm = 213 # for pngs (manually read scalebar)

# Lazy way to iterate through files
perturb = ['V0_N0_', 'V.2_N0_', 'V0_N.2_', ]
geom = ['H','S','T']

pull = path+ 'BlueLight/' # Pull files from here 
push = pull               # Push files to here   

headers = ['name', 'beam_thickness', 'length_average', 
          'edge_connectivity','shape_anisotropy']
traits = []

for i in geom:
    for j in perturb:
        break
    name = j+i
    c= float(tif_scale.conversion.loc[tif_scale.geom==name])
    cell_im = process(pull+name+'.TIF')
    print(name)
    node_im, segment_im = segmentize(cell_im)
    
    fig = draw_min_maj(np.logical_not(cell_im))
    plt.savefig(f'aniso_{name}.png', bbox_inches='tight', 
                dpi=300)
    
    if False:
        # Segment and Node Thickness
        bt, nt = thickness(cell_im, node_im, segment_im)
        bt, nt = round(bt/c, 6), round(nt, 2)
        
        dil_nodes = expand_nodes(node_im, int(nt))
        # Determine Node Connectivity 
        con = connectivity(segment_im, node_im, dil_nodes)
        avg_connect = np.mean(con)
        # Length of Struts
        _, seg_table =  shape_analyze(segment_im, minval=1,
                                               table=True)
        avg_length = np.mean(seg_table['major_axis_length'])/c
        
        # Shape Anisotropy 
        _, cell_table =  shape_analyze(cell_im, bg=1, 
                                       minval=20, table=True)
        cell_table['anisotropy'] = cell_table['major_axis_length'] / cell_table['minor_axis_length']
        anisotropy = np.mean(cell_table['anisotropy'] )
        
        traits.append([name, bt, avg_length, avg_connect,anisotropy])
    if False:
        # Just storing this for later
        colored, feats = color_by(label_im, props, feat='major_axis_length')
        fig = color_struts(colored/c)
        
        stack = see_nodes(cell_im, dil_nodes, node_im, segment_im)
        #cv2.imwrite(f'seenodes_{name}.png', 255*stack)
        fig = node_connectivity_map(cell_im, c)
        #plt.savefig('labnodes_'+name+'.png', dpi=300)

#df = pd.DataFrame(traits, columns=headers)
#df.set_index('name', inplace=True)
#plt.imshow(c2(segment_im, node_im, dt_edt(cell_im)))
