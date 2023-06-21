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
    

headers = ['name', 'beam_thickness', 'length_average', 
          'edge_connectivity','shape_anisotropy']
traits = []

fig, axes = plt.subplots(3,3, figsize=(12.5,10))
    
m = 0
axeses = axes.flatten()
def y_comp(r):
    y1, x1, y2, x2 = r.bbox
    return y2-y1

def aniso(r):
    return r['major_axis_length']/r['minor_axis_length']

def orient(r):
    return np.degrees(r['orientation'])

row0 = True
df = pd.DataFrame(columns=['name', 'cells', 'struct. area', 'neg. area', 'nodes'])
for i in geom:
    for j in perturb:
        if row0: # Set Top Titles 
            axeses[m].set_title(j[:-1])
        
        name = j+i
        print(name)
        c, pull, end = choose(name, 'tif') # Get info 
    
        if False: # Resegment images
            bin_im = process(f'{impath}{pull}{name}{end}')
            node_im, seg_im, skel_im = segmentize(bin_im)
            #cv2.imwrite(f'{impath}node_{name}.png', 255*node_im)
            #cv2.imwrite(f'{impath}sk_{name}.png', 255*skel_im)
            #cv2.imwrite(f'{impath}bin_{name}.png', 255*bin_im)
            #cv2.imwrite(f'{impath}seg_{name}.png', 255*seg_im)
            
        if True: # Load from saved 
            bin_im, skel_im, node_im, seg_im = [cv2.imread(f'{impath}{pull}{code}_{name}.png', cv2.IMREAD_GRAYSCALE).astype(bool) for code in ['bin', 'skel', 'node', 'seg']]

        bt, nt = thickness(bin_im, node_im, seg_im)
        dil_nodes = expand_nodes(node_im, int(nt))
        short_seg = logical_minus(seg_im, dil_nodes)
        #f = see_nodes(bin_im, dil_nodes, node_im, seg_im)
        #cv2.imwrite(f'see_nodes_{name}.png', 255*f)
        
        if False: # Node Connectivity 
            dil_nodes = expand_nodes(node_im, int(nt))
            node_connectivity, connect_map = connectivity(seg_im, node_im, dil_nodes,
                                          return_map=True, return_all=False)

            # Node connectivity visualization
            f = node_connectivity_map(bin_im, connect_map)
            plt.savefig('labnodes_'+name+'.png', dpi=300, bbox_inches='tight')
            
        if False: # Characterize Struts
            len_im, seg_table =  shape_analyze(seg_im, minval=1,table=False)
            
            colored, feats = color_by(len_im, seg_table, 
                                      feat='major_axis_length')
            #axeses[m], im  = color_struts(colored/c,  ax = axeses[m])
            f = color_struts(colored/c)
            plt.savefig(f'maxaxis_{name}.png', dpi=200, bbox_inches='tight')
            
        if True:
            # Shape Anisotropy                                                                                                     
            cell_im, cell_table =  shape_analyze(bin_im, bg=1,
                                                 minval=20, table=True)
            nodes, _ = cv2.connectedComponents(dil_nodes)
            cell_table = pd.DataFrame(cell_table)
            df.loc[len(df.index)] = [name, len(cell_table), 
                                     round(np.sum(bin_im)/(c**2),2), 
                                     round(np.sum(cell_table.area)/(c**2),2),
                                     nodes] 
            #colored, _ = color_by(cell_im, cell_table, func=aniso)
            #axeses[m], im   = color_struts(colored, isax = axeses[m],
            #                               struts=False, bounds=(1,2))
            #f = color_struts(colored, struts=False)
            #plt.savefig(f'area_{name}.png', bbox_inches='tight', dpi=300)
        if False:
            len_thick_im = length_thick(seg_im, dt_edt(bin_im), listed=False)
            f = color_struts(len_thick_im/2, struts=False)
            plt.savefig(f'length_thick_{name}.png', dpi=300, bbox_inches='tight')
        m +=1
    row0=False
    
if False:
    fig.subplots_adjust(right=0.8, wspace=0.05, hspace=0.025)
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.8)
    for idx, a in enumerate([0,3,6]):
        axeses[a].set_ylabel(geom[idx])
    
    #plt.suptitle('Geometries by Y-Strut Length')
    plt.savefig('length_by_thickness_grid.png', dpi=500, bbox_inches='tight')