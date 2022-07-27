# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 19:43:26 2021

@author: nicol
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 09:48:33 2021

@author: nicol
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 09:58:28 2021

@author: nicol
"""

import pickle
import os, sys, tqdm
import matplotlib.pyplot as plt
import numpy as np
from skimage.io import imread, imsave
import pandas as pd
import pyclesperanto_prototype as cle
import napari
import scanpy as sc

import vedo


'''
TODO:
    
    COLOR GRAY CELLS WITH THE GOPI GERM LAYER COLOR!!!
'''
df = pd.read_csv(os.path.join('..','live_mapped.csv'))

# remove outlier cells
df = df[np.sqrt(df.x**2+df.y**2+df.z**2)<1.2]
# tp=12.
# df = df[df.t==tp]
df.TissueName = [str(i) for i in tqdm.tqdm(df.TissueName)]

c1 = df.germ_layer

df['germ_name'] = [['Ectoderm','Mesendoderm','Endoderm'][i-1] for i in c1]

newname = []
for gn, tn in zip(df.germ_name, df.TissueName):
    if tn == 'nan':
        newname.append(gn)
    else:
        newname.append(tn)
df['TissueName_corr'] = newname

tissuenames = ['Endoderm', 'Mesoderm', 'Epidermal', 'Forebrain / Optic',
               'Hindbrain / Spinal Cord', 'Midbrain', 'Neural Crest', 'Mesendoderm','Ectoderm']
# c2 = [ tissuenames.index(i) for i in df.TissueName ]

palette1 = ['cyan', 'red', 'yellow',] #ecto,meso,endo

palette2=[
    '#CC7722',
    '#FF0000',
    '#8B0000',
    '#008080',
    '#9467bd',
    '#8c564b',
    '#e377c2',
    'red',
    'cyan'
]

# color = [palette[i] for i in df.TissueName]

# video with vedo

t_idxs = list(set(df.t_idx))
t_idxs.sort()
elevation = 0
azimuth = 0


###########################
'''
make stacked barplot
'''
cells_tl = df
fig, ax = plt.subplots()

t_idxs = list(set(cells_tl.t_idx))

palette2_hist=[
    '#CC7722',
    '#FF0000',
    '#FF0000',
    '#8B0000',
    '#008080',
    '#9467bd',
    '#8c564b',
    '#e377c2',
    'cyan'
]
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])
lighten_degree = [0.7,0.7,0.7,1.,0.7,0.7,0.7,0.7,0.7]
palette2_hist = [lighten_color(i,j) for i,j in zip(palette2_hist, lighten_degree)]

tns = ['Endoderm', 'Mesoderm', 'Mesendoderm', 'Epidermal', 'Forebrain / Optic',
               'Hindbrain / Spinal Cord', 'Midbrain', 'Neural Crest','Ectoderm']

for t_idx in tqdm.tqdm(range(150, 259)):
    cells = cells_tl[cells_tl.t_idx == t_idx]
    comp_tns = []
    for tn in tns:    
        cells_tn = cells[cells.TissueName_corr==tn]
        comp_tns.append(len(cells_tn))
    comp_tns = np.array(comp_tns)/np.sum(comp_tns)
    ##################
    bottom = 0
    for i in range(len(tns)):
        # if t_idx==150:
        #     print(palette2_hist[i])
        ax.bar([t_idx], comp_tns[i], bottom=bottom, label=tns[i], 
               color = palette2_hist[i], width=1, alpha=1., edgecolor = "none")
        bottom += comp_tns[i]
        
ax.set_xlim((150,258))
ax.set_ylim(0,1)
ax.set_xticks([168,189,210,231,252])
ax.set_xticklabels(2*np.array([168,189,210,231,252])/42)
# ax.legend()

###################################
'''
make video
'''

# vd = vedo.Video("movie_final.mp4", duration=15, backend='ffmpeg')
# plt = vedo.Plotter(N=1)

# t_min=150
# ###################
# df_tp = df[df.t_idx==t_min]
     
# # c1_tp = df_tp.germ_layer
# # points1 = vedo.Points( df_tp[['x','y','z']].to_numpy() )
# # points1.alpha(0.4).pointSize(20)
# # points1.cmap(palette1, c1_tp)
# # elevation1=0
# # if t_idx == 150:
# #     elevation1=90
# # plt.show([points1], at=0, bg="black", azimuth=azimuth, elevation=elevation1)
 
# palette = []
# tnames = []
# for color, tname in zip(palette2,tissuenames):
#     if tname in list(df_tp.TissueName_corr):
#         palette.append(color)
#         tnames.append(tname)
 
# c2_tp = [ tnames.index(i) for i in df_tp.TissueName_corr ]   
# # print(s)
 
# points2 = vedo.Points( df_tp[['x','y','z']].to_numpy() )
# points2.alpha(0.4).pointSize(20)
# points2.cmap(palette, c2_tp)
# plt.show([points2], at=0, bg="black", 
#       azimuth=60, 
#       elevation=90, 
#       interactive=False)
# #################################


# for t_idx in np.arange(t_min,260):
 
#     df_tp = df[df.t_idx==t_idx]
     
#     # c1_tp = df_tp.germ_layer
#     # points1 = vedo.Points( df_tp[['x','y','z']].to_numpy() )
#     # points1.alpha(0.4).pointSize(20)
#     # points1.cmap(palette1, c1_tp)
#     # elevation1=0
#     # if t_idx == 150:
#     #     elevation1=90
#     # plt.show([points1], at=0, bg="black", azimuth=azimuth, elevation=elevation1)
     
#     palette = []
#     tnames = []
#     for color, tname in zip(palette2,tissuenames):
#         if tname in list(df_tp.TissueName_corr):
#             palette.append(color)
#             tnames.append(tname)
     
#     c2_tp = [ tnames.index(i) for i in df_tp.TissueName_corr ]   
#     # print(s)
     
#     elevation=0
#     azimuth=0
#     if t_idx == t_min:
#         elevation=0
#         azimuth=-145
#     points2 = vedo.Points( df_tp[['x','y','z']].to_numpy() )
#     points2.alpha(0.8).pointSize(20)
#     points2.cmap(palette, c2_tp)
#     plt.show([points2], at=0, bg="black", 
#              axes=0,
#       azimuth=azimuth, 
#       elevation=elevation, 
#       interactive=False, zoom=1.)
    
#     vd.addFrame()
    
# vd.close()
