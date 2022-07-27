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
'''
df = pd.read_csv(os.path.join('..','live_mapped.csv'))
tp=12.
df = df[df.t==tp]
df.TissueName = [str(i) for i in tqdm.tqdm(df.TissueName)]


c1 = df.germ_layer

tissuenames = ['Endoderm', 'Mesoderm', 'Epidermal', 'Forebrain / Optic',
               'Hindbrain / Spinal Cord', 'Midbrain', 'Neural Crest', 'nan']
c2 = [ tissuenames.index(i) for i in df.TissueName ]

palette1 = ['cyan', 'red', 'yellow',]

palette2=[
    '#CC7722',
    '#FF0000',
    '#8B0000',
    '#008080',
    '#9467bd',
    '#8c564b',
    '#e377c2',
    'k'
]

# color = [palette[i] for i in df.TissueName]

# video with vedo

t_idxs = list(set(df.t_idx))
t_idxs.sort()
elevation = 0
azimuth = 0

plt = vedo.Plotter(N=2)

points1 = vedo.Points( df[['x','y','z']].to_numpy() )
points1.alpha(0.75).pointSize(20)
points1.cmap(palette1, c1)
plt.show([points1], at=0, bg="black", azimuth=azimuth, elevation=90)


points2 = vedo.Points( df[['x','y','z']].to_numpy() )
points2.alpha(0.75).pointSize(20)
points2.cmap(palette2, c2)
plt.show([points2], at=1, bg="black", azimuth=azimuth, elevation=elevation)
