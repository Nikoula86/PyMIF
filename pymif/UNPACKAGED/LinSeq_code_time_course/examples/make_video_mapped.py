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
adata = sc.read_h5ad(os.path.join('..','..','HCR_time_course','mapped_transfer_hpf12.h5ad'))


df = adata.obs
tissuenames = ['Endoderm', 'Mesoderm', 'Epidermal', 'Forebrain / Optic',
               'Hindbrain / Spinal Cord', 'Midbrain', 'Neural Crest']
c = [ tissuenames.index(i) for i in df.TissueName ]

print(tissuenames)
print(c)

# ax.imshow(img_gauss[z])
palette=[
    '#CC7722',
    '#FF0000',
    '#8B0000',
    '#008080',
    '#9467bd',
    '#8c564b',
    '#e377c2',
]
# color = [palette[i] for i in df.TissueName]

genes = list(adata.var_names)
# video with vedo
down = 1

points = vedo.Points( adata.obs[['x','y','z']].to_numpy()[::down] )
points.alpha(.5).pointSize(10)
points.cmap(palette, c)
# cam = dict(pos=(5.598, -3.741, -2.091),
#             focalPoint=(0.01097, 0.05061, 0.01802),
#             viewup=(0.6059, 0.6065, 0.5148),
#             distance=7.074,
#             clippingRange=(3.581, 11.49))

vd = vedo.Video("mapped.mp4", duration=15, backend='ffmpeg')
for i in range(360):
    elevation = 0
    azimuth = 1
    if i == 0:
        elevation = 90
        azimuth = 0
        vedo.show([points], bg="black", new=False, interactive=False, azimuth=azimuth, elevation=elevation)
    else:
        vedo.show([points], bg="black", new=False, interactive=False, azimuth=azimuth, elevation=elevation)
        vd.addFrame()
    
vd.close()

    
# # setup figure
# i=-1
# vedo.show([
#     points,
#     ], axes=4, pos=(0,0),
#     camera=cam
#     # new=False, interactive=False, 
#     # # resetcam=True, 
#     # zoom=1.2, 
#     # # azimuth=200,
#     # # elevation=90
#     # azimuth=0,
#     # elevation=45,
#     # # roll=-60,
#     # title='my title'
#     )

# # plt = vedo.Plotter()
# # plt.show([points], camera=cam)
# # # start video
# # vd = vedo.Video("mapped.mp4", duration=30, backend='ffmpeg')

# # # show only reference
# # for i in range(120):
# #     plt.show([
# #             points,
# #         ], 
# #         # new=False,
# #         axes=4, #pos=(0,0), 
# #         interactive=False, resetcam=True,
# #         zoom=1.2,
# #         azimuth=3,
# #         elevation=0,
# #         # text=g
# #         )
    
# #     vd.addFrame()

# # vd.close()
