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
adata = sc.read_h5ad(os.path.join('..','..','HCR_time_course','hpf10_synthetic_HCR.h5'))

genes = list(adata.var_names)
# video with vedo
down = 1

points = vedo.Points( adata.obs[['x','y','z']].to_numpy()[::down] )
points.alpha(.5).pointSize(10)
cam = dict(pos=(5.598, -3.741, -2.091),
           focalPoint=(0.01097, 0.05061, 0.01802),
           viewup=(0.6059, 0.6065, 0.5148),
           distance=7.074,
           clippingRange=(3.581, 11.49))
# setup figure
i=-1
vedo.show([
    points,
    ], axes=4, pos=(0,0),
    camera=cam
    # new=False, interactive=False, 
    # # resetcam=True, 
    # zoom=1.2, 
    # # azimuth=200,
    # # elevation=90
    # azimuth=0,
    # elevation=45,
    # # roll=-60,
    # title='my title'
    )

# start video
vd = vedo.Video("movie2.mp4", duration=30, backend='ffmpeg')

# show only reference
t = vedo.Text2D(s=5, font='Arial')
j=0
for g in genes:
    for i in range(26):
        c = adata[:,[g]].X[::down,0]
        percs = np.percentile(c, (50,99))
        c = np.clip(c, percs[0], percs[1])
        points.cmap("hot", c) # any matplotlib color map name or custom
        t.text(g)
        a = vedo.show([
                t,
                points,
            ], 
            new=False,
            axes=4, pos=(0,0), 
            interactive=False, resetcam=True,
            zoom=1.2,
            azimuth=3,
            elevation=0,
            # text=g
            )
        
        vd.addFrame()
        j+=1

vd.close()
