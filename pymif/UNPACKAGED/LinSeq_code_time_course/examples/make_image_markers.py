# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 16:24:19 2021

@author: nicol
"""

import pickle, glob
import os, sys, tqdm
import numpy as np
from skimage.io import imread, imsave
import pandas as pd
import pyclesperanto_prototype as cle
import napari
import scanpy as sc

import vedo

hpf = 12
hpf_scrna = hpf
if hpf == 12:
    hpf_scrna = 14
num_locs_cells=4000
alpha_str = str(1)
tol=1e-15
epsilon=1e-4
num_neighbors_source = 10
num_neighbors_target = 10
missing_gene = 'ta'


hcr_file = os.path.join('..','..','HCR_time_course','hpf'+str(hpf)+'_synthetic_HCR.h5')
mapped_file = os.path.join('..','..','HCR_time_course','mapped_hpf'+str(hpf)+\
                           '_nlocs-'+str(num_locs_cells)+'_alpha-'+alpha_str+\
                           '_tol-'+str(tol)+'_eps-'+str(epsilon)+\
                           '_nS-'+str(num_neighbors_source)+\
                           '_nT-'+str(num_neighbors_target)+
                           '_NO-'+missing_gene+'.h5')
scrnaseq_file = glob.glob(os.path.join('..','..','zebrafish_Wagner_data','*_hvg_%02dhpf.h5ad'%hpf_scrna))[0]

print(hcr_file)
print(mapped_file)

adata_hcr = sc.read_h5ad( hcr_file )

genes = list(adata_hcr.var_names)

print(adata_hcr)

down = 2

points = adata_hcr.obs[['x','y','z']].to_numpy()[::down]
points_hcr = []
c_hcr = []
for g in genes:
    c = np.array(adata_hcr[:,[g]].X[:,0].flatten()[::down])
    percs = np.percentile(c, (10,99))
    c = np.clip(c, percs[0], percs[1])
    c_hcr.append( c )
    p = vedo.Points(points).clone()
    p.pointSize(6).cmap("hot", c).alpha(0.9)
    points_hcr.append(p)

texts = [vedo.Text2D(s=2, font='Arial') for g in genes]
texts = [t.text(g) for t, g in zip(texts,genes)]

info = [[t,p] for t,p in zip(texts, points_hcr)]

# setup figure
plt = vedo.Plotter(N=len(info))
for i in range(len(info)):
    elevation = 0
    azimuth = 0
    if i == 0:
        elevation = 90
        azimuth = 130
    plt.show(info[i][0],info[i][1], at=i, azimuth=azimuth, elevation=elevation)

    
