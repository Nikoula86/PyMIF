import pickle, glob
import os, sys, tqdm
import matplotlib.pyplot as plt
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
scrnaseq_file = glob.glob(os.path.abspath(os.path.join('..','..','zebrafish_Wagner_data','*_hvg_%02dhpf.h5ad'%hpf_scrna)))[0]

print(hcr_file)
print(mapped_file)

adata_hcr = sc.read_h5ad( hcr_file )
adata_mapped = sc.read_h5ad( mapped_file )
adata_scrna = sc.read_h5ad( scrnaseq_file )

print(adata_hcr)
print(adata_mapped)
print(adata_scrna)

down = 1

points_mapped = adata_mapped.obs[['x','y','z']].to_numpy()[::down]
points_hcr = adata_hcr.obs[['x','y','z']].to_numpy()[::down]
c_mapped = np.array(adata_mapped[:,[missing_gene]].X[:,0].flatten()[::down])
c_hcr = np.array(adata_hcr[:,[missing_gene]].X[:,0].flatten()[::down])
c_scrna = np.array(adata_scrna[:,[missing_gene]].X[:,0].flatten()[::down])

print(np.min(c_mapped), np.max(c_mapped))
print(np.min(c_hcr), np.max(c_hcr))
print(np.min(c_scrna), np.max(c_scrna))

# import matplotlib.pyplot as plt
# fig,ax = plt.subplots(1,3)
# ax[0].hist([c_mapped])
# ax[1].hist([c_hcr])
# ax[2].hist([c_scrna])
# ax[0].set_yscale('log')
# ax[1].set_yscale('log')
# ax[2].set_yscale('log')


print('color hcr')
points_hcr = vedo.Points(points_hcr)
points_hcr.pointSize(6).cmap("hot", c_hcr)
print('color mapped')
points_mapped = vedo.Points(points_mapped)
points_mapped.pointSize(6).cmap("hot", c_mapped)

t1 = vedo.Text2D('HCR pattern', s=2, font='Arial')
t2 = vedo.Text2D('novosparc pattern',s=2, font='Arial')
# setup figure
plt = vedo.Plotter(N=2)

plt.show([t1, points_hcr], at=0)
plt.show([t2, points_mapped], at=1)
# vedo.show([
#     [t1, points_hcr],
#     [t2, points_mapped]
#     ], N=2, axes=0, 
#     # pos=(0,0),
#     # new=True, interactive=True, 
#     # resetcam=True, 
#     # zoom=1., 
#     # azimuth=200,
#     # elevation=90
#     # azimuth=100,
#     # elevation=0
#     )

