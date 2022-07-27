import pickle
import os, sys, tqdm
import matplotlib.pyplot as plt
import numpy as np
from skimage.io import imread, imsave
import pandas as pd
import pyclesperanto_prototype as cle
import napari

import vedo

folders = [
#     os.path.join('..','HCR_time_course','hpf08','Panel1','2021-09-20_084304'), # reference for hpf08
#     os.path.join('..','HCR_time_course','hpf08','Panel1','2021-09-20_084954'), # 
#     os.path.join('..','HCR_time_course','hpf08','Panel1','2021-09-20_085804'), # 
    

    # os.path.join('..','HCR_time_course','hpf10','Panel1','2021-09-20_104647'), # reference for hpf10
    # os.path.join('..','HCR_time_course','hpf10','Panel1','2021-09-20_104949'), # 
    # os.path.join('..','HCR_time_course','hpf10','Panel1','2021-09-20_105356'), # 
    
#     os.path.join('..','HCR_time_course','hpf10','Panel2','2021-10-26_072116'),
#     os.path.join('..','HCR_time_course','hpf10','Panel2','2021-10-26_072414'),
    
#     os.path.join('..','HCR_time_course','hpf10','Panel3','2021-10-26_101501'),
#     os.path.join('..','HCR_time_course','hpf10','Panel3','2021-10-26_101721'),
#     os.path.join('..','HCR_time_course','hpf10','Panel3','2021-10-26_102004'),

#     os.path.join('..','HCR_time_course','hpf10','Panel4','2021-10-26_104907'),
#     os.path.join('..','HCR_time_course','hpf10','Panel4','2021-10-26_105128'),


    os.path.join('..','..','HCR_time_course','hpf12','Panel1','2021-09-20_110602'), # reference for hpf12
    os.path.join('..','..','HCR_time_course','hpf12','Panel1','2021-09-20_105648'), # 

#     os.path.join('..','HCR_time_course','hpf12','Panel2','2021-10-26_104331'), 

#     os.path.join('..','HCR_time_course','hpf12','Panel3','2021-10-26_110047'),  
#     os.path.join('..','HCR_time_course','hpf12','Panel3','2021-10-26_110714'),  

    # os.path.join('..','HCR_time_course','hpf12','Panel4','2021-10-26_111407'),  
    # os.path.join('..','HCR_time_course','hpf12','Panel4','2021-10-26_111929'),  

]


markers_dict = {'Panel1':'sp5l','Panel2':'foxa2','Panel3':'ta','Panel4':'prdm1'}
markers = [markers_dict[folder.split(os.sep)[-2]] for folder in folders]
print(markers)

'''
'''
# fluos_dfs_original = [
#     pd.read_csv(os.path.join(i,'cells_segmented','cell_fluo_reg.csv')) for i in folders    
# ]

# fluos_dfs = [i[i.density>50] for i in fluos_dfs_original]
# # fluos_dfs = [i.sample(5000) for i in fluos_dfs]
# n=len(fluos_dfs)


# # define anchor points as average on a grid

# point1 = fluos_dfs[0][['x_unit','y_unit','z_unit']].to_numpy()
# point1_color = fluos_dfs[0][[markers[0]]].to_numpy().flatten()
# point2 = [fluos_dfs[1][['x_unit','y_unit','z_unit']].to_numpy()]
# point2_color = fluos_dfs[1][[markers[1]]].to_numpy().flatten()

# point1_o = fluos_dfs_original[0][['x_unit','y_unit','z_unit']].to_numpy()
# point1_o_color = fluos_dfs_original[0][[markers[0]]].to_numpy()
# point2_o = [fluos_dfs_original[1][['x_unit','y_unit','z_unit']].to_numpy()]
# point2_o_color = fluos_dfs_original[1][[markers[1]]].to_numpy()

# sources = []
# targets =[]

# grid_size = .5
# pbar = tqdm.tqdm(range(20))
# for n in pbar:
#     pbar.set_description("grid_size: %.3f"%grid_size)
# #     print(grid_size)
#     bins = np.arange(-1,1.1,grid_size)
#     grid_size = grid_size*0.95
#     n_min = 20

#     source = []
#     target = []

#     a = []
#     for i in range(len(bins)-1):
#         for j in range(len(bins)-1):
#             for k in range(len(bins)-1):
#                 p_s = point2_o[-1]
#                 s = p_s[(p_s[:,0]>bins[i])&(p_s[:,0]<=bins[i+1])&
#                         (p_s[:,1]>bins[j])&(p_s[:,1]<=bins[j+1])&
#                         (p_s[:,2]>bins[k])&(p_s[:,2]<=bins[k+1])]

#                 p_t = point1_o
#                 t = p_t[(p_t[:,0]>bins[i])&(p_t[:,0]<=bins[i+1])&
#                         (p_t[:,1]>bins[j])&(p_t[:,1]<=bins[j+1])&
#                         (p_t[:,2]>bins[k])&(p_t[:,2]<=bins[k+1])]

#                 if (s.shape[0]>n_min)&(t.shape[0]>n_min):
#                     source.append(np.mean(s,0))
#                     target.append(np.mean(t,0))
#     source = np.array(source)
#     target = np.array(target)

#     sources.append(source)
#     targets.append(target)

#     p2 = vedo.Points(point2_o[-1]).c('green4')
#     newpt2 = p2.clone().warp(source, target, mode='3d', sigma=1)
#     point2_o.append(newpt2.points())

#     p2 = vedo.Points(point2[-1]).c('green4')
#     newpt2 = p2.clone().warp(source, target, mode='3d', sigma=1)
#     point2.append(newpt2.points())

# point2 = np.array(point2)
# point2_o = np.array(point2_o)
# sources = np.array(sources)
# targets = np.array(targets)

# pickle.dump([point1_o,point2_o,sources,targets, point1_o_color, point2_o_color],open("example_warp.pickle","wb"))




'''
'''
point1_o, point2_o, sources, targets, point1_o_color, point2_o_color = pickle.load(open('example_warp.pickle','rb'))

# video with vedo
down = 2

p2 = []
for p in point2_o:
    p2.append(p)
p2 = [p[::down] for p in p2]

p1 = vedo.Points(point1_o[::down]).c('red4').pointSize(6)
p1_c = point1_o_color[::down]

p2 = [vedo.Points(p).c('cyan4').pointSize(6) for p in p2]
p2_c = point2_o_color[::down]

arrows = [vedo.Arrows(s, t, s=2).c('k') for s,t in zip(sources, targets)]

# setup figure
i=-1
vedo.show([
    p1,
    p2[i],
    # arrows[i]
    ], axes=0, pos=(0,0),
    new=True, interactive=False, 
    # resetcam=True, 
    zoom=1., 
    azimuth=200,
    elevation=90
    # azimuth=100,
    # elevation=0
    )

# start video
vd = vedo.Video("movie.mp4", duration=30, backend='ffmpeg')

# show only reference
for i in range(10):
    p = p1.clone().alpha(1.)
    vedo.show([
        p,
        ], axes=0, pos=(0,0), interactive=False, resetcam=False)
    vd.addFrame()
for i in range(10):
    p = p1.clone().alpha(1.)
    p.pointSize(6).cmap("hot", p1_c) # any matplotlib color map name or custom
    t = vedo.Text2D("sp5l - reference fish", s=3, font='Arial')
    a = vedo.show([
        t,
        p,
        ], axes=0, pos=(0,0), interactive=False, resetcam=False)
    vd.addFrame()
    a.remove(t)

# show only source
for i in range(10):
    p = p2[0].clone().alpha(1.)
    vedo.show([
        p,
        ], axes=0, pos=(0,0), interactive=False, resetcam=False)
    vd.addFrame()
for i in range(10):
    p = p2[0].clone().alpha(1.)
    p.pointSize(6).cmap("hot", p2_c) # any matplotlib color map name or custom
    t = vedo.Text2D("sp5l - warped fish", s=3, font='Arial')
    a = vedo.show([
        t,
        p,
        ], axes=0, pos=(0,0), interactive=False, resetcam=False)
    vd.addFrame()
    a.remove(t)

# overlay the two
for i in range(20):
    vedo.show([
        p1.alpha(.3),
        p2[0].alpha(.3),
        ], axes=0, pos=(0,0), interactive=False, resetcam=False)
    vd.addFrame()
for i in range(15):
    vedo.show([
        p1.alpha(0.3),
        p2[0].alpha(0.3),
        arrows[0]
        ], axes=0, pos=(0,0), interactive=False, resetcam=False)
    vd.addFrame()
for i in range(15):
    vedo.show([
        p1.alpha(.3),
        p2[1].alpha(.3)
        ], axes=0, pos=(0,0), interactive=False, resetcam=False)
    vd.addFrame()
for i in range(10):
    vedo.show([
        p1.alpha(0.3),
        p2[1].alpha(0.3),
        arrows[1]
        ], axes=0, pos=(0,0), interactive=False, resetcam=False)
    vd.addFrame() 
for i in range(10):
    vedo.show([
        p1.alpha(.3),
        p2[2].alpha(.3)
        ], axes=0, pos=(0,0), interactive=False, resetcam=False)
    vd.addFrame()

# show iterations   
for i in range(len(p2)-2):
    vedo.show([
        p1.alpha(.3),
        p2[i+2].alpha(.3)
        ], axes=0, pos=(0,0), interactive=False, resetcam=False)
    vd.addFrame()
for i in range(10):
    vedo.show([
        p1.alpha(.3),
        p2[-1].alpha(.3),
        ], axes=0, pos=(0,0), interactive=False, resetcam=False)
    vd.addFrame()    

# show marker
for i in range(10):
    p = p1.clone().alpha(1.)
    p.pointSize(6).cmap("hot", p1_c) # any matplotlib color map name or custom
    t = vedo.Text2D("sp5l - reference fish", s=3, font='Arial')
    a = vedo.show([
        t,
        p.alpha(.5),
        # p2[-1].alpha(.1),
        ], axes=0, pos=(0,0), interactive=False, resetcam=False)
    vd.addFrame()    
    a.remove(t)
for i in range(10):
    p = p2[-1].clone().alpha(1.)
    p.pointSize(6).cmap("hot", p2_c) # any matplotlib color map name or custom
    t = vedo.Text2D("sp5l - warped fish", s=3, font='Arial')
    a = vedo.show([
        t,
        # p1.alpha(.1),
        p.alpha(.5),
        ], axes=0, pos=(0,0), interactive=False, resetcam=False)
    vd.addFrame()    
    a.remove(t)
vd.close()
