# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 15:23:19 2022

@author: gritti
"""

import numpy as np
from skimage.io import imread, imsave
import os, string, tqdm
import pandas as pd
from ...imagej_funs._make_lut import make_lut
from ...imagej_funs._imagej_metadata_tags import imagej_metadata_tags

def compile_conditions_multifields_timelapse(
        path, 
        channel_order, 
        luts_name, 
        df,
        ffs
    ):
    
    ffs = [ff/np.median(ff) if ff is not None else 1. for ff in ffs]

    # find out all wells
    wells = df.groupby(['row','col']).size().reset_index()

    # define well id to convert e.g. r01c01 into A01
    d = dict(enumerate(string.ascii_uppercase, 1))

    pbar = tqdm.tqdm(wells.iterrows())
    for i, p in pbar:
        r = int(p.row)
        c = int(p.col)
        well = d[r]+'%02d'%c
        
        conversion = pd.DataFrame({})

        pbar.set_description(well)
        pbar.update()
        
        outpath = os.path.join(os.path.split(path)[0], 'compiled', well)
        if not os.path.exists(outpath):
            os.makedirs(outpath)
            
        df_well = df[(df.row==r)&(df.col==c)]
        
        timepoints = list(set(df_well.timepoint))
        timepoints.sort()
        
        for timepoint in timepoints:
            df_tp = df_well[df_well.timepoint==timepoint]
        
            # find all fields inside this well
            fields = df_tp.groupby(['Ypos','Xpos']).size().reset_index()
            fields = fields.sort_values(by=['Ypos','Xpos'])
            l = list(set(fields.Ypos))
            l.sort()
            fields['Yidx'] = [l.index(v) for v in fields.Ypos]
            l = list(set(fields.Xpos))
            l.sort()
            fields['Xidx'] = [l.index(v) for v in fields.Xpos]
            
            for j, f in tqdm.tqdm(fields.iterrows(), total = len(fields)):
                x = f.Xpos
                y = f.Ypos
                xidx = f.Xidx
                yidx = f.Yidx
                
                df_pos = df_tp[(df_tp.Xpos==x)&(df_tp.Ypos==y)]
                
                # print('-'*50)
                # print(df_pos)
                
                if len(df_pos)>0:
    
                    # print('Images foudn')
                    stack = []
                    for k, ch in enumerate(channel_order):
                        df_pos_ch = df_pos[df_pos.channel==(ch+1)]
                        df_pos_ch = df_pos_ch.sort_values(by='Zpos')
                        print('-'*25,'x:',xidx,'y:',yidx,'ch:',ch)
                        print(df_pos_ch)
                        # [print(img_file) for img_file in df_pos_ch.filename]
                        # print([os.path.join(folder_raw,exp_folder,'Images',img_file) for img_file in df_pos_ch.filename])
                        stack_ch = np.stack([imread(os.path.join(path,img_file))/ffs[k] for img_file in df_pos_ch.filename])
                        stack.append(stack_ch)
        
                    # order channels
                    stacks = np.array(stack).astype(np.uint16)
                    stacks = np.swapaxes(stacks, 0, 1)
        
                    # create imagej metadata with LUTs
                    luts_dict = make_lut(luts_name)
                    # luts_dict = make_lut_old()
                    ijtags = imagej_metadata_tags({'LUTs': [luts_dict[lut_name] for lut_name in luts_name]}, '>')
                    
                    outname = 'field%03d_tp%05d.tif'%(j,timepoint)
        
                    raw = pd.DataFrame({'tile_idx':[j],
                                        'filename':[outname],
                                        'row_idx':[yidx],
                                        'col_idx':[xidx],
                                        'timepoint':[timepoint]})
                    conversion = pd.concat([conversion,raw], ignore_index=True)
        
                    # print(outname)
                    imsave(os.path.join(outpath,outname),stacks, byteorder='>', imagej=True,
                                    metadata={'mode': 'composite'}, extratags=ijtags, check_contrast=False)
                
        conversion.to_csv(os.path.join(outpath, 'metadata.csv'))
            


