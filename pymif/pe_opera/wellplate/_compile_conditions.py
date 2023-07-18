

import numpy as np
from skimage.io import imread, imsave
import glob, os,  string, tqdm, time
import pandas as pd
from ...imagej_funs._make_lut import make_lut
from ...imagej_funs._imagej_metadata_tags import imagej_metadata_tags
from ..pe_io._extract_ffc_info import extract_ffc_info

def compile_conditions(
        path, 
        conditions, 
        channel_order, 
        luts_name,
        df,
        ffs,
        ff_mode = 'PE', 
        outfolder = 'compiled',
        image_folder = os.path.join("Images"),
        ):

    '''This function combines images of a 96WP acquired by PE.

    Parameters
    ----------
    path: string
            a string containing the path t the experiment folder. Has to point to the "Images" folder
    conditions: iterable, string
            list of conditions for every well (shape 12x8)
    channels_order: iterable, int
            the order in which the channels have to be arranged in the output images.
            E.g. if PE acquires GFP first and BF next, and you want BF-GFP, then the order is [1,0]
    luts_name: iterable, string
            list of colors to be used to show the channels

    Returns
    -------
    a "compiled" folder, containing the conditions subfolders, containing 1 multichannel tif file for every well in the experiment

    NOTE:
    This script assume the experiment contains just one FOV per well!
    '''
    # ff_mode: 'PE' for PE FF correction, use 'slide' for autofluorescence slide, use 'None' for no correction
    ffs = [1. for ff in ffs]
    if ff_mode == 'slide':
        ffs = [ff/np.median(ff) if ff is not None else 1. for ff in ffs]
    elif ff_mode == 'PE':
        ffs_info = extract_ffc_info(path, channel_order)
        ffs = [ff_info['ff_profile'] for ff_info in ffs_info]

    # find out all wells
    wells = df.groupby(['row','col']).size().reset_index()    

    # define well id to convert e.g. r01c01 into A01
    d = dict(enumerate(string.ascii_uppercase, 1))


    pbar = tqdm.tqdm(wells.iterrows())
    conversion = pd.DataFrame({})
    for p in pbar:
        # print(p)
        r = int(p[1].row)
        c = int(p[1].col)
        well = d[r]+'%02d'%c

        
        cond = conditions[int(p[1].col)-1][int(p[1].row)-1]
        
        pbar.set_description(well + ' ' + cond)
        pbar.update()
        
        outpath = os.path.join(path,outfolder,cond)
        if not os.path.exists(outpath):
            os.makedirs(outpath)    

        df_well = df[(df.row==r)&(df.col==c)]
        
        if len(df_well)>0:
            
            stack = []
            for k, ch in enumerate(channel_order):
                df_pos_ch = df_well[df_well.channel==(ch+1)]
                df_pos_ch = df_pos_ch.sort_values(by='Zpos')
                
                # print('-'*25,'ch:',ch)
                # print(df_pos_ch)
                
                stack_ch = np.stack([imread(os.path.join(path,image_folder,img_file))/ffs[k] for img_file in df_pos_ch.filename])
                stack.append(stack_ch)

            # order channels
            stacks = np.array(stack).astype(np.uint16)
            stacks = np.swapaxes(stacks, 0, 1)

            # create imagej metadata with LUTs
            luts_dict = make_lut(luts_name)
            # luts_dict = make_lut_old()
            ijtags = imagej_metadata_tags({'LUTs': [luts_dict[lut_name] for lut_name in luts_name]}, '>')

            outname = well+'.tif'

            raw = pd.DataFrame({'filename':[outname],
                                'row_idx':[r],
                                'col_idx':[c]})
            conversion = pd.concat([conversion,raw], ignore_index=True)
            
            imsave(os.path.join(outpath,outname),stacks, byteorder='>', imagej=True,
                            metadata={'mode': 'composite'}, extratags=ijtags)

    conversion.to_csv(os.path.join(outpath, 'metadata.csv'))
