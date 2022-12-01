

import numpy as np
from skimage.io import imread, imsave
import glob, os,  string, tqdm, time
import pandas as pd
from ...imagej_funs._make_lut import make_lut
from ...imagej_funs._imagej_metadata_tags import imagej_metadata_tags

def compile_conditions(path, 
                        conditions, 
                        channel_order, 
                        luts_name,
                        df,
                        ffs):

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
    ffs = [ff/np.median(ff) if ff is not None else 1. for ff in ffs]

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
        
        outpath = os.path.join(os.path.split(path)[0],'compiled',cond)
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
                
                stack_ch = np.stack([imread(os.path.join(path,img_file))/ffs[k] for img_file in df_pos_ch.filename])
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


##########################################################################################

if __name__ == '__main__':

    paths = [
                os.path.join(
                    '2021-07-21_BraGFP_150-300-450-cells_inhibSB43_120h',
                    'Images'
                    ),
            ]

    # this is to, for instance, arrange BF in first channel and GFP in second channel
    # available LUTS: gray, red, green, blue, magenta, cyan, yellow
    channel_orders = [
                        [1,0],
                    ]
    luts_names = [ 
                    ['gray', 'green'],
                ]

    # build the conditions map for every well
    c1 = [
            ['init_150cells_inhib_48-72' for i in range(4)]+['init_150cells_inhib_72-96' for i in range(4)],
            ['init_150cells_inhib_48-72' for i in range(4)]+['init_150cells_inhib_72-96' for i in range(4)],
            ['init_150cells_inhib_48-72' for i in range(4)]+['init_150cells_inhib_72-96' for i in range(4)],
            ['init_150cells_inhib_48-72' for i in range(4)]+['init_150cells_inhib_72-96' for i in range(4)],
            ['init_300cells_inhib_48-72' for i in range(4)]+['init_300cells_inhib_72-96' for i in range(4)],
            ['init_300cells_inhib_48-72' for i in range(4)]+['init_300cells_inhib_72-96' for i in range(4)],
            ['init_300cells_inhib_48-72' for i in range(4)]+['init_300cells_inhib_72-96' for i in range(4)],
            ['init_300cells_inhib_48-72' for i in range(4)]+['init_300cells_inhib_72-96' for i in range(4)],
            ['init_450cells_inhib_48-72' for i in range(4)]+['init_450cells_inhib_72-96' for i in range(4)],
            ['init_450cells_inhib_48-72' for i in range(4)]+['init_450cells_inhib_72-96' for i in range(4)],
            ['init_450cells_inhib_48-72' for i in range(4)]+['init_450cells_inhib_72-96' for i in range(4)],
            ['init_450cells_inhib_48-72' for i in range(4)]+['init_450cells_inhib_72-96' for i in range(4)],
        ]

    conditions_all = [
            # c1,
            c1,
        ]

    for i in range(len(paths)):
        path = paths[i]
        print(path)
        channel_order = channel_orders[i]
        luts_name = luts_names[i]
        conditions = conditions_all[i]
        
        time.sleep(1)

        compile_condition(path, conditions, channel_order, luts_name)
