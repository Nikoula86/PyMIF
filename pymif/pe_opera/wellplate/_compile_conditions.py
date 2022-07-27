

import numpy as np
from skimage.io import imread, imsave
import glob, os,  string, tqdm, time
from ...imagej_funs._make_lut import make_lut
from ...imagej_funs._imagej_metadata_tags import imagej_metadata_tags

def compile_conditions(path, conditions, channel_order, luts_name):
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
    # find all tiff files in the folder
    flist = glob.glob(os.path.join(path,'*.tiff'))
    flist.sort()

    # find out all positions (the first 6 characters, e.g.: r01c01 )
    pos = list(set( [os.path.split(f)[-1][:6] for f in flist] ))
    pos.sort()

    # define well id to convert e.g. r01c01 into A01
    d = dict(enumerate(string.ascii_uppercase, 1))

    pbar = tqdm.tqdm(pos)
    for p in pbar:
        well = p[4:6]+d[int(p[1:3])]
        cond = conditions[int(p[4:6])-1][int(p[1:3])-1]

        pbar.set_description(p + ' ' + well + ' ' + cond)
        pbar.update()
        
        outpath = os.path.join(os.path.split(path)[0],'compiled',cond)
        if not os.path.exists(outpath):
            os.makedirs(outpath)

        # extract all files from this well
        flist = glob.glob(os.path.join(path,p+'*.tiff'))
        flist.sort()

        # extract all files from this timepoint
        channels_list = glob.glob(os.path.join(path,p+'*'+'sk*fk*.tiff'))
        channels_list.sort()

        # find channels
        channels = list(set([f[f.index('-ch')+3:f.index('-ch')+4] for f in channels_list]))
        channels = np.array([int(ch) for ch in channels])
        channels.sort()

        stacks = []
        for channel in channels:
            # extract all files from this channel
            stack_list = glob.glob(os.path.join(path,p+'*-ch'+str(channel)+'*sk*fk*.tiff'))
            stack_list.sort()

            stack = []
            for f in stack_list:
                stack.append(imread(f))
            stack = np.array(stack)
            if stack.shape[0]==1:
                stack = stack[0]
            stacks.append(stack)

        stacks = np.array(stacks).astype(np.uint16)

        # order channels
        stacks = np.array([stacks[ch] for ch in channel_order]).astype(np.uint16)
        if stacks.ndim == 4:
            stacks = np.moveaxis(stacks,1,0)
        # print(stacks.shape)

        # create imagej metadata with LUTs
        luts_dict = make_lut(luts_name)
        ijtags = imagej_metadata_tags({'LUTs': [luts_dict[i] for i in luts_name]}, '>')

        outname = well+'.tif'
        imsave(os.path.join(outpath,outname),stacks, byteorder='>', imagej=True,
                        metadata={'mode': 'composite'}, extratags=ijtags)


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
