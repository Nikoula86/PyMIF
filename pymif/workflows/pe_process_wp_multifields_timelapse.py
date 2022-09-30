# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 15:29:59 2022

@author: gritti
"""

import sys, os
import pandas as pd
from skimage.io import imread
sys.path.append('W:\\people\\gritti\\code\\pymif')
from pymif import pe_opera

exp_folder = 'Y:\\Nick Marschlich\\EMBL_Barcelona\\Imaging\\Opera PE\\P3_Metabolism\\220929_mezzo-GFP_BF_2-DG_injections'
# ff_folder = 'Y:\\Nicola_Gritti\hannah\\20220824_HF_HBMEC-RM-exp2_FF__2022-08-24T15_45_50-Measurement 1'

'''
First, create the csv file for image bookeeping
'''
# pe_opera.xml2csv(exp_folder)
# pe_opera.xml2csv(ff_folder)

'''
Next, create the FF images with ImageJ manually....

Finally, compile experimental images using flat field correction
'''
df = pd.read_csv(os.path.join(exp_folder,'metadata.csv'))

ffs = [
    None,
    None
    # imread('DAPI_ff.tif'),
    # imread('AF488_ff.tif'),
    # imread('AF561_ff.tif'),
    # imread('AF647_ff.tif'),
]

pe_opera.wellplate.compile_conditions_multifields_timelapse(
    path = os.path.join(exp_folder,'Images'), 
    channel_order = [1, 0], 
    luts_name = ['gray', 'green'], 
    df = df,
    ffs = ffs
)
