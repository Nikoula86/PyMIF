import sys, os
import pandas as pd
from skimage.io import imread
sys.path.append('Z:\\people\\gritti\\code\\pymif')
import pymif

exp_folder = 'Y:\\Nicola_Gritti\hannah\\20220824_HF_HBMEC-RM-exp2__2022-08-24T15_02_41-Measurement 1'
ff_folder = 'Y:\\Nicola_Gritti\hannah\\20220824_HF_HBMEC-RM-exp2_FF__2022-08-24T15_45_50-Measurement 1'

'''
First, create the csv file for image bookeeping
'''
pymif.pe_opera.xml2csv(exp_folder)
pymif.pe_opera.xml2csv(ff_folder)

'''
Next, create the FF images with ImageJ manually....

Finally, compile experimental images using flat field correction
'''
df = pd.read_csv(os.path.join(exp_folder,'metadata.csv'))

ffs = [
    None,
    imread('DAPI_ff.tif'),
    imread('AF488_ff.tif'),
    imread('AF561_ff.tif'),
    imread('AF647_ff.tif'),
]

pymif.pe_opera.wellplate.compile_conditions_multifields(
    exp_folder, 
    [1,2,3,0,4], 
    ['gray','blue','green','orange','red'], 
    df,
    ffs
)
