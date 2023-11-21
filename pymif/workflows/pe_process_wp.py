import sys, os
import pandas as pd
sys.path.append('Z:\\people\\gritti\\code\\pymif')
from pymif import pe_opera

exp_folder = 'Y:\\Krisztina Arato\\results\\PE\\20221110\\KrA-20221110-alldoxy__2022-11-10T09_39_30-Measurement 1'
# ff_folder = 'Y:\\Nicola_Gritti\hannah\\20220824_HF_HBMEC-RM-exp2_FF__2022-08-24T15_45_50-Measurement 1'

'''
First, create the csv file for image bookeeping
'''
# pe_opera.xml2csv(exp_folder)

'''
Finally, compile experimental images using flat field correction
'''
df = pd.read_csv(os.path.join(exp_folder,'metadata.csv'))

conditions = [['alldoxy' for i in range(8)] for j in range(12)]

pe_opera.wellplate.compile_conditions(
        path=exp_folder, 
        conditions=conditions, 
        channel_order=[1,2,0], 
        luts_name=["gray","green","orange"],
        df=df,
        ffs=None,
        downsample=1.,
        ff_mode = 'PE', 
        outfolder = 'compiled',
        image_folder = os.path.join("Images"),
        which_proj = 'none',
)

