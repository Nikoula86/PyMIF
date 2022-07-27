import numpy as np
import pyclesperanto_prototype as cle
from skimage.util import img_as_float32
import pandas as pd
import math, tqdm
from skimage.io import imread
import matplotlib.pyplot as plt

def normalize(
            img, 
            percentiles=[3,99.99], 
            percentiles_downsample=[2,2,2]
            ):

    percentiles = np.array(percentiles)
    percentiles_downsample = np.array(percentiles_downsample)
    
    assert (percentiles.ndim==1)&(len(percentiles)==2), 'Detected incompatible percentiles definition!'
    assert (percentiles_downsample.ndim==1)&(len(percentiles_downsample)==3), 'Detected incompatible percentiles_downsample definition!'
    
    
    # to make it more robust, normalise values before finding peaks.
    print('To float...')
    img_float = img_as_float32(img)

    # find 3 and 97 percentiles (typically used by image analysis, otherwise try 10 and 90 percent)
    print('Percs...')
    img_down = img_float[::percentiles_downsample[0], ::percentiles_downsample[1], ::percentiles_downsample[2]]
    percs_vals = np.percentile(img_down, [percentiles[0], percentiles[1]])

    print('Normalize...')
    img_float = (img_float - percs_vals[0]) / (percs_vals[1] - percs_vals[0])

    print('Clip...')
    img_float = np.clip(img_float, 0, 1)
    
    print('Done.')
    return img_float

def runDoG(
        img, 
        cell_diameter=[5.,20.,20.], 
        top_hat_radius=[5.,20.,20.],
        gpu_downsample=[1,2,2], 
        ):

    gpu_downsample = np.array(gpu_downsample)
    cell_diameter = np.array(cell_diameter)
    top_hat_radius = np.array(top_hat_radius)
    
    assert (gpu_downsample.ndim==1)&(len(gpu_downsample)==3), 'Detected incompatible gpu_downsample definition!'
    assert (cell_diameter.ndim==1)&(len(cell_diameter)==3), 'Detected incompatible cell_diameter definition!'
    assert (top_hat_radius.ndim==1)&(len(top_hat_radius)==3), 'Detected incompatible top_hat_rad definition!'

    diameter_gpu = cell_diameter/gpu_downsample
    top_hat_radius_gpu = top_hat_radius/gpu_downsample
    
    print('Push to GPU...')
    img_gpu = cle.push(img[::gpu_downsample[0], ::gpu_downsample[1], ::gpu_downsample[2]])

    # correct for scattering and increased autofluorescence in some regions of the sample
    print('TopHat...')
    img_gpu = cle.top_hat_box(img_gpu, 
                            radius_x=top_hat_radius_gpu[2], 
                            radius_y=top_hat_radius_gpu[1], 
                            radius_z=top_hat_radius_gpu[0])

    d_gpu = cell_diameter/2.

    print('Sigma1...')
    d1 = diameter_gpu/(1.+np.sqrt(2))
    blurred1 = cle.gaussian_blur(img_gpu, sigma_x=d1[2], sigma_y=d1[1], sigma_z=d1[0])

    print('Sigma2...')
    d2 = d1*np.sqrt(2)
    blurred2 = cle.gaussian_blur(img_gpu, sigma_x=d2[2], sigma_y=d2[1], sigma_z=d2[0])

    print('DoG...')
    img_DoG = blurred1-blurred2
    
    print('Done.')
    return img_DoG

def runLocalMax(
            img_DoG, 
            thr_DoG=0.01, 
            gpu_downsample=[1,2,2],
            lims=np.array([[50,50],[100,100],[100,100]])
            ):

    gpu_downsample = np.array(gpu_downsample)

    assert (gpu_downsample.ndim==1)&(len(gpu_downsample)==3), 'Detected incompatible gpu_downsample definition!'

    print('Detect local maxima...')
    detected_spots = cle.detect_maxima_box(img_DoG, radius_x=0, radius_y=0, radius_z=0)
    selected_spots = cle.binary_and(img_DoG>thr_DoG, detected_spots)
    
    print('Detect coordinates...')
    p = np.where(selected_spots)

    df = pd.DataFrame({
                        'z': p[0].astype(float)*gpu_downsample[0],
                        'y': p[1].astype(float)*gpu_downsample[1],
                        'x': p[2].astype(float)*gpu_downsample[2],
                        })

    print('Filter coordinates...')
    df = df[(df.z>lims[0,0])&(df.z<(img_DoG.shape[0]*gpu_downsample[0]-lims[0,1]))&
            (df.y>lims[1,0])&(df.y<(img_DoG.shape[1]*gpu_downsample[1]-lims[1,1]))&
            (df.x>lims[2,0])&(df.x<(img_DoG.shape[2]*gpu_downsample[2]-lims[2,1]))]

    print('Done.')
    return df

def unit_sphere(df, pixel_size):
    """
    
    Parameters
    ----------
    df : dataframe
        DESCRIPTION.
    pixel_size : list, int
        DESCRIPTION.
    save : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    df : dataframe
        dataframe with positions normalized on the unit sphere.

    """

    def sphereFit(spX,spY,spZ):
        #   Assemble the A matrix
        spX = np.array(spX)
        spY = np.array(spY)
        spZ = np.array(spZ)
        A = np.zeros((len(spX),4))
        A[:,0] = spX*2
        A[:,1] = spY*2
        A[:,2] = spZ*2
        A[:,3] = 1

        #   Assemble the f matrix
        f = np.zeros((len(spX),1))
        f[:,0] = (spX*spX) + (spY*spY) + (spZ*spZ)
        C, residules, rank, singval = np.linalg.lstsq(A,f)

        #   solve for the radius
        t = (C[0]*C[0])+(C[1]*C[1])+(C[2]*C[2])+C[3]
        radius = math.sqrt(t)

        return radius, C[0], C[1], C[2]

    df1 = df.copy()

    # normalize position on unit sphere
    df1['z_unit'] = df1.z*pixel_size[0]
    df1['y_unit'] = df1.y*pixel_size[1]
    df1['x_unit'] = df1.x*pixel_size[2]
    
    r, c0, c1, c2 = sphereFit(df1.z_unit, df1.y_unit, df1.x_unit)

    df1.z_unit = (df1.z_unit-c0)/r
    df1.y_unit = (df1.y_unit-c1)/r
    df1.x_unit = (df1.x_unit-c2)/r

    return df1

def compute_loc_fluo(df, imgs, cell_diameter=[2.,10.,10.], ch_names=None):
    """
    

    Parameters
    ----------
    df : dataframe
        pandas dataframe containing location info.
    imgs : list, str
        list of images.
    ch_names : list, optional
        name of channels. Must match number and order of files given. The default is None.
        
    Returns
    -------
    df : dataframe
        return dataframe populated with fluorescence intensity values.

    """
    ### load fluo images
    # imgs = np.stack([imread(f) for f in tqdm.tqdm(files)])
    
    n_cells = len(df)
    n_ch = len(imgs)

    df1 = df.copy()
    
    if ch_names == None:
        ch_names=['ch_'+str(i) for i in range(len(n_ch))]
    else:
        ch_names = [str(i) for i in ch_names]
        
    assert len(ch_names)==n_ch

    # initialize with empty fluorescence intensities
    for ch_name in ch_names:
        df1[ch_name] = 0.

    ### compute fluo intensity
    for i in tqdm.tqdm(range(n_cells), total=n_cells):
        cell = df1.loc[i]
        x,y,z = int(cell.x), int(cell.y), int(cell.z)
        r_z = int(np.round(cell_diameter[0]/2.))
        r_y = int(np.round(cell_diameter[1]/2.))
        r_x = int(np.round(cell_diameter[2]/2.))
        
        fluos = [np.mean(img[z-r_z:z+r_z,
                             y-r_y:y+r_y,
                             x-r_x:x+r_x].astype(float)) for img in imgs]
        
        for ch_name, fluo in zip(ch_names, fluos):
            df1.loc[i,ch_name] = fluo
            
    return df1

def rotate_animation_df(df, col_names, ncols=4, elev=15, azim_init=-60, 
            figsize=(15,15), save=False, interval=10, folder='', show_ticks=False):
    """
    

    Parameters
    ----------
    df : dataframe
        DESCRIPTION.
    col_names : list, str
        should contain a valid column value.
    ncols : int
    elev : float
    azim : float
    figsize : tuple
    save : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    
    # check that all keys provided are valid
    assert len(col_names)>0
    assert all([i in df.keys() for i in col_names])
    
    from matplotlib.animation import FuncAnimation
    fig = plt.figure(figsize=figsize)
    axs = []
    nrows = (len(col_names)-1)//ncols + 1
    
    i = 0
    for col_name in col_names:
        ax = fig.add_subplot(nrows,ncols,i+1,projection='3d')
        ax.scatter(
            df.x_unit, df.y_unit, df.z_unit, 
            c = df[col_name], 
            cmap = 'plasma',
            alpha=0.2
           )
            
        ax.title.set_text(col_name)
        if not show_ticks:
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
            ax.zaxis.set_ticklabels([])
        axs.append(ax)
        
        i+=1
        
    plt.tight_layout()
        
    def update(frame):
        for ax in axs:
            ax.view_init(elev=elev, azim=azim_init+frame)
    
    ani = FuncAnimation(fig, update, frames=np.linspace(0, 365, 36), interval=interval)
    
    return ani
