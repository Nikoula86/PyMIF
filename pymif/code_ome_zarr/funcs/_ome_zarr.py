import glob, os, sys, math
from typing import Any, List, Optional, Sequence, Tuple, Union
import numpy as np
from tqdm import tqdm
import pandas as pd
from ome_zarr.io import parse_url
import zarr
import dask.array as da
try:
    from cucim import __version__ as sk_version
    from cucim.skimage.transform import downscale_local_mean

    DOWNSCALE_METHOD = "cucim.skimage.transform.downscale_local_mean"

except ImportError:
    from skimage import __version__ as sk_version
    from skimage.transform import downscale_local_mean

    DOWNSCALE_METHOD = "skimage.transform.downscale_local_mean"

class OmeZarr():
    """_summary_
    """
    
    def __init__(self, 
                 zarrpath: str,
                 mode: str = "a"
                 ):
        """Instantiate a zarr dataset.

        Args:
            zarrpath (str): Path to ome.zarr file
            mode (str, optional): Access mode. Defaults to "a" (read/write).
        """
        
        self.zarrpath = zarrpath
        self.store = parse_url(zarrpath, mode=mode).store
        self.root = zarr.group(store=self.store)
        
    def create_dataset(self,
                meta: dict,
                n_scales: int = 1,
                zyx_chunks: Optional[Sequence[int]] = None,
                dtype = np.uint16,
                ):
        """Generate a group with datasets and metadata. Requires mode="a" or "w".

        Args:
            meta (dict): Dictionary with relevant metadata.
            n_scales (int, optional): Number of resolution layers. Defaults to 1.
            zyx_chunks (Optional[Sequence[int]], optional): 3D chunk size (Z-Y-X). Defaults to None.

        """

        zyx_shape = (int(meta["size_z"]),int(meta["size_y"]),int(meta["size_x"]))
        if zyx_chunks is None:
            zyx_chunks = zyx_shape
            
        datasets = []
        try:
            val_max = np.iinfo(dtype).max
        except:
            val_max = np.finfo(dtype).max

        for i in range(n_scales):
            factor = 2**i
            array_path = f"{i}"
            shape = (meta["size_t"],meta["size_c"],) + tuple(int(math.ceil(s / factor)) for s in zyx_shape)
            chunks = (1,1,) + tuple([math.ceil(c/factor) for c in zyx_chunks])
            self.root.full(array_path, shape=shape, dtype=dtype, chunks=chunks, 
                                            compressor = zarr.Blosc(cname="zstd", clevel=5, 
                                                                    shuffle=zarr.Blosc.BITSHUFFLE),
                                            fill_value = val_max)
            datasets.append(
                {
                    "path": array_path,
                    "coordinateTransformations": [
                        {
                            "type": "scale",
                            "scale": [ float(meta["dt"]), 1., 
                                      float(meta["dz"])*factor, 
                                      float(meta["dy"])*factor, 
                                      float(meta["dx"])*factor ],
                        }
                    ],
                }
            )

        self.root.attrs["multiscales"] = [
                {
                    "datasets": datasets,
                    "axes": [
                        {"name": "t", "type": "time", "unit": meta["unit_t"]},
                        {"name": "c", "type": "channel"},
                        {"name": "z", "type": "space", "unit": meta["unit_z"]},
                        {"name": "y", "type": "space", "unit": meta["unit_y"]},
                        {"name": "x", "type": "space", "unit": meta["unit_x"]},
                    ],
                    "metadata": {
                        "method": DOWNSCALE_METHOD,
                        "version": sk_version,
                    },
                }
            ] 
        
        def default_omero_metadata(name, meta):
            channels = range(int(meta["size_c"]))
            try:
                val_max = np.iinfo(dtype).max
            except:
                val_max = np.finfo(dtype).max
            
            return {
                "id": 1,
                "name": name,
                "channels": [
                    {
                        "active": True,
                        "coefficient": 1,
                        "color": meta[f"channel_{ch}_color"],
                        "family": "linear",
                        "inverted": False,
                        "label": meta[f"channel_{ch}_name"],
                        "window": {"start": 0, "end": 1500, "min": 0, "max": val_max},  # just a guess
                    }
                    for ch in channels
                ],
            }          
        
        self.root.attrs["omero"] = default_omero_metadata(self.zarrpath, 
                                                          meta, 
                                                          )

    def get_n_scales(self):
        
        return len(sorted(self.root.array_keys()))
    
    def get_axes_info(self,
                  res_level: Optional[int]=0):
        
        meta = self.root.attrs.asdict()
        ax_name = tuple( [ a["name"] for a in meta["multiscales"][0]["axes"] ] )
        ax_unit = tuple( [ a["unit"] if "unit" in a.keys() else "" for a in meta["multiscales"][0]["axes"] ] )
        ax_scale = meta["multiscales"][0]["datasets"][res_level]["coordinateTransformations"][0]["scale"]
        ax_shape = self.root[res_level].shape
        
        return ax_name, ax_unit, ax_scale, ax_shape
    
    def get_channels_info(self):
        
        channels_meta = self.root.attrs.asdict()["omero"]["channels"]
        
        ch_idx = [i for i in range(len(channels_meta))]
        ch_name = [channels_meta[i]["label"] for i in range(len(channels_meta))]
        ch_color = [channels_meta[i]["color"] for i in range(len(channels_meta))]
        
        return ch_idx, ch_name, ch_color

    def get_meta(self):
        
        meta_out = {}
        meta_out["axes_order"] = ""

        axes_name, axes_unit, axes_scale, axes_shape = self.get_axes_info()
        for i, ax_name in enumerate(axes_name):
            meta_out[f"size_{ax_name}"] = axes_shape[i]
            meta_out[f"d{ax_name}"] = axes_scale[i]
            meta_out[f"unit_{ax_name}"] = axes_unit[i]
            meta_out[f"chunk_size_{ax_name}"] = self.root[0].chunks[i]
            meta_out["axes_order"] += ax_name
            
        chs_idx, chs_name, chs_color = self.get_channels_info()
        for ch_idx in chs_idx:
            meta_out[f"channel_{ch_idx}_ID"] = ch_idx
            meta_out[f"channel_{ch_idx}_name"] = chs_name[ch_idx]
            meta_out[f"channel_{ch_idx}_color"] = chs_color[ch_idx]
            
        meta_out["n_scales"] = self.get_n_scales()
            
        return meta_out
        
    def write_stack(self,
                    stack,
                    timepoint: int,
                    channel: int):
        """Write a stack into the dataset.

        Args:
            stack (np.ndarray or da.Array): Image stack (3D, ZYX).
            timepoint (int): Time index.
            channel (int): Channel index.
        """
        
        ###
        # https://forum.image.sc/t/writing-tile-wise-ome-zarr-with-pyramid-size/85063/2
        # https://zarr.readthedocs.io/en/stable/api/hierarchy.html#zarr.hierarchy.Group.require_dataset
        ###
        n_scales = self.get_n_scales()

        zyx_shapes = [ self.root[i].shape[2:] for i in range(n_scales)]
        zyx_chunk_sizes = [ self.root[i].chunks[2:] for i in range(n_scales)]

        zarrays = [ self.root.require_dataset(
            str(i),
            shape = self.root[i].shape,
            exact = True,
            chunks = self.root[i].chunks,
            dtype = self.root[i].dtype,
        ) for i in range(n_scales) ]
        
        # print([type(zarray), zarray.shape) for zarray in zarrays])
        
        tile_idx = [math.ceil(s/c) for s,c in zip( zyx_shapes[0], zyx_chunk_sizes[0] )]
        
        for z_idx in range(tile_idx[0]):
            
            for y_idx in range(tile_idx[1]):
                
                for x_idx in range(tile_idx[2]):
                    
                    # print("z:%d/%d - y:%d/%d - x:%d/%d"%(z_idx,tile_idx[0],y_idx,tile_idx[1],x_idx,tile_idx[2]))
                    z1 = z_idx * zyx_chunk_sizes[0][0]
                    z2 = z1 + zyx_chunk_sizes[0][0]
                    y1 = y_idx * zyx_chunk_sizes[0][1]
                    y2 = y1 + zyx_chunk_sizes[0][1]
                    x1 = x_idx * zyx_chunk_sizes[0][2]
                    x2 = x1 + zyx_chunk_sizes[0][2]
                    
                    if isinstance(stack, np.ndarray):
                        tile = stack[z1:z2, y1:y2, x1:x2]
                    elif isinstance(stack, da.Array):
                        tile = stack[z1:z2, y1:y2, x1:x2].compute()
                        
                    zarrays[0][timepoint, channel, z1:z2, y1:y2, x1:x2] = tile
            
                    for i in range(1, n_scales):
                        z1 = z_idx * zyx_chunk_sizes[i][0]
                        z2 = z1 + zyx_chunk_sizes[i][0]
                        y1 = y_idx * zyx_chunk_sizes[i][1]
                        y2 = y1 + zyx_chunk_sizes[i][1]
                        x1 = x_idx * zyx_chunk_sizes[i][2]
                        x2 = x1 + zyx_chunk_sizes[i][2]

                        # print("Resolution level:", i)
                        factors = (2**i,) * stack.ndim
                        tile_down = downscale_local_mean(tile, factors)
                        # print(zarrays[1][timepoint, channel, z1:z2, y1:y2, x1:x2].shape)
                        # print(tile_down.shape)
                        zarrays[i][timepoint, channel, z1:z2, y1:y2, x1:x2] = tile_down
                    # print(type(stack), stack.shape)
                    # print(type(zarray), zarray.shape)
                    # print(type(tile), tile.shape)
                    # print(type(zarray[timepoint, channel, z1:z2, y1:y2, x1:x2]), zarray[timepoint, channel, z1:z2, y1:y2, x1:x2].shape)
            
    def read_stack(self, 
                   time_point: int, 
                   channel: int, 
                   res_level: Optional[int] = 0):
        """Read a single timepoint/channel as a numpy ndarray.

        Args:
            time_point (int): Index of the timepoint to open.
            channel (int): Index of the channel to open.
            res_level (Optional[int], optional): Resolution layer to open. Defaults to 0.

        Returns:
            ndarray: Stack as numpy array.
        """
        
        zarray = self.root.require_dataset(
                    str(res_level),
                    shape = self.root[res_level].shape,
                    exact = True,
                    chunks = self.root[res_level].chunks,
                    dtype = self.root[res_level].dtype,
                )
        
        ndarray = zarray[time_point, channel, ...]
           
        return ndarray
                
    def close(self):
        
        self.store.close()       
                