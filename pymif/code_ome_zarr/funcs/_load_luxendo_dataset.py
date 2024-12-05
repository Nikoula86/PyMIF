import tqdm
from _import_metadata_luxendo import import_metadata_luxendo
from _open_h5 import open_h5
from typing import Any, List, Optional, Sequence, Tuple, Union
import dask.array as da

def load_luxendo_dataset(
                        path: str,
                        zyx_chunks: Optional[Sequence[int]] = None,
                        ):
    """_summary_

    Args:
        path (str): _description_
        zyx_chunks (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """

    meta = import_metadata_luxendo( path )

    arrays = [ [ None for ch in range(meta["size_c"]) ] for tp in range(meta["size_t"]) ]
    
    for ch in range(meta["size_c"]):
    
        for tp in tqdm.tqdm( range(meta["size_t"]), total=meta["size_t"] ):

            
            file_template = meta["file_template"]
            fname = f"{path}/{file_template}"
            fname = fname.replace("{t}", f"{tp}")
            fname = fname.replace( "{c}", meta[f"channel_{ch}_ID"] )

            arrays[tp][ch] = open_h5(fname, chunks=zyx_chunks)

    arrays = da.stack(arrays, axis=0)
        
    return arrays, meta
