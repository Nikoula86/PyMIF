import tqdm
from _import_metadata_viventis import import_metadata_viventis
from _open_tif import open_tif
import dask.array as da

def load_viventis_dataset(
                        path: str,
                        zyx_chunks = None,
                        ):
    """_summary_

    Args:
        path (str): _description_
        zyx_chunks (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """

    meta = import_metadata_viventis( path )

    arrays = [ [ None for ch in range(meta["size_c"]) ] for tp in range(meta["size_t"]) ]
    
    for ch in range(meta["size_c"]):
    
        for tp in tqdm.tqdm( range(meta["size_t"]), total=meta["size_t"] ):

            
            file_template = meta["file_template"]
            fname = f"{path}/{file_template}"
            fname = fname.replace("{t}", f"{(tp+1):04}")
            fname = fname.replace( "{c}", meta[f"channel_{ch}_name"] )

            arrays[tp][ch] = open_tif(fname, chunks=zyx_chunks)

    arrays = da.stack(arrays, axis=0)
        
    return arrays, meta
