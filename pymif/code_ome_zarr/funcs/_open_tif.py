from skimage.io import imread
import dask.array as da
from typing import Any, List, Optional, Sequence, Tuple, Union

def open_tif(
            tif_file: str,
            chunks: Optional[Sequence[int]] = None,
            ):
    """_summary_

    Args:
        tif_file (str): _description_
        chunks (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """

    image = imread( tif_file, aszarr=True )
    image = da.from_zarr( image, chunks=chunks )
    
    return image
