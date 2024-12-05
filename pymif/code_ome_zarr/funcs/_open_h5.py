import h5py
import dask.array as da
from typing import Any, List, Optional, Sequence, Tuple, Union

def open_h5(
            h5_file: str,
            chunks: Optional[Sequence[int]] = None,
            ):
    """_summary_

    Args:
        h5_file (str): _description_

    Returns:
        _type_: _description_
    """

    dataset = h5py.File(h5_file, "r")["Data"]
    if chunks is None:
        chunks = list( dataset.shape )
        chunks[0] = 1
    image = da.from_array( dataset, chunks=chunks )
    
    return image
