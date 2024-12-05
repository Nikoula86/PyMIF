import tqdm
from typing import Any, List, Optional, Sequence, Tuple, Union
from _ome_zarr import OmeZarr

def save2omezarr(
    array,
    meta: dict,
    output_zarr: str,
    n_scales: Optional[int] = 3,
    zyx_chunks: Optional[Sequence[int]] = None,
    ):

    if zyx_chunks is None:
        zyx_chunks = list( array.shape[2:] )
        zyx_chunks[0] = 1

    ds_raw = OmeZarr(output_zarr, mode="a")
    ds_raw.create_dataset(meta, n_scales=n_scales, zyx_chunks=zyx_chunks)

    for ch in range(meta["size_c"]):
        ch_name = meta[f"channel_{ch}_name"]
        print(f"\t Converting channel \"{ch_name}\"")
        
        for tp in tqdm.tqdm(range(meta["size_t"]), total=meta["size_t"]):
            
            ds_raw.write_stack(array[tp, ch, :], tp, ch)

    ds_raw.close()
    
    return
