import os, sys, shutil
sys.path.append("../funcs")
from _load_luxendo_dataset import load_luxendo_dataset
from _load_viventis_dataset import load_viventis_dataset
from _save2omezarr import save2omezarr
from _ome_zarr import OmeZarr

### test luxendo

path = "../../test_data/luxendo/2024-10-30_141806_Exp3_S40/processed/task_C"
output_zarr = "../../test_data/zarr/luxendo/2024-10-30_141806_Exp3_S40"

dataset, meta = load_luxendo_dataset(
                            path,
                            zyx_chunks = (256,256,256),
                            )

if os.path.exists(output_zarr):
    shutil.rmtree(output_zarr)

save2omezarr(dataset, 
             meta, 
             output_zarr,
             n_scales = 3,
             zyx_chunks = (256,256,256),
             )

ds_raw = OmeZarr(output_zarr, mode="r")
meta = ds_raw.get_meta()
image = ds_raw.read_stack(time_point=1, channel=0, res_level=0)
ds_raw.close()

print(image)

