import sys
sys.path.append("../funcs")
from _load_luxendo_dataset import load_luxendo_dataset

path = "../../test_data/luxendo/2024-10-30_141806_Exp3_S40/processed/task_C"

dataset, meta = load_luxendo_dataset(
                            path,
                            zyx_chunks = None,#(256,256,256),
                            )

meta_out = {
    "size_t": 3,
    "size_c": 2,
    "size_z": 2895,
    "size_y": 2076,
    "size_x": 2860,
    "dt": 1.,
    "dz": 0.39000002,
    "dy": 0.39000002,
    "dx": 0.39000002,
    "unit_t": "s",
    "unit_z": "micrometer",
    "unit_y": "micrometer",
    "unit_x": "micrometer",
    "channel_0_ID": "0",
    "channel_0_name": "ch:0",
    "channel_0_color": "0",
    "channel_1_ID": "1",
    "channel_1_name": "ch:1",
    "channel_1_color": "1",
    "file_template": "uni_tp-{t}_ch-{c}_st-0_obj-right_cam-right_etc.lux.h5"
}

for key in meta:
    assert meta[key] == meta_out[key], f"{key} is different! \
                Expected {meta[key]}, obtained {meta_out[key]} instead."
    
expected_shape = (3,2,2895,2076,2860)
assert dataset.shape == expected_shape, f"Output shape is different! \
                Expected {expected_shape}, obtained {dataset.shape} instead."

pxl_val = dataset[2,1,1400,1000,1000].compute()
expected_pxl_val = 565
assert pxl_val == expected_pxl_val, f"Pixel value different! \
                Expected {expected_pxl_val}, obtained {pxl_val} instead."

