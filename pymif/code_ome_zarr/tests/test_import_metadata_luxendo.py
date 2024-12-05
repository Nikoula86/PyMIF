import sys
sys.path.append("../funcs")
from _import_metadata_luxendo import import_metadata_luxendo

path = "../../test_data/luxendo/2024-10-30_141806_Exp3_S40/processed/task_C"

meta = import_metadata_luxendo(path)

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
