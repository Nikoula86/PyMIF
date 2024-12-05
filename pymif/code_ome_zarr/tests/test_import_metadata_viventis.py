import sys
sys.path.append("../funcs")
from _import_metadata_viventis import import_metadata_viventis

path = "../../test_data/viventis/20241104_162954_Experiment/Position 1_Settings 1"

meta = import_metadata_viventis(path)

meta_out = {
    "size_t": 5,
    "size_c": 2,
    "size_z": 81,
    "size_y": 2304,
    "size_x": 2304,
    "dt": 900.,
    "dz": 1.,
    "dy": 0.173,
    "dx": 0.173,
    "unit_t": "s",
    "unit_z": "um",
    "unit_y": "um",
    "unit_x": "um",
    "channel_0_ID": "Channel:1",
    "channel_0_name": "Hoechst",
    "channel_0_color": "16711935",
    "channel_1_ID": "Channel:2",
    "channel_1_name": "FM4-64",
    "channel_1_color": "-16776961",
    "file_template": "t{t}_{c}.tif"
}

for key in meta:
    assert meta[key] == meta_out[key], f"{key} is different! \
                Expected {meta[key]}, obtained {meta_out[key]} instead."
