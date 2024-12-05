import sys
sys.path.append("../funcs")
from _open_tif import open_tif

tif_file = "../../test_data/viventis/20241104_162954_Experiment/Position 1_Settings 1/t0001_FM4-64.tif"

image = open_tif(tif_file)

assert image.shape == (81,2304,2304)
assert image[10,10,10].compute() == 99
