import sys
sys.path.append("../funcs")
from _open_h5 import open_h5

h5_file = "../../test_data/luxendo/2024-10-30_141806_Exp3_S40/processed/task_C/uni_tp-0_ch-0_st-0_obj-right_cam-right_etc.lux.h5"

image = open_h5(h5_file)

assert image.shape == (2895,2076,2860)
assert image[1300,1000,1000].compute() == 33
