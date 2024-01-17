from __future__ import absolute_import, division, print_function

from cctbx.maptbx.qscore import qscore_np as qscore
from iotbx.data_manager import DataManager


dm = DataManager()
dm.process_real_map_file(
    "/Users/user/Desktop/data/data_large/7UAE/emd_26422.map.gz")
dm.process_model_file(
    "/Users/user/Desktop/data/data_large/7UAE/pdb7uae.ent.gz")
mmm = dm.get_map_model_manager()
result = qscore(mmm)
print(result)
