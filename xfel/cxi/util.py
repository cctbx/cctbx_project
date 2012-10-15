from __future__ import division
import os
def is_odd_numbered(file_name):
      if (file_name.endswith("_00000.pickle")):
        return int(os.path.basename(file_name).split("_00000.pickle")[0][-1])%2==1
      elif (file_name.endswith(".pickle")):
        return int(os.path.basename(file_name).split(".pickle")[0][-1])%2==1

