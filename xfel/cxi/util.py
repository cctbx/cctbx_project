from __future__ import division
import os

#for integration pickles:
allowable_basename_endings = ["_00000.pickle",
                              ".pickle",
                              "_refined_experiments.json",
                              "_experiments.json"
                             ]
def is_odd_numbered(file_name):
  for allowable in allowable_basename_endings:
    if (file_name.endswith(allowable)):
      return int(os.path.basename(file_name).split(allowable)[0][-1])%2==1
  raise ValueError
