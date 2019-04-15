from __future__ import division
from libtbx import adopt_init_args
from cctbx import maptbx

class run(object):
  def __init__(self, model, map_data):
    adopt_init_args(self, locals())
    # Find blob
    co = maptbx.connectivity(map_data = self.map_data, threshold = 5.)
    connectivity_map = co.result()
    sorted_by_volume = sorted(
      zip(co.regions(), range(0, co.regions().size())), key=lambda x: x[0],
        reverse=True)
    blob_indices = []
    for p in sorted_by_volume:
      v, i = p
      print v, i
      if(i>0):
        blob_indices.append(i)
