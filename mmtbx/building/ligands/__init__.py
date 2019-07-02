from __future__ import absolute_import, division, print_function
from libtbx import adopt_init_args
from cctbx import maptbx
from scitbx.array_family import flex

class run(object):
  def __init__(self, model, map_data):
    adopt_init_args(self, locals())
    # Find blob
    co = maptbx.connectivity(map_data = self.map_data, threshold = 5.)
    #connectivity_map = co.result()
    #sorted_by_volume = sorted(
    #  zip(co.regions(), range(0, co.regions().size())), key=lambda x: x[0],
    #    reverse=True)
    #blob_indices = []
    #for p in sorted_by_volume:
    #  v, i = p
    #  print v, i
    #  if(i>0):
    #    blob_indices.append(i)
    #######
    # You get everything you need:
    map_result = co.result()
    volumes = co.regions()
    print(volumes)
    coors = co.maximum_coors()
    vals = co.maximum_values()
    minb, maxb = co.get_blobs_boundaries_tuples()
    # This will give you the order
    i_sorted_by_volume = flex.sort_permutation(data=volumes, reverse=True) # maybe co.regions() should go there
    for i in i_sorted_by_volume:
      print("blob #", i)
      print(coors[i])
      print(vals[i])
      print(maxb[i], minb[i])
