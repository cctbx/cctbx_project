from __future__ import absolute_import, division, print_function

from cctbx.development import random_structure
from cctbx.sgtbx import space_group_info
from six.moves import range
from cctbx_asymmetric_map_ext import *
from cctbx.array_family import flex
from cctbx import maptbx
import random
from cctbx.maptbx import is_periodic

if (1):
  kk=0
  random.seed(kk)
  flex.set_random_seed(kk)

def run(single_side=False, quick=True):
  result = flex.bool()
  periodic_list=[]
  aperiodic_list=[]
  if quick:
    rf_list=[0.33]
  else:
    rf_list=[0.33, 0.48, ]
  for resolution_factor in rf_list:
   for sgn in range(1,231):
    kk=random.randint(1,99999)
    random.seed(kk)
    flex.set_random_seed(kk)
    group = space_group_info(sgn)
    xrs = random_structure.xray_structure(
      space_group_info       = group,
      volume_per_atom        = 25.,
      general_positions_only = False,
      elements               = ('C', 'N', 'O', 'H')*30,
      min_distance           = 1.0)
    sgt = xrs.space_group().type()
    fc = xrs.structure_factors(d_min=2).f_calc()
    fft_map = fc.fft_map(symmetry_flags = maptbx.use_space_group_symmetry,
      resolution_factor=resolution_factor)
    map_data = fft_map.real_map_unpadded()
    xx=is_periodic(map_data)
    periodic_list.append(xx)
    assert xx in [True,None]

    # Now cut off edges and should not work:
    min_size=int(flex.double(map_data.all()).min_max_mean().min)
    assert min_size >= 12
    if quick:
      range_to_use=[1]
      k_list=[1]
    else:
      range_to_use=list(range(1,max(2,min(10,1+min_size//12))))
      k_list=[0,1]
    for i in range_to_use:
     for k in k_list:
      if single_side:
        lower_bounds=(k,k,k)
        upper_bounds=tuple([x - i - 1 for x in map_data.all()])
      else:
        lower_bounds=(k,k,k)
        upper_bounds=(map_data.all()[0]-i-1,map_data.all()[1],map_data.all()[2])

      map_box = maptbx.copy(map_data,
       lower_bounds,upper_bounds)
      new_map_box=map_box.shift_origin()
      xx=is_periodic(new_map_box)
      aperiodic_list.append(xx)
      assert xx in [False,None]


  print ("For periodic.. True:",periodic_list.count(True),
     "None:",periodic_list.count(None)," False:",periodic_list.count(False))
  print ("For aperiodic.. True:",aperiodic_list.count(True),
     "None:",aperiodic_list.count(None)," False:",aperiodic_list.count(False))
  assert periodic_list.count(False) == 0
  assert periodic_list.count(None) <= 0.05 * len(periodic_list)
  assert aperiodic_list.count(True) == 0
  assert aperiodic_list.count(None) <= 0.05 * len(aperiodic_list)

if (__name__ == "__main__"):
  run()
  print ("OK")
