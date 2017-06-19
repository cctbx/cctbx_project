from __future__ import absolute_import, division
#!/usr/bin/env python
# FormatPilatusHelpers.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Helper methods for class for working with Pilatus images, for instance for
# identifying the regions to be masked.

def pilatus_6M_mask():
  '''Hard coded mask regions for a Pilatus 6M instrument.'''
  # FIX me, the paramters are listed here as f0, f1, s0, s1 but the prototype specifies f0, s0, f1, s1
  return [[488, 494, 1, 2527],
          [982, 988, 1, 2527],
          [1476, 1482, 1, 2527],
          [1970, 1976, 1, 2527],
          [1, 2463, 196, 212],
          [1, 2463, 408, 424],
          [1, 2463, 620, 636],
          [1, 2463, 832, 848],
          [1, 2463, 1044, 1060],
          [1, 2463, 1256, 1272],
          [1, 2463, 1468, 1484],
          [1, 2463, 1680, 1696],
          [1, 2463, 1892, 1908],
          [1, 2463, 2104, 2120],
          [1, 2463, 2316, 2332]]

def pilatus_2M_mask():
  '''Hard coded mask regions for a Pilatus 2M detector.'''

  return [[488, 494, 1, 1679],
          [982, 988, 1, 1679],
          [1, 1475, 196, 212],
          [1, 1475, 408, 424],
          [1, 1475, 620, 636],
          [1, 1475, 832, 848],
          [1, 1475, 1044, 1060],
          [1, 1475, 1256, 1272],
          [1, 1475, 1468, 1484]]

def pilatus_300K_mask():
  '''Hard coded mask regions for a Pilatus 300K instrument.'''

  return [[1, 487, 196, 212],
          [1, 487, 408, 424]]

def determine_pilatus_mask(xdetector):
  '''Return an appropriate pixel mask for a Pilatus detector.'''

  size = xdetector[0].get_image_size()

  # Hardcoded module size and gap size
  module_size_fast, module_size_slow = (487, 195)
  gap_size_fast, gap_size_slow = (7, 17)

  # Edge dead areas not included, only gaps between modules matter
  n_fast, remainder = divmod(size[0], module_size_fast)
  assert (n_fast-1) * gap_size_fast == remainder

  n_slow, remainder = divmod(size[1], module_size_slow)
  assert (n_slow-1) * gap_size_slow == remainder

  # Specify the dead areas between the modules, i.e. the rows and columns
  # where there are no active pixels
  mask = []
  for i_fast in range(n_fast-1):
    mask.append([
      (i_fast+1) * module_size_fast + i_fast * gap_size_fast + 1,
      (i_fast+1) * module_size_fast + i_fast * gap_size_fast + gap_size_fast,
      1, size[1]
    ])
  for i_slow in range(n_slow-1):
    mask.append([
      1, size[0],
      (i_slow+1) * module_size_slow + i_slow * gap_size_slow + 1,
      (i_slow+1) * module_size_slow + i_slow * gap_size_slow + gap_size_slow
    ])

  return mask

def get_vendortype(xdetector):
  array = xdetector[0].get_image_size()
  if array == (2463,2527): return "Pilatus-6M"
  elif array == (2463,195): return "Pilatus-6M" # special treatment of Pilatus-12M, treat as -6M for the viewer
  elif array == (1475,1679): return "Pilatus-2M"
  elif array == (487,619): return "Pilatus-300K"
  return "Undetermined Pilatus size"

def get_vendortype_eiger(xdetector):
  array = xdetector[0].get_image_size()
  #print array,
  if array == (4150,4371): return "Eiger-16M"
  elif array == (3110,3269): return "Eiger-9M"
  elif array == (2070,2167): return "Eiger-4M"
  elif array == (1030,1065): return "Eiger-1M"
  return "Undetermined EigerX size"
