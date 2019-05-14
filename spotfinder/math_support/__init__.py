from __future__ import absolute_import, division, print_function
def pixels_to_mmPos(x,y,pixel_size):
  return [pixel_size*x,pixel_size*y]

def scitbx_stats(data):
  from scitbx.math import basic_statistics
  bs = basic_statistics(values = data)
  return bs.mean, bs.bias_corrected_standard_deviation

def stats_profile(data):
  from scitbx.array_family import flex
  fdata = flex.double()
  for item in data:
    fdata.append(item)
  perm = flex.sort_permutation(fdata)
  percentile05 = int(0.05*len(fdata))
  percentile95 = int(0.95*len(fdata))

  '''The return value is the ratio of the 95%ile value to the 5%ile value.
  This gives some measure of how braod the spread is between big and small
  spots.'''
  return fdata[perm[percentile95]]/fdata[perm[percentile05]]
