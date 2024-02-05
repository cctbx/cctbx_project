from __future__ import absolute_import, division, print_function
import sys
from libtbx.utils import Sorry
import cctbx.eltbx

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def run(filename, log = None):
  '''
    Read custom scattering factors from external file
  '''
  if log is None:
    log = sys.stdout
  print('Opening file with custom scattering factors: ', filename, file=log)
  with open(filename,'r') as f:
     all_lines = f.readlines()

  new_scattering_dictionary = dict()

  for l in all_lines:
    fields = l.split()
    el = fields[0]
    parameters = fields[1:]
    # Make sure the first field is a string
    # assert(isinstance(el, basestring))
    # Python 3
    assert (isinstance(el, str))

    # How many gaussians?
    n_gauss = (len(parameters))//2

#    # double check what cctbx allows
#    if n_gauss < 4 and n_gauss > 5:
#      raise Sorry('Only 4 or 5 gaussians are supported')

    # Make sure values are floats
    for value in parameters:
      if not isfloat(value):
        raise Sorry('Values in %s are not floats.' % filename)

    # Is there a constant?
    has_constant = False
    if (len(parameters))%2 == 1:
      has_constant = True

    # Get arrays of a and b
    def vals(i,j): return [float(s) for s in parameters[i:j]]
    array_of_a = vals(0,n_gauss)
    array_of_b = vals(n_gauss,2*n_gauss)

    if has_constant:
      constant = float(parameters[-1])

    if has_constant:
      new_scattering_dictionary[el] = cctbx.eltbx.xray_scattering.gaussian(
        tuple(array_of_a),
        tuple(array_of_b),
        constant)
    else:
      new_scattering_dictionary[el] = cctbx.eltbx.xray_scattering.gaussian(
        tuple(array_of_a),
        tuple(array_of_b))

  return new_scattering_dictionary

