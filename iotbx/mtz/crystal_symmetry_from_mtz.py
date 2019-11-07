"Extracts crystal symmetry from MTZ file."
from __future__ import absolute_import, division, print_function

from iotbx import mtz
from cctbx import crystal

def extract_from(file_name):
  mtz_object = mtz.object(file_name=file_name)
  assert mtz_object.n_symmetry_matrices() > 0
  return mtz_object.crystals()[0].crystal_symmetry()
