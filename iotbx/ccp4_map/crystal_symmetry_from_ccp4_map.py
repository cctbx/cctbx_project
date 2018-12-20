from __future__ import division
import iotbx.ccp4_map
from cctbx import crystal
import sys, os

SHIFT_TO_MRCFILE=False

def extract_from(file_name):

  if SHIFT_TO_MRCFILE:
    from iotbx.mrcfile.crystal_symmetry_from_ccp4_map \
        import extract_from as extract_using_mrcfile
    return extract_using_mrcfile(file_name)

  # XXX This is to hide stdout from ccp4io code (ccp4_cmap_open) about input
  # XXX being non-ccp4 map. That's why we never have io in c/c++ code!
  #
  # redirect sys.stdout
  save_stdout = sys.stdout
  oldstdout_fno = os.dup(sys.stdout.fileno())
  sys.stdout.flush()
  newstdout = os.dup(1)
  devnull = os.open(os.devnull, os.O_WRONLY)
  os.dup2(devnull, 1)
  os.close(devnull)
  sys.stdout = os.fdopen(newstdout, 'w')
  ####
  try: # run application
    m = iotbx.ccp4_map.map_reader(file_name=file_name)
  except: # intentional
    # restore sys.stdout
    sys.stdout = save_stdout
    sys.stdout.flush()
    os.dup2(oldstdout_fno, 1)
  # restore sys.stdout
  sys.stdout = save_stdout
  sys.stdout.flush()
  os.dup2(oldstdout_fno, 1)
  # continue as normal
  ucp = m.unit_cell_parameters
  sgn = max(1, m.space_group_number)
  return crystal.symmetry(ucp, sgn)
