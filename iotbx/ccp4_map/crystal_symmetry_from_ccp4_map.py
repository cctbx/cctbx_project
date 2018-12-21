from __future__ import division
import iotbx.ccp4_map
from cctbx import crystal
import sys, os

def extract_from(file_name):

  from iotbx.mrcfile.crystal_symmetry_from_ccp4_map \
        import extract_from as extract_using_mrcfile
  return extract_using_mrcfile(file_name)

  # BELOW HERE IS NOT USED ANY MORE

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
    m = iotbx.mrcfile.map_reader(file_name=file_name)
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
