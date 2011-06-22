import libtbx.load_env

from iotbx.scalepack import crystal_symmetry_from_hkl as from_scalepack_hkl
from iotbx.xds import crystal_symmetry_from_hkl as from_xds_hkl
from iotbx.dtrek import crystal_symmetry_from_ref as from_dtrek_ref
if (libtbx.env.has_module("ccp4io")):
  from iotbx.mtz import crystal_symmetry_from_mtz as from_mtz
else:
  from_mtz = None
from iotbx.shelx import crystal_symmetry_from_ins as from_shelx_ins
from iotbx.cns import crystal_symmetry_from_inp as from_cns_inp
from iotbx.cns import crystal_symmetry_from_sdb as from_cns_sdb
from iotbx.pdb import crystal_symmetry_from_pdb as from_pdb
from iotbx.solve import crystal_symmetry_from_inp as from_solve_inp
from iotbx.xplor import crystal_symmetry_from_map as from_xplor_map
from cctbx import crystal
from libtbx.path import canonical_path
import os

def from_string(string):
  '''
     Interprets a symmetry-string object in the following format:
     a, b, c, alpha, beta, gamma, space-group
     Please be aware that this is comma delimited, so you may not
     have commas within the space-group substring.
  '''
  parts = string.split(",")
  unit_cell = None
  space_group = None
  if len(parts) == 7:
    unit_cell = parts[:-1]
    space_group = parts[-1]
  elif len(parts) == 6:
    unit_cell = parts
  elif len(parts) == 1:
    space_group = parts[0]
  else:
    return None
  if unit_cell is not None:
    try:
      unit_cell = [float(number) for number in unit_cell]
    except KeyboardInterrupt: raise
    except Exception:
      return None
  try:
    return crystal.symmetry(unit_cell=unit_cell,space_group=space_group)
  except KeyboardInterrupt: raise
  except Exception:
    return None

def extract_from(file_name):
  '''
     extract_from takes a file name or a string description.
     If given a file name, it attempts to open the file with all the known
     symmetry-containing file-types and extract the symmetry information.
     If all of these fail, it attempts to interpret the input string
     according to the function cyrstal_symmetry_from_any.from_string.
     If given a crystal.symmetry it returns the crystal.symmetry.
  '''
  if type(file_name) == type(crystal.symmetry()):
    return file_name
  for fmt in (from_scalepack_hkl,
              from_xds_hkl,
              from_dtrek_ref,
              from_mtz,
              from_shelx_ins,
              from_cns_inp,
              from_cns_sdb,
              from_pdb,
              from_solve_inp,
              from_xplor_map):
    if (fmt is None): continue
    try: return fmt.extract_from(file_name)
    except KeyboardInterrupt: raise
    except Exception: pass
  return from_string(file_name)

def extract_and_append(file_names, target_list, extract_function=extract_from):
  file_names_done = set()
  for file_name in file_names:
    if (file_name is None): continue
    if (not os.path.isfile(file_name)): continue
    file_name = canonical_path(file_name)
    if (file_name in file_names_done): continue
    file_names_done.add(file_name)
    try:
      crystal_symmetry = extract_function(file_name=file_name)
    except KeyboardInterrupt: raise
    except Exception: pass
    else:
      if (crystal_symmetry is not None):
        target_list.append(crystal_symmetry)
