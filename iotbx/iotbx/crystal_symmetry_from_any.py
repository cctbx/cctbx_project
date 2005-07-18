from cctbx import crystal
from iotbx.scalepack import crystal_symmetry_from_hkl as from_scalepack_hkl
from iotbx.xds import crystal_symmetry_from_hkl as from_xds_hkl
from iotbx.dtrek import crystal_symmetry_from_ref as from_dtrek_ref
from iotbx.mtz import crystal_symmetry_from_mtz as from_mtz
from iotbx.shelx import crystal_symmetry_from_ins as from_shelx_ins
from iotbx.cns import crystal_symmetry_from_inp as from_cns_inp
from iotbx.cns import crystal_symmetry_from_sdb as from_cns_sdb
from iotbx.pdb import crystal_symmetry_from_pdb as from_pdb
from iotbx.solve import crystal_symmetry_from_inp as from_solve_inp
from iotbx.xplor import crystal_symmetry_from_map as from_xplor_map

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
    except:
      return None
  try:
    return crystal.symmetry(unit_cell=unit_cell,space_group=space_group)
  except:
    return None

def extract_from(file_name):
  '''
     extract_from takes a file name or a string description.
     If given a file name, it attempts to open the file with all the known
     symmetry-containing file-types and extract the symmetry information.
     If all of these fail, it attempts to interpret the input string
     according to the function cyrstal_symmetry_from_any.from_string.
  '''
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
    try: return fmt.extract_from(file_name)
    except: pass
  return from_string(file_name)
