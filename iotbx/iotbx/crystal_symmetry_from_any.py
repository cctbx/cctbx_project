from iotbx.shelx import crystal_symmetry_from_ins
from iotbx.cns import crystal_symmetry_from_inp
from iotbx.scalepack import crystal_symmetry_from_hkl
from iotbx.mtz import crystal_symmetry_from_mtz

def extract_from(file_name):
  crystal_symmetry = None
  if (crystal_symmetry == None):
    try: crystal_symmetry = crystal_symmetry_from_ins.extract_from(
      open(file_name))
    except: pass
  if (crystal_symmetry == None):
    try: crystal_symmetry = crystal_symmetry_from_inp.extract_from(
      open(file_name))
    except: pass
  if (crystal_symmetry == None):
    try: crystal_symmetry = crystal_symmetry_from_hkl.extract_from(
      open(file_name))
    except: pass
  if (crystal_symmetry == None):
    try: crystal_symmetry = crystal_symmetry_from_mtz.extract_from(file_name)
    except: pass
  return crystal_symmetry
