from iotbx.scalepack import crystal_symmetry_from_hkl as from_scalepack_hkl
from iotbx.xds import crystal_symmetry_from_hkl as from_xds_hkl
from iotbx.mtz import crystal_symmetry_from_mtz as from_mtz
from iotbx.shelx import crystal_symmetry_from_ins as from_shelx_ins
from iotbx.cns import crystal_symmetry_from_inp as from_cns_inp

def extract_from(file_name):
  crystal_symmetry = None
  if (crystal_symmetry is None):
    try: crystal_symmetry = from_scalepack_hkl.extract_from(file_name)
    except: pass
  if (crystal_symmetry is None):
    try: crystal_symmetry = from_xds_hkl.extract_from(file_name)
    except: pass
  if (crystal_symmetry is None):
    try: crystal_symmetry = from_mtz.extract_from(file_name)
    except: pass
  if (crystal_symmetry is None):
    try: crystal_symmetry = from_shelx_ins.extract_from(open(file_name))
    except: pass
  if (crystal_symmetry is None):
    try: crystal_symmetry = from_cns_inp.extract_from(open(file_name))
    except: pass
  return crystal_symmetry
