from iotbx.scalepack import crystal_symmetry_from_hkl as from_scalepack_hkl
from iotbx.xds import crystal_symmetry_from_hkl as from_xds_hkl
from iotbx.dtrek import crystal_symmetry_from_ref as from_dtrek_ref
from iotbx.mtz import crystal_symmetry_from_mtz as from_mtz
from iotbx.shelx import crystal_symmetry_from_ins as from_shelx_ins
from iotbx.cns import crystal_symmetry_from_inp as from_cns_inp
from iotbx.cns import crystal_symmetry_from_sdb as from_cns_sdb

def extract_from(file_name):
  for fmt in (from_scalepack_hkl,
              from_xds_hkl,
              from_dtrek_ref,
              from_mtz,
              from_shelx_ins,
              from_cns_inp,
              from_cns_sdb):
    try: return fmt.extract_from(file_name)
    except: pass
  return None
