from __future__ import absolute_import, division, print_function

import iotbx

def is_monomer_library_file(file_name):
  try:
    cif_model = iotbx.cif.reader(file_path=file_name).model()
    for cif_block in cif_model.values():
      if '_chem_comp_atom' in cif_block:
        return True
  except Exception as e:
    return False
