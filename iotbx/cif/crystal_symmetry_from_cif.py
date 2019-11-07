from __future__ import absolute_import, division, print_function
import iotbx.cif
from iotbx.cif.builders import crystal_symmetry_builder

def extract_from(file_name=None, file=None):
  assert [file_name, file].count(None) == 1
  cif_block = list(iotbx.cif.reader(file_path=file_name,
                               file_object=file).model().values())[0] # FIXME
  builder = crystal_symmetry_builder(cif_block)
  return builder.crystal_symmetry
