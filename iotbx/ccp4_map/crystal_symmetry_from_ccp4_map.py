from __future__ import absolute_import, division, print_function

def extract_from(file_name):

  from iotbx.mrcfile.crystal_symmetry_from_ccp4_map \
        import extract_from as extract_using_mrcfile
  return extract_using_mrcfile(file_name)

