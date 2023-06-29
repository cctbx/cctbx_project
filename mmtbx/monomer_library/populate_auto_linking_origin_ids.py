from __future__ import absolute_import, division, print_function

from mmtbx.monomer_library.server import mon_lib_list_cif, geostd_list_cif
from mmtbx.monomer_library.server import merge_and_overwrite_cifs

# ['SS BOND', # short desc.
#      # complete desc.
#      'Disulphide bond for CYS-like sulphur atoms within 3A (default) using '
#      'values determined from hi-res structures and published in CCN. '
#      'Some bonds are automatically excluded based on distance from metals.',
#      # citation
#      'Comput. Cryst. Newsl. (2015), 6, 13-13.',
#      # geo file header - bond, angle, dihedral (None will suppress output)
#      ['Disulphide bridge']*3,
#      # internals
#      [0,1,2], # does not seem to be used much...

def main():
  list_cif = mon_lib_list_cif()
  geostd_list_cif_obj = geostd_list_cif()
  list_cif = merge_and_overwrite_cifs(geostd_list_cif_obj,
                                      list_cif,
                                      )
  print(dir(list_cif))
  print(dir(list_cif.cif))
  cif = list_cif.cif
  print(list(cif.blocks.keys()))
  j=0
  tmp = []
  for i, block in enumerate(cif.blocks):
    if block in ['link_list']: continue
    if block.startswith('link_'):
      j+=1
      tmp.append([])
      tmp[-1].append(block)
      tmp[-1].append('%s from Monomer Library or GeoStd' % block.replace('link_', ''))
      #tmp[-1].append('')
      #tmp[-1].append([block]*6)
      #tmp[-1].append([0,1,2,3,4,5])
    print(i, j, block)

  outl = 'standard_cif_links = [\n'
  for block in sorted(tmp):
    outl += '  %s,\n' % block
  outl += ']\n'
  print(outl)

if __name__ == '__main__':
  main()
