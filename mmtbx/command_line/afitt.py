# LIBTBX_SET_DISPATCHER_NAME mmtbx.afitt

from __future__ import absolute_import, division, print_function
import sys
from mmtbx.geometry_restraints import afitt

if __name__ == "__main__":
  if len(sys.argv[1:])<3:
    print('''
usage: mmtbx.afitt [-h] [-ff FF] pdb_file cif_file ligand_names

positional arguments:
  pdb_file      pdb file
  cif_file      cif file
  ligand_names  3-letter ligand names separated by commas

optional arguments:
  -h, --help    show this help message and exit
  -ff FF        afitt theory: mmff94, mmff94s pm3 or am1
    ''')
  else:
    afitt.run2()
