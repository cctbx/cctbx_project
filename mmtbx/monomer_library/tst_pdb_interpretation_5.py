from __future__ import absolute_import, division, print_function
from mmtbx.monomer_library import pdb_interpretation, server
from libtbx.utils import format_cpu_times
import sys

def exercise_unknown_ligand_type_energies_gap(mon_lib_srv, ener_lib):
  """
  Test for unknown multi-atom ligand resulting in shorter type_energies arrays.

  Bug: OFO ligand (unknown, multi-atom) from 1a7e causes conformer_i_seq.convert()
  to produce arrays shorter than n_atoms, breaking VDW radius lookups.
  """

  raw_records = """
CRYST1   41.885   81.170   36.685  90.00  90.00  90.00 P 21 21 21    4
ATOM   1172  N   LEU A 118      20.710  33.287  16.579  1.00  9.14           N
ATOM   1173  CA  LEU A 118      21.556  33.906  17.579  1.00 10.59           C
ATOM   1174  C   LEU A 118      21.343  33.277  18.938  1.00 11.40           C
ATOM   1175  O   LEU A 118      21.996  33.665  19.911  1.00 12.83           O
ATOM   1176  CB  LEU A 118      21.251  35.412  17.648  1.00  9.85           C
ATOM   1177  CG  LEU A 118      21.389  36.235  16.387  1.00 10.49           C
ATOM   1178  CD1 LEU A 118      21.341  37.707  16.784  1.00  9.34           C
ATOM   1179  CD2 LEU A 118      22.719  35.916  15.657  1.00 11.20           C
ATOM   1180  OXT LEU A 118      20.451  32.437  19.021  1.00 12.76           O
ATOM   1181  H   LEU A 118      19.746  33.465  16.577  1.00  0.00           H
TER    1182      LEU A 118
HETATM 1185 FE1  OFO A 119      21.944  45.048  11.891  1.00  7.00          FE
HETATM 1186  O   OFO A 119      20.987  46.213  12.935  1.00  6.74           O
HETATM 1187 FE2  OFO A 119      19.598  45.781  14.015  1.00  7.34          FE
HETATM 1188  OH  OFO A 119      18.761  47.693  13.352  1.00  8.54           O
HETATM 1189  HO  OFO A 119      19.241  48.506  13.509  1.00  0.00           H
HETATM 1190  O   HOH A 121      28.495  28.279   8.620  1.00 17.65           O
HETATM 1191  H1  HOH A 121      28.094  28.549   7.794  1.00  0.00           H
HETATM 1192  H2  HOH A 121      28.827  27.398   8.449  1.00  0.00           H
HETATM 1193  O   HOH A 122       1.095  61.206  12.261  1.00 12.09           O
END
"""

  params = pdb_interpretation.master_params.extract()
  params.clash_guard.nonbonded_distance_threshold = None

  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=params,
    raw_records=raw_records)

  all_chain_proxies = processed_pdb_file.all_chain_proxies
  n_atoms = all_chain_proxies.pdb_atoms.size()
  te_size = len(all_chain_proxies.nonbonded_energy_type_registry.symbols)
  th_size = len(all_chain_proxies.type_h_bonds)

  assert te_size == n_atoms, \
    "type_energies size=%d but n_atoms=%d" % (te_size, n_atoms)
  assert th_size == n_atoms, \
    "type_h_bonds size=%d but n_atoms=%d" % (th_size, n_atoms)

def run(args):
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  exercise_unknown_ligand_type_energies_gap(mon_lib_srv, ener_lib)
  print(format_cpu_times())

if (__name__ == "__main__"):
  run(sys.argv[1:])
