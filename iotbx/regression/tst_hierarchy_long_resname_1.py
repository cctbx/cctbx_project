from __future__ import absolute_import, division, print_function
import time
import iotbx.pdb
from libtbx.test_utils import assert_lines_in_text

# ------------------------------------------------------------------------------

# from https://github.com/wwPDB/extended-wwPDB-identifier-examples
# https://github.com/wwPDB/extended-wwPDB-identifier-examples/blob/main/Models/7fgz-extended_CCD_code-model.cif
mmcif_str = '''
data_XXXX
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   2140 N  N   . LYS A 1 261 ? 0.399   -10.171 39.802 1.00 40.11 ? 279 LYS A N   1
ATOM   2141 C  CA  . LYS A 1 261 ? -0.169  -9.988  41.173 1.00 43.86 ? 279 LYS A CA  1
ATOM   2142 C  C   . LYS A 1 261 ? 0.687   -9.011  41.991 1.00 41.94 ? 279 LYS A C   1
ATOM   2143 O  O   . LYS A 1 261 ? 1.044   -7.920  41.556 1.00 39.32 ? 279 LYS A O   1
ATOM   2144 C  CB  . LYS A 1 261 ? -0.260  -11.336 41.902 1.00 46.47 ? 279 LYS A CB  1
ATOM   2145 C  CG  . LYS A 1 261 ? -1.583  -12.074 41.713 1.00 49.13 ? 279 LYS A CG  1
ATOM   2146 C  CD  . LYS A 1 261 ? -1.611  -13.468 42.315 1.00 51.03 ? 279 LYS A CD  1
ATOM   2147 C  CE  . LYS A 1 261 ? -2.923  -13.799 42.993 1.00 52.86 ? 279 LYS A CE  1
ATOM   2148 N  NZ  . LYS A 1 261 ? -3.209  -12.856 44.100 1.00 54.19 ? 279 LYS A NZ  1
HETATM 2149 CA CA  . CA  B 2 .   ? -17.362 -22.385 28.047 1.00 15.20 ? 301 CA  A CA  1
HETATM 2150 N  N1  . 7ZTVU C 3 .   ? -7.743  -6.355  8.243  1.00 18.72 ? 302 7ZTVU A N1  1
HETATM 2151 C  C2  . 7ZTVU C 3 .   ? -8.462  -5.534  9.265  1.00 16.68 ? 302 7ZTVU A C2  1
HETATM 2152 C  C3  . 7ZTVU C 3 .   ? -8.092  -5.865  10.711 1.00 17.35 ? 302 7ZTVU A C3  1
HETATM 2153 N  N4  . 7ZTVU C 3 .   ? -7.975  -7.334  10.767 1.00 17.04 ? 302 7ZTVU A N4  1
HETATM 2154 C  C5  . 7ZTVU C 3 .   ? -6.781  -7.689  10.027 1.00 17.11 ? 302 7ZTVU A C5  1
HETATM 2155 C  C6  . 7ZTVU C 3 .   ? -7.363  -7.730  8.633  1.00 16.97 ? 302 7ZTVU A C6  1
HETATM 2156 C  C7  . 7ZTVU C 3 .   ? -8.071  -7.951  12.064 1.00 16.77 ? 302 7ZTVU A C7  1
HETATM 2157 C  C8  . 7ZTVU C 3 .   ? -8.490  -9.440  12.025 1.00 16.60 ? 302 7ZTVU A C8  1
HETATM 2158 O  O8  . 7ZTVU C 3 .   ? -8.391  -10.140 10.781 1.00 14.14 ? 302 7ZTVU A O8  1
HETATM 2159 C  C9  . 7ZTVU C 3 .   ? -8.438  -6.311  6.935  1.00 20.13 ? 302 7ZTVU A C9  1
HETATM 2160 C  C10 . 7ZTVU C 3 .   ? -7.646  -6.965  5.796  1.00 22.62 ? 302 7ZTVU A C10 1
HETATM 2161 S  S   . 7ZTVU C 3 .   ? -8.372  -6.708  4.279  1.00 24.90 ? 302 7ZTVU A S   1
HETATM 2162 O  O1S . 7ZTVU C 3 .   ? -7.579  -7.371  3.224  1.00 25.47 ? 302 7ZTVU A O1S 1
HETATM 2163 O  O2S . 7ZTVU C 3 .   ? -8.533  -5.211  4.078  1.00 25.01 ? 302 7ZTVU A O2S 1
HETATM 2164 O  O3S . 7ZTVU C 3 .   ? -9.691  -7.326  4.282  1.00 26.42 ? 302 7ZTVU A O3S 1
HETATM 2165 N  N1  . 7ZTVU D 3 .   ? -11.233 1.732   10.446 1.00 74.61 ? 303 7ZTVU A N1  1
HETATM 2166 C  C2  . 7ZTVU D 3 .   ? -11.682 2.495   9.264  1.00 77.03 ? 303 7ZTVU A C2  1
HETATM 2167 C  C3  . 7ZTVU D 3 .   ? -11.907 3.979   9.617  1.00 76.82 ? 303 7ZTVU A C3  1
HETATM 2168 N  N4  . 7ZTVU D 3 .   ? -12.335 4.282   11.022 1.00 73.77 ? 303 7ZTVU A N4  1
HETATM 2169 C  C5  . 7ZTVU D 3 .   ? -12.654 3.123   11.908 1.00 69.90 ? 303 7ZTVU A C5  1
HETATM 2170 C  C6  . 7ZTVU D 3 .   ? -12.333 1.719   11.409 1.00 69.07 ? 303 7ZTVU A C6  1
HETATM 2171 C  C7  . 7ZTVU D 3 .   ? -11.415 5.232   11.725 1.00 71.46 ? 303 7ZTVU A C7  1
HETATM 2172 C  C8  . 7ZTVU D 3 .   ? -12.126 5.976   12.871 1.00 70.74 ? 303 7ZTVU A C8  1
HETATM 2173 O  O8  . 7ZTVU D 3 .   ? -11.360 5.921   14.096 1.00 61.89 ? 303 7ZTVU A O8  1
HETATM 2174 C  C9  . 7ZTVU D 3 .   ? -10.756 0.357   10.142 1.00 74.40 ? 303 7ZTVU A C9  1
HETATM 2175 C  C10 . 7ZTVU D 3 .   ? -10.121 -0.289  11.396 1.00 75.92 ? 303 7ZTVU A C10 1
HETATM 2176 S  S   . 7ZTVU D 3 .   ? -8.757  -1.214  11.091 1.00 77.53 ? 303 7ZTVU A S   1
HETATM 2177 O  O1S . 7ZTVU D 3 .   ? -7.735  -0.305  10.525 1.00 74.01 ? 303 7ZTVU A O1S 1
HETATM 2178 O  O2S . 7ZTVU D 3 .   ? -9.042  -2.313  10.141 1.00 83.78 ? 303 7ZTVU A O2S 1
HETATM 2179 O  O3S . 7ZTVU D 3 .   ? -8.264  -1.864  12.321 1.00 70.44 ? 303 7ZTVU A O3S 1
HETATM 2180 O  O   . HOH E 4 .   ? 10.567  -6.539  28.064 1.00 30.11 ? 401 HOH A O   1
HETATM 2181 O  O   . HOH E 4 .   ? 2.856   -8.848  40.692 1.00 32.14 ? 402 HOH A O   1
HETATM 2182 O  O   . HOH E 4 .   ? -10.577 -29.762 22.833 1.00 22.07 ? 403 HOH A O   1
HETATM 2183 O  O   . HOH E 4 .   ? 14.705  -5.437  13.988 1.00 28.06 ? 404 HOH A O   1
'''

def test1():
  """
  Test correct reading of long residue name mmCIF file into hierarchy.
  Testing output to pdb and cif.
  """
  inp = iotbx.pdb.input(lines=mmcif_str.split("\n"), source_info=None)
  h = inp.construct_hierarchy()
  chains = list(h.overall_counts().chain_ids.keys())
  chains.sort()
  answer = ['A']
  assert (chains == answer), '%s %s' % (chains, answer)
  # print(h.overall_counts().resnames)
  resnames = sorted(h.overall_counts().resnames.keys())
  assert resnames == [' CA', '7ZTVU', 'HOH', 'LYS'], resnames

  o_pdb_str = h.as_pdb_string()
  # print(o_pdb_str)
  # Note here incorrect/trimmed residue names
  # There's no way to correctly output resnames longer than 3 char in PDB format
  # print(o_pdb_str)
  assert o_pdb_str == ""
#   assert_lines_in_text(o_pdb_str, """\
# ATOM      8  CE  LYS A 279      -2.923 -13.799  42.993  1.00 52.86           C
# ATOM      9  NZ  LYS A 279      -3.209 -12.856  44.100  1.00 54.19           N
# TER
# HETATM   10 CA    CA A 301     -17.362 -22.385  28.047  1.00 15.20          CA
# HETATM   11  C10 7ZTVU A 302      -7.646  -6.965   5.796  1.00 22.62           C
# HETATM   12  C2  7ZTVU A 302      -8.462  -5.534   9.265  1.00 16.68           C
#     """)

  o_cif_str = "%s" % h.as_cif_block()
  # print(o_cif_str)
  assert_lines_in_text(o_cif_str, """\
ATOM     8  CE   .  LYS    A  279  ?   -2.92300  -13.79900  42.99300  1.000  52.86000  C   ?  A  ?  1  CE   1
ATOM     9  NZ   .  LYS    A  279  ?   -3.20900  -12.85600  44.10000  1.000  54.19000  N   ?  A  ?  1  NZ   1
HETATM  10  CA   .  CA     A  301  ?  -17.36200  -22.38500  28.04700  1.000  15.20000  CA  ?  B  ?  .  CA   1
HETATM  11  C10  .  7ZTVU  A  302  ?   -7.64600   -6.96500   5.79600  1.000  22.62000  C   ?  C  ?  .  C10  1
HETATM  12  C2   .  7ZTVU  A  302  ?   -8.46200   -5.53400   9.26500  1.000  16.68000  C   ?  C  ?  .  C2   1
HETATM  13  C3   .  7ZTVU  A  302  ?   -8.09200   -5.86500  10.71100  1.000  17.35000  C   ?  C  ?  .  C3   1
""")

def test2():
  """
  Testing selections
  """
  inp = iotbx.pdb.input(lines=mmcif_str.split("\n"), source_info=None)
  h = inp.construct_hierarchy()
  asc = h.atom_selection_cache()
  s = asc.iselection("resname 7ZTVU")
  assert list(s) == list(range(10,40)), list(s)
  s = asc.iselection("resname 7ZTVU and name N1")
  assert list(s) == [18,33], list(s)

if (__name__ == "__main__"):
  t0 = time.time()
  test1()
  test2()
  print("OK. Time: %8.3f"%(time.time()-t0))
