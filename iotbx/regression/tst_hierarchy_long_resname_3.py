from __future__ import absolute_import, division, print_function
import time
# from libtbx.test_utils import assert_lines_in_text
import iotbx.pdb
import iotbx.cif
import mmtbx.model
from mmtbx import monomer_library
from libtbx.utils import Sorry

# ------------------------------------------------------------------------------

# from https://github.com/wwPDB/extended-wwPDB-identifier-examples
# https://github.com/wwPDB/extended-wwPDB-identifier-examples/blob/main/Models/7fgz-extended_CCD_code-model.cif
model_cif = '''
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
'''

ligand_cif = """\
data_7ZTVU
#

_chem_comp.id                                   7ZTVU
_chem_comp.name                                 "4-(2-HYDROXYETHYL)-1-PIPERAZINE ETHANESULFONIC ACID"
_chem_comp.type                                 NON-POLYMER
_chem_comp.pdbx_type                            HETAIN
_chem_comp.formula                              "C8 H18 N2 O4 S"
_chem_comp.mon_nstd_parent_comp_id              ?
_chem_comp.pdbx_synonyms                        HEPES
_chem_comp.pdbx_formal_charge                   0
_chem_comp.pdbx_initial_date                    1999-07-08
_chem_comp.pdbx_modified_date                   2020-06-17
_chem_comp.pdbx_ambiguous_flag                  N
_chem_comp.pdbx_release_status                  REL
_chem_comp.pdbx_replaced_by                     ?
_chem_comp.pdbx_replaces                        ?
_chem_comp.formula_weight                       238.305
_chem_comp.one_letter_code                      ?
_chem_comp.three_letter_code                    7ZTVU
_chem_comp.pdbx_model_coordinates_details       ?
_chem_comp.pdbx_model_coordinates_missing_flag  N
_chem_comp.pdbx_ideal_coordinates_details       ?
_chem_comp.pdbx_ideal_coordinates_missing_flag  N
_chem_comp.pdbx_model_coordinates_db_code       1CXQ
_chem_comp.pdbx_subcomponent_list               ?
_chem_comp.pdbx_processing_site                 RCSB
#   #
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.pdbx_align
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
_chem_comp_atom.model_Cartn_x
_chem_comp_atom.model_Cartn_y
_chem_comp_atom.model_Cartn_z
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
_chem_comp_atom.pdbx_component_atom_id
_chem_comp_atom.pdbx_component_comp_id
_chem_comp_atom.pdbx_ordinal
7ZTVU  N1    N1    N  0  1  N  N  N  61.652  47.337  51.616  -0.820  -0.330   0.258  N1    7ZTVU   1
7ZTVU  C2    C2    C  0  1  N  N  N  63.051  46.905  51.967  -1.669   0.557   1.061  C2    7ZTVU   2
7ZTVU  C3    C3    C  0  1  N  N  N  63.962  47.257  50.796  -1.693   0.072   2.511  C3    7ZTVU   3
7ZTVU  N4    N4    N  0  1  N  N  N  63.551  46.544  49.578  -0.342   0.147   3.079  N4    7ZTVU   4
7ZTVU  C5    C5    C  0  1  N  N  N  62.151  46.910  49.248   0.507  -0.740   2.276  C5    7ZTVU   5
7ZTVU  C6    C6    C  0  1  N  N  N  61.218  46.554  50.402   0.531  -0.255   0.826  C6    7ZTVU   6
7ZTVU  C7    C7    C  0  1  N  N  N  64.368  46.893  48.385  -0.416  -0.441   4.423  C7    7ZTVU   7
7ZTVU  C8    C8    C  0  1  N  N  N  65.854  46.721  48.613   0.966  -0.395   5.076  C8    7ZTVU   8
7ZTVU  O8    O8    O  0  1  N  N  N  66.240  45.400  48.979   1.402   0.961   5.174  O8    7ZTVU   9
7ZTVU  C9    C9    C  0  1  N  N  N  60.685  46.992  52.709  -0.746   0.258  -1.085  C9    7ZTVU  10
7ZTVU  C10   C10   C  0  1  N  N  N  60.870  47.974  53.893   0.128  -0.621  -1.980  C10   7ZTVU  11
7ZTVU  S     S     S  0  1  N  N  N  59.559  47.692  55.090   0.219   0.104  -3.640  S     7ZTVU  12
7ZTVU  O1S   O1S   O  0  1  N  N  N  59.732  48.694  56.115  -1.019  -0.001  -4.327  O1S   7ZTVU  13
7ZTVU  O2S   O2S   O  0  1  N  N  N  59.737  46.327  55.596   0.986   1.301  -3.646  O2S   7ZTVU  14
7ZTVU  O3S   O3S   O  0  1  N  N  N  58.331  47.863  54.358   1.110  -0.870  -4.397  O3S   7ZTVU  15
7ZTVU  H21   1H2   H  0  1  N  N  N  63.405  47.334  52.932  -2.683   0.551   0.659  H21   7ZTVU  16
7ZTVU  H22   2H2   H  0  1  N  N  N  63.111  45.828  52.251  -1.271   1.571   1.024  H22   7ZTVU  17
7ZTVU  H31   1H3   H  0  1  N  N  N  64.013  48.358  50.630  -2.044  -0.959   2.543  H31   7ZTVU  18
7ZTVU  H32   2H3   H  0  1  N  N  N  65.034  47.072  51.038  -2.367   0.702   3.093  H32   7ZTVU  19
7ZTVU  H51   1H5   H  0  1  N  N  N  62.060  47.982  48.958   1.520  -0.734   2.677  H51   7ZTVU  20
7ZTVU  H52   2H5   H  0  1  N  N  N  61.817  46.449  48.288   0.109  -1.754   2.313  H52   7ZTVU  21
7ZTVU  H61   1H6   H  0  1  N  N  N  61.169  45.456  50.592   1.204  -0.885   0.244  H61   7ZTVU  22
7ZTVU  H62   2H6   H  0  1  N  N  N  60.143  46.711  50.150   0.881   0.775   0.794  H62   7ZTVU  23
7ZTVU  H71   1H7   H  0  1  N  N  N  64.033  46.316  47.491  -0.748  -1.476   4.348  H71   7ZTVU  24
7ZTVU  H72   2H7   H  0  1  N  N  N  64.136  47.925  48.032  -1.123   0.124   5.029  H72   7ZTVU  25
7ZTVU  H81   1H8   H  0  1  N  N  N  66.427  47.061  47.719   1.673  -0.961   4.470  H81   7ZTVU  26
7ZTVU  H82   2H8   H  0  1  N  N  N  66.220  47.457  49.365   0.911  -0.832   6.073  H82   7ZTVU  27
7ZTVU  HO8   HO8   H  0  1  N  N  N  67.172  45.292  49.122   2.275   0.945   5.590  HO8   7ZTVU  28
7ZTVU  H91   1H9   H  0  1  N  N  N  59.631  46.964  52.345  -1.748   0.323  -1.509  H91   7ZTVU  29
7ZTVU  H92   2H9   H  0  1  N  N  N  60.772  45.927  53.027  -0.313   1.256  -1.022  H92   7ZTVU  30
7ZTVU  H101  1H10  H  0  0  N  N  N  61.885  47.904  54.349   1.131  -0.687  -1.557  H101  7ZTVU  31
7ZTVU  H102  2H10  H  0  0  N  N  N  60.921  49.037  53.560  -0.304  -1.620  -2.044  H102  7ZTVU  32
7ZTVU  HOS3  3HOS  H  0  0  N  N  N  57.638  47.714  54.990   1.194  -0.529  -5.298  HOS3  7ZTVU  33
#   #
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
_chem_comp_bond.pdbx_ordinal
7ZTVU  N1   C2    SING  N  N   1
7ZTVU  N1   C6    SING  N  N   2
7ZTVU  N1   C9    SING  N  N   3
7ZTVU  C2   C3    SING  N  N   4
7ZTVU  C2   H21   SING  N  N   5
7ZTVU  C2   H22   SING  N  N   6
7ZTVU  C3   N4    SING  N  N   7
7ZTVU  C3   H31   SING  N  N   8
7ZTVU  C3   H32   SING  N  N   9
7ZTVU  N4   C5    SING  N  N  10
7ZTVU  N4   C7    SING  N  N  11
7ZTVU  C5   C6    SING  N  N  12
7ZTVU  C5   H51   SING  N  N  13
7ZTVU  C5   H52   SING  N  N  14
7ZTVU  C6   H61   SING  N  N  15
7ZTVU  C6   H62   SING  N  N  16
7ZTVU  C7   C8    SING  N  N  17
7ZTVU  C7   H71   SING  N  N  18
7ZTVU  C7   H72   SING  N  N  19
7ZTVU  C8   O8    SING  N  N  20
7ZTVU  C8   H81   SING  N  N  21
7ZTVU  C8   H82   SING  N  N  22
7ZTVU  O8   HO8   SING  N  N  23
7ZTVU  C9   C10   SING  N  N  24
7ZTVU  C9   H91   SING  N  N  25
7ZTVU  C9   H92   SING  N  N  26
7ZTVU  C10  S     SING  N  N  27
7ZTVU  C10  H101  SING  N  N  28
7ZTVU  C10  H102  SING  N  N  29
7ZTVU  S    O1S   DOUB  N  N  30
7ZTVU  S    O2S   DOUB  N  N  31
7ZTVU  S    O3S   SING  N  N  32
7ZTVU  O3S  HOS3  SING  N  N  33
#   #
loop_
_pdbx_chem_comp_descriptor.comp_id
_pdbx_chem_comp_descriptor.type
_pdbx_chem_comp_descriptor.program
_pdbx_chem_comp_descriptor.program_version
_pdbx_chem_comp_descriptor.descriptor
7ZTVU  SMILES            ACDLabs               10.04  "O=S(=O)(O)CCN1CCN(CCO)CC1"
7ZTVU  SMILES_CANONICAL  CACTVS                3.341  "OCCN1CCN(CC1)CC[S](O)(=O)=O"
7ZTVU  SMILES            CACTVS                3.341  "OCCN1CCN(CC1)CC[S](O)(=O)=O"
7ZTVU  SMILES_CANONICAL  "OpenEye OEToolkits"  1.5.0  "C1CN(CCN1CCO)CCS(=O)(=O)O"
7ZTVU  SMILES            "OpenEye OEToolkits"  1.5.0  "C1CN(CCN1CCO)CCS(=O)(=O)O"
7ZTVU  InChI             InChI                 1.03   "InChI=1S/C8H18N2O4S/c11-7-5-9-1-3-10(4-2-9)6-8-15(12,13)14/h11H,1-8H2,(H,12,13,14)"
7ZTVU  InChIKey          InChI                 1.03   JKMHFZQWWAIEOD-UHFFFAOYSA-N
#   #
loop_
_pdbx_chem_comp_identifier.comp_id
_pdbx_chem_comp_identifier.type
_pdbx_chem_comp_identifier.program
_pdbx_chem_comp_identifier.program_version
_pdbx_chem_comp_identifier.identifier
7ZTVU  "SYSTEMATIC NAME"  ACDLabs               10.04  "2-[4-(2-hydroxyethyl)piperazin-1-yl]ethanesulfonic acid"
7ZTVU  "SYSTEMATIC NAME"  "OpenEye OEToolkits"  1.5.0  "2-[4-(2-hydroxyethyl)piperazin-1-yl]ethanesulfonic acid"
#   #
loop_
_pdbx_chem_comp_audit.comp_id
_pdbx_chem_comp_audit.action_type
_pdbx_chem_comp_audit.date
_pdbx_chem_comp_audit.processing_site
7ZTVU  "Create component"   2023-01-01  RCSB
#
_pdbx_chem_comp_synonyms.ordinal     1
_pdbx_chem_comp_synonyms.comp_id     7ZTVU
_pdbx_chem_comp_synonyms.name        HEPES
_pdbx_chem_comp_synonyms.provenance  ?
_pdbx_chem_comp_synonyms.type        ?
"""

def test1():
  """
  Creating restraints for long residue name
  Not finished yet.
  """
  #dumping to disk if needed:
  # for name, s in [('model', model_cif), ('restr', ligand_cif)]:
  #   with open('%s.cif' % name, 'w') as f:
  #     f.write(s)
  inp = iotbx.pdb.input(lines=model_cif.split("\n"), source_info=None)
  cif_object = iotbx.cif.reader(input_string = ligand_cif).model()
  cif_objects = [('bla.cif', cif_object)]
  model = mmtbx.model.manager(model_input = inp, restraint_objects = cif_objects)
  try:
    model.process(make_restraints=True)
    geo_str = model.restraints_as_geo()
  except Sorry as e:
    pass
  # print(geo_str)
  # for l in [
  #     'bond pdb=" NZ  LYS A 279 "',
  #     '     pdb=" O   longHOH A 401 "',
  #     'nonbonded pdb=" CD  LYS A 279 "',
  #     '          pdb=" O   longHOH A 401 "']:
  #   assert_lines_in_text(geo_str, l)
  # model_cif = model.model_as_mmcif()
  # # print(model_cif)
  # for l in [
  #     'longHOH H2 O H1 103.91 3.000',
  #     'data_comp_longHOH',
  #     '   HETATM 10 O . longHOH A 401 ? -3.20900 -12.85600 46.10000 1.000 30.11000 O ? B ? . 1']:
  #   assert_lines_in_text(model_cif, l)
  # model_pdb = model.model_as_pdb()
  # # print(model_pdb)
  # for l in [
  #     'LINK         NZ  LYS A 279                 O   longHOH A 401 ',
  #     'HETATM   10  O   lon A 401      -3.209 -12.856  46.100  1.00 30.11           O']:
  #   assert_lines_in_text(model_pdb, l)

if (__name__ == "__main__"):
  t0 = time.time()
  try:
    mon_lib_srv = monomer_library.server.server()
  except: # intentional
    print("Can not initialize monomer_library, skipping test.")
  if mon_lib_srv is not None:
    test1()
  print("OK. Time: %8.3f"%(time.time()-t0))
