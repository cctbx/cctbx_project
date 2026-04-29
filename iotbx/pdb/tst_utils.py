"""Test model utilities"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb.utils

pdb_str_to_be_cif="""
CRYST1   40.339   36.116   46.266  90.00  90.00  90.00 P 1
ATOM      1  CA  ASP AXYB2      34.633  18.762  20.254  1.00 22.59           C
ATOM      2  CA  LYS AXYB3      36.047  17.704  23.610  1.00 19.79           C
ATOM      3  CA  ILE AXYB4      35.551  19.482  26.886  1.00 19.33           C
ATOM      4  CA AHIS AXYB5      38.649  21.223  28.218  0.50 19.79           C
ATOM      5  CA BHIS AXYB6      38.583  21.270  28.209  0.50 20.43           C
ATOM      6  CA  GLY A 138      38.261  15.285  27.690  1.00  6.80           C
ATOM      7  CA  ALA A 139      34.607  14.241  27.428  1.00  4.76           C
ATOM      8  CA ALEU A 140      33.091  14.490  23.937  0.50  5.08           C
ATOM      9  CA BLEU A 140      33.072  14.565  23.972  0.50  5.41           C
ATOM     10  CA  ASN A 141      30.271  17.061  23.474  1.00  5.65           C
"""

pdb_str_other="""
CRYST1   40.339   36.116   46.266  90.00  90.00  90.00 P 1
ATOM      1  CA  ASP DXYB2      34.633  18.762  20.254  1.00 22.59           C
ATOM      2  CA  LYS DXYB3      36.047  17.704  23.610  1.00 19.79           C
ATOM      3  CA  ILE DXYB4      35.551  19.482  26.886  1.00 19.33           C
ATOM      4  CA AHIS DXYB5      38.649  21.223  28.218  0.50 19.79           C
ATOM      5  CA BHIS DXYB6      38.583  21.270  28.209  0.50 20.43           C
ATOM      6  CA  GLY D 138      38.261  15.285  27.690  1.00  6.80           C
ATOM      7  CA  ALA D 139      34.607  14.241  27.428  1.00  4.76           C
ATOM      8  CA ALEU D 140      33.091  14.490  23.937  0.50  5.08           C
ATOM      9  CA BLEU D 140      33.072  14.565  23.972  0.50  5.41           C
ATOM     10  CA  ASN D 141      30.271  17.061  23.474  1.00  5.65           C
"""

pseudo_as_pdb="""
ATOM     10  CA  ASN A 141      30.271  17.061  23.474  1.00  5.65
ATOM     24 ZC1'  GC U  11      10.024   9.813   3.777  1.00 36.50      UNK
"""
spacing_as_pdb="""
ATOM     10  I   ASN A 141      30.271  17.061  23.474  1.00  5.65
ATOM     11 I    ASN A 141      30.271  17.061  23.474  1.00  5.65
HETATM   12  I   LIG A 141      30.271  17.061  23.474  1.00  5.65
HETATM   13 I    LIG A 141      30.271  17.061  23.474  1.00  5.65
"""

as_pdb = pdb_str_to_be_cif
other_as_pdb = pdb_str_other

# Convert pdb_str_hybrid_residues to mmcif:
from libtbx.test_utils import convert_pdb_to_cif_for_pdb_str
convert_pdb_to_cif_for_pdb_str(locals(), chain_addition="ZXLONG")
as_cif = pdb_str_to_be_cif # now it is cif
other_as_cif = pdb_str_other# now it is cif

def exercise_all_chain_ids():
  ids = iotbx.pdb.utils.all_chain_ids()
  assert len(ids)==3906
  assert len(set(ids))==3906

def exercise_add_models_and_hierarchies():
  from iotbx.pdb.utils import get_pdb_info, add_models, add_hierarchies
  ph1 = get_pdb_info(as_pdb).hierarchy
  ph2 = get_pdb_info(other_as_pdb).hierarchy
  m1 = ph1.as_model_manager(crystal_symmetry = None)
  m2 = ph2.as_model_manager(crystal_symmetry = None)
  ph = add_hierarchies(hierarchy_list = [ph1,ph2])
  ph2a = ph.apply_atom_selection("chain D")
  assert ph2.is_similar_hierarchy(ph2a)
  m = add_models(model_list = [m1,m2])
  ph1a = m.get_hierarchy().apply_atom_selection("chain A")
  assert ph1a.is_similar_hierarchy(ph1)


def exercise_set_element_ignoring_spacings():
  from iotbx.pdb.utils import get_pdb_info
  pdb_info = get_pdb_info(spacing_as_pdb, allow_incorrect_spacing = True)
  from iotbx.pdb.utils import set_element_ignoring_spacings
  set_element_ignoring_spacings(pdb_info.hierarchy)

def exercise_check_pseudo_atoms():
  from iotbx.pdb.utils import get_pdb_info
  pdb_info = get_pdb_info(pseudo_as_pdb, check_pseudo = True)
  from iotbx.pdb.utils import check_for_pseudo_atoms
  check_for_pseudo_atoms(pdb_info.hierarchy)

def exercise_get_pdb_info():
  from iotbx.pdb.utils import get_pdb_info
  pdb_info_from_pdb = get_pdb_info(as_pdb)
  pdb_info_from_cif = get_pdb_info(as_cif)
  assert not pdb_info_from_pdb.hierarchy.is_similar_hierarchy(
    pdb_info_from_cif.hierarchy)
  for model in pdb_info_from_cif.hierarchy.models():
    for chain in model.chains():
      chain.id = chain.id.replace("ZXLONG","") # make it short again
  assert pdb_info_from_pdb.hierarchy.is_similar_hierarchy(
    pdb_info_from_cif.hierarchy)
  assert pdb_info_from_pdb.crystal_symmetry.is_similar_symmetry(
    pdb_info_from_cif.crystal_symmetry)

def exercise_interleave_alt_confs():
  from scitbx.array_family import flex
  import iotbx.pdb
  pdb_inp_lines_1 = flex.split_lines("""\
CRYST1   14.600   26.100   29.200  90.00  90.00  90.00 P 21 21 21
SCALE1      0.068493 -0.000000 -0.000000        0.00000
SCALE2      0.000000  0.038314 -0.000000        0.00000
SCALE3      0.000000  0.000000  0.034247        0.00000
ATOM     17  N  ASER     4      -0.155   3.125   4.014  1.00 17.55           N
ATOM     18  CA ASER     4      -0.175   1.896   4.797  1.00 15.51           C
ATOM     19  C  ASER     4       1.158   1.683   5.505  1.00 16.48           C
ATOM     20  O  ASER     4       1.508   0.560   5.868  1.00  6.79           O
ATOM     21  CB ASER     4      -0.490   0.698   3.899  1.00 17.79           C
ATOM     22  OG ASER     4      -1.752   0.849   3.272  0.50 11.67           O
ATOM     24  N  ALEU     5       1.898   2.771   5.698  1.00 10.29           N
ATOM     25  CA ALEU     5       3.208   2.705   6.333  1.00  3.27           C
ATOM     26  C  ALEU     5       3.407   3.868   7.301  1.00  7.66           C
ATOM     27  O  ALEU     5       3.458   5.024   6.884  1.00 14.56           O
ATOM     28  CB ALEU     5       4.310   2.711   5.272  1.00  3.87           C
ATOM     29  CG ALEU     5       5.736   2.474   5.769  1.00 14.70           C
ATOM     30  CD1ALEU     5       5.818   1.155   6.516  1.00 19.83           C
ATOM     31  CD2ALEU     5       6.719   2.499   4.611  1.00 13.84           C
  """)

  pdb_inp_lines_2 = flex.split_lines("""\
CRYST1   14.600   26.100   29.200  90.00  90.00  90.00 P 21 21 21
SCALE1      0.068493 -0.000000 -0.000000        0.00000
SCALE2      0.000000  0.038314 -0.000000        0.00000
SCALE3      0.000000  0.000000  0.034247        0.00000
ATOM     17  N  BSER     4      -0.155   3.125   4.014  1.00 17.55           N
ATOM     18  CA BSER     4      -0.175   1.896   4.797  1.00 15.51           C
ATOM     19  C  BSER     4       1.158   1.683   5.505  1.00 16.48           C
ATOM     20  O  BSER     4       1.508   0.560   5.868  1.00  6.79           O
ATOM     21  CB BSER     4      -0.490   0.698   3.899  1.00 17.79           C
ATOM     23  OG BSER     4       0.484   0.555   2.880  0.50  2.71           O
ATOM     24  N  BLEU     5       1.898   2.771   5.698  1.00 10.29           N
ATOM     25  CA BLEU     5       3.208   2.705   6.333  1.00  3.27           C
ATOM     26  C  BLEU     5       3.407   3.868   7.301  1.00  7.66           C
ATOM     27  O  BLEU     5       3.458   5.024   6.884  1.00 14.56           O
ATOM     28  CB BLEU     5       4.310   2.711   5.272  1.00  3.87           C
ATOM     29  CG BLEU     5       5.736   2.474   5.769  1.00 14.70           C
ATOM     30  CD1BLEU     5       5.818   1.155   6.516  1.00 19.83           C
ATOM     31  CD2BLEU     5       6.719   2.499   4.611  1.00 13.84           C
""")

  h1 = iotbx.pdb.input(source_info=None, lines=pdb_inp_lines_1).construct_hierarchy()
  h2 = iotbx.pdb.input(source_info=None, lines=pdb_inp_lines_2).construct_hierarchy()

  # Interleave the alt confs
  from iotbx.pdb.utils import interleave_alt_confs
  new_h = interleave_alt_confs(h1,h2)

  # Make sure we got both confs
  new_h1 = new_h.apply_atom_selection('altloc A')
  new_h2 = new_h.apply_atom_selection('altloc B')
  assert new_h1.as_pdb_string() == h1.as_pdb_string()
  assert new_h2.as_pdb_string() == h2.as_pdb_string()

  assert new_h1.is_similar_hierarchy(h1)
  assert new_h2.is_similar_hierarchy(h2)

  assert not new_h1.is_similar_hierarchy(h2)
  assert not new_h2.is_similar_hierarchy(h1)

def run():
  exercise_add_models_and_hierarchies()
  exercise_set_element_ignoring_spacings()
  exercise_check_pseudo_atoms()
  exercise_get_pdb_info()
  exercise_all_chain_ids()
  exercise_interleave_alt_confs()
  print("OK")

if __name__ == '__main__':
  run()
