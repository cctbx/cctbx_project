
from __future__ import absolute_import, division, print_function
import time
from libtbx.test_utils import assert_lines_in_text
import iotbx.pdb
import iotbx.cif
import mmtbx.model
from mmtbx import monomer_library
from libtbx.utils import Sorry
from six.moves import cStringIO as StringIO

# ------------------------------------------------------------------------------

# from 1yjp, Gln4 replaced with longGLY, restraints adapted form monomer library
mm_cif = '''
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
ATOM   1  N N   . GLY A 1 1 ? -9.009  4.612  6.102  1.00 16.77 ? 1  GLY A N   1
ATOM   2  C CA  . GLY A 1 1 ? -9.052  4.207  4.651  1.00 16.57 ? 1  GLY A CA  1
ATOM   3  C C   . GLY A 1 1 ? -8.015  3.140  4.419  1.00 16.16 ? 1  GLY A C   1
ATOM   4  O O   . GLY A 1 1 ? -7.523  2.521  5.381  1.00 16.78 ? 1  GLY A O   1
ATOM   5  N N   . ASN A 1 2 ? -7.656  2.923  3.155  1.00 15.02 ? 2  ASN A N   1
ATOM   6  C CA  . ASN A 1 2 ? -6.522  2.038  2.831  1.00 14.10 ? 2  ASN A CA  1
ATOM   7  C C   . ASN A 1 2 ? -5.241  2.537  3.427  1.00 13.13 ? 2  ASN A C   1
ATOM   8  O O   . ASN A 1 2 ? -4.978  3.742  3.426  1.00 11.91 ? 2  ASN A O   1
ATOM   9  C CB  . ASN A 1 2 ? -6.346  1.881  1.341  1.00 15.38 ? 2  ASN A CB  1
ATOM   10 C CG  . ASN A 1 2 ? -7.584  1.342  0.692  1.00 14.08 ? 2  ASN A CG  1
ATOM   11 O OD1 . ASN A 1 2 ? -8.025  0.227  1.016  1.00 17.46 ? 2  ASN A OD1 1
ATOM   12 N ND2 . ASN A 1 2 ? -8.204  2.155  -0.169 1.00 11.72 ? 2  ASN A ND2 1
ATOM   13 N N   . ASN A 1 3 ? -4.438  1.590  3.905  1.00 12.26 ? 3  ASN A N   1
ATOM   14 C CA  . ASN A 1 3 ? -3.193  1.904  4.589  1.00 11.74 ? 3  ASN A CA  1
ATOM   15 C C   . ASN A 1 3 ? -1.955  1.332  3.895  1.00 11.10 ? 3  ASN A C   1
ATOM   16 O O   . ASN A 1 3 ? -1.872  0.119  3.648  1.00 10.42 ? 3  ASN A O   1
ATOM   17 C CB  . ASN A 1 3 ? -3.259  1.378  6.042  1.00 12.15 ? 3  ASN A CB  1
ATOM   18 C CG  . ASN A 1 3 ? -2.006  1.739  6.861  1.00 12.82 ? 3  ASN A CG  1
ATOM   19 O OD1 . ASN A 1 3 ? -1.702  2.925  7.072  1.00 15.05 ? 3  ASN A OD1 1
ATOM   20 N ND2 . ASN A 1 3 ? -1.271  0.715  7.306  1.00 13.48 ? 3  ASN A ND2 1
ATOM   21 N N   . longGLY A 1 4 ? -1.005  2.228  3.598  1.00 10.29 ? 4  longGLY A N   1
ATOM   22 C CA  . longGLY A 1 4 ? 0.384   1.888  3.199  1.00 10.53 ? 4  longGLY A CA  1
ATOM   23 C C   . longGLY A 1 4 ? 1.435   2.606  4.088  1.00 10.24 ? 4  longGLY A C   1
ATOM   24 O O   . longGLY A 1 4 ? 1.547   3.843  4.115  1.00 8.86  ? 4  longGLY A O   1
ATOM   30 N N   . GLN A 1 5 ? 2.154   1.821  4.871  1.00 10.38 ? 5  GLN A N   1
ATOM   31 C CA  . GLN A 1 5 ? 3.270   2.361  5.640  1.00 11.39 ? 5  GLN A CA  1
ATOM   32 C C   . GLN A 1 5 ? 4.594   1.768  5.172  1.00 11.52 ? 5  GLN A C   1
ATOM   33 O O   . GLN A 1 5 ? 4.768   0.546  5.054  1.00 12.05 ? 5  GLN A O   1
ATOM   34 C CB  . GLN A 1 5 ? 3.056   2.183  7.147  1.00 11.96 ? 5  GLN A CB  1
ATOM   35 C CG  . GLN A 1 5 ? 1.829   2.950  7.647  1.00 10.81 ? 5  GLN A CG  1
ATOM   36 C CD  . GLN A 1 5 ? 1.344   2.414  8.954  1.00 13.10 ? 5  GLN A CD  1
ATOM   37 O OE1 . GLN A 1 5 ? 0.774   1.325  9.002  1.00 10.65 ? 5  GLN A OE1 1
ATOM   38 N NE2 . GLN A 1 5 ? 1.549   3.187  10.039 1.00 12.30 ? 5  GLN A NE2 1
ATOM   39 N N   . ASN A 1 6 ? 5.514   2.664  4.856  1.00 11.99 ? 6  ASN A N   1
ATOM   40 C CA  . ASN A 1 6 ? 6.831   2.310  4.318  1.00 12.30 ? 6  ASN A CA  1
ATOM   41 C C   . ASN A 1 6 ? 7.854   2.761  5.324  1.00 13.40 ? 6  ASN A C   1
ATOM   42 O O   . ASN A 1 6 ? 8.219   3.943  5.374  1.00 13.92 ? 6  ASN A O   1
ATOM   43 C CB  . ASN A 1 6 ? 7.065   3.016  2.993  1.00 12.13 ? 6  ASN A CB  1
ATOM   44 C CG  . ASN A 1 6 ? 5.961   2.735  2.003  1.00 12.77 ? 6  ASN A CG  1
ATOM   45 O OD1 . ASN A 1 6 ? 5.798   1.604  1.551  1.00 14.27 ? 6  ASN A OD1 1
ATOM   46 N ND2 . ASN A 1 6 ? 5.195   3.747  1.679  1.00 10.07 ? 6  ASN A ND2 1
ATOM   47 N N   . TYR A 1 7 ? 8.292   1.817  6.147  1.00 14.70 ? 7  TYR A N   1
ATOM   48 C CA  . TYR A 1 7 ? 9.159   2.144  7.299  1.00 15.18 ? 7  TYR A CA  1
ATOM   49 C C   . TYR A 1 7 ? 10.603  2.331  6.885  1.00 15.91 ? 7  TYR A C   1
ATOM   50 O O   . TYR A 1 7 ? 11.041  1.811  5.855  1.00 15.76 ? 7  TYR A O   1
ATOM   51 C CB  . TYR A 1 7 ? 9.061   1.065  8.369  1.00 15.35 ? 7  TYR A CB  1
ATOM   52 C CG  . TYR A 1 7 ? 7.665   0.929  8.902  1.00 14.45 ? 7  TYR A CG  1
ATOM   53 C CD1 . TYR A 1 7 ? 6.771   0.021  8.327  1.00 15.68 ? 7  TYR A CD1 1
ATOM   54 C CD2 . TYR A 1 7 ? 7.210   1.756  9.920  1.00 14.80 ? 7  TYR A CD2 1
ATOM   55 C CE1 . TYR A 1 7 ? 5.480   -0.094 8.796  1.00 13.46 ? 7  TYR A CE1 1
ATOM   56 C CE2 . TYR A 1 7 ? 5.904   1.649  10.416 1.00 14.33 ? 7  TYR A CE2 1
ATOM   57 C CZ  . TYR A 1 7 ? 5.047   0.729  9.831  1.00 15.09 ? 7  TYR A CZ  1
ATOM   58 O OH  . TYR A 1 7 ? 3.766   0.589  10.291 1.00 14.39 ? 7  TYR A OH  1
ATOM   59 O OXT . TYR A 1 7 ? 11.358  2.999  7.612  1.00 17.49 ? 7  TYR A OXT 1
'''

ligand_cif = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
_chem_comp.initial_date
_chem_comp.modified_date
_chem_comp.source
 longGLY longGLY 'GLYCINE                             ' L-peptide 7 4 . 2009-08-12 2012-12-06
;
Copy of CCP4 Monomer Library entry.
Added neutron distances
;
#
data_comp_longGLY
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
 longGLY           N      N    NH1      -0.204
 longGLY           H      H    HNH1      0.204
 longGLY           CA     C    CH2       0.002
 longGLY           HA1    H    HCH2      0.051
 longGLY           HA2    H    HCH2      0.051
 longGLY           C      C    C         0.318
 longGLY           O      O    O        -0.422
loop_
_chem_comp_tree.comp_id
_chem_comp_tree.atom_id
_chem_comp_tree.atom_back
_chem_comp_tree.atom_forward
_chem_comp_tree.connect_type
 longGLY      N      n/a    CA     START
 longGLY      H      N      .      .
 longGLY      CA     N      C      .
 longGLY      HA1    CA     .      .
 longGLY      HA2    CA     .      .
 longGLY      C      CA     .      END
 longGLY      O      C      .      .
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 longGLY      N      H         coval       0.860    0.020    1.020
 longGLY      N      CA        coval       1.451    0.016    1.451
 longGLY      CA     HA1       coval       0.970    0.020    1.090
 longGLY      CA     HA2       coval       0.970    0.020    1.090
 longGLY      CA     C         coval       1.516    0.018    1.516
 longGLY      C      O         coval       1.231    0.020    1.231
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 longGLY      H      N      CA      114.000    3.000
 longGLY      HA1    CA     HA2     109.000    3.000
 longGLY      HA2    CA     C       109.000    3.000
 longGLY      HA1    CA     C       109.000    3.000
 longGLY      N      CA     HA1     110.000    3.000
 longGLY      N      CA     HA2     110.000    3.000
 longGLY      N      CA     C       113.300    2.900
 longGLY      CA     C      O       120.800    2.100
"""

def test1():
  """
  Creating restraints for long residue name
  """
  #dumping to disk if needed:
  # for name, s in [('model', mm_cif), ('restr', ligand_cif)]:
  #   with open('%s.cif' % name, 'w') as f:
  #     f.write(s)
  inp = iotbx.pdb.input(lines=mm_cif.split("\n"), source_info=None)
  cif_object = iotbx.cif.reader(input_string = ligand_cif).model()
  cif_objects = [('bla.cif', cif_object)]
  model = mmtbx.model.manager(
      model_input = inp,
      restraint_objects = cif_objects)
  try:
    model.process(make_restraints=True)
    geo_str = model.restraints_as_geo()
  except Sorry as e:
    pass
  # print(geo_str)
  for l in [
    'bond pdb=" C   ASN A   3 "',
    '     pdb=" N   longGLY A   4 "',
    'nonbonded pdb=" C   longGLY A   4 "',
    '          pdb=" CD  GLN A   5 "',
    ]:
    # print(geo_str.find(l),l)
    assert_lines_in_text(geo_str, l)
  model_cif = model.model_as_mmcif()
  # print(model_cif)
  for l in [
    '   ATOM 21 N . longGLY A 4 ? -1.00500 2.22800 3.59800 1.000 10.29000 N ? B ? . N 1',
    'longGLY N H coval 0.860 0.020 1.020',
    'data_comp_longGLY',
    ]:
    assert_lines_in_text(model_cif, l)
  model_pdb = model.model_as_pdb()
  # print(model_pdb)
  assert not model.can_be_output_as_pdb()

def test2():
  """
  Try creating restraints for long residue name without cif.
  Check error message formatting in
  Number of atoms with unknown nonbonded energy
  (atom.quote() function)
  """
  #dumping to disk if needed:
  # for name, s in [('model', mm_cif), ('restr', ligand_cif)]:
  #   with open('%s.cif' % name, 'w') as f:
  #     f.write(s)
  inp = iotbx.pdb.input(lines=mm_cif.split("\n"), source_info=None)
  # cif_object = iotbx.cif.reader(input_string = ligand_cif).model()
  # cif_objects = [('bla.cif', cif_object)]
  mlog = StringIO()
  model = mmtbx.model.manager(
      model_input = inp,
      log = mlog)
  try:
    model.process(make_restraints=True)
  except Sorry as e:
    mlog_txt = mlog.getvalue()
    # print(str(e))
    print(mlog_txt)
    assert_lines_in_text(mlog_txt, """ "ATOM     22  CA  longGLY A   4 .*.     C  " """)
    assert_lines_in_text(mlog_txt, """ Unknown residues: {'longGLY': 1} """)


if (__name__ == "__main__"):
  t0 = time.time()
  mon_lib_srv = None
  try:
    mon_lib_srv = monomer_library.server.server()
  except: # intentional
    print("Can not initialize monomer_library, skipping test.")
  if mon_lib_srv is not None:
    test1()
    test2()
  print("OK. Time: %8.3f"%(time.time()-t0))
