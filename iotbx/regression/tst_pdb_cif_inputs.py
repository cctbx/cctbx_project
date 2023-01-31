from __future__ import absolute_import, division, print_function
import sys
import iotbx.pdb
from iotbx.pdb.mmcif import cif_input
import inspect
import six

t_pdb_str = """\
ATOM      1  N   ALA     2       1.643  -2.366  -1.408  1.00
ATOM      3  CA  ALA     2       1.280  -3.608  -2.069  1.00
ATOM      6  CB  ALA     2       1.361  -4.762  -1.068  1.00
ATOM     10  C   ALA     2      -0.114  -3.466  -2.684  1.00
ATOM     11  O   ALA     2      -0.327  -3.827  -3.840  1.00
"""
t_cif_str = """\
data_1UCS
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
ATOM   1    N N    . ASN A 1 1  ? 12.388 -5.203 29.298 1.00 13.23 ? 1   ASN A N    1
ATOM   2    C CA   . ASN A 1 1  ? 12.398 -3.931 28.500 1.00 6.85  ? 1   ASN A CA   1
ATOM   3    C C    A ASN A 1 1  ? 12.384 -4.463 27.091 0.45 6.81  ? 1   ASN A C    1
ATOM   4    C C    B ASN A 1 1  ? 12.891 -4.205 27.064 0.55 9.54  ? 1   ASN A C    1
ATOM   5    O O    A ASN A 1 1  ? 12.645 -5.635 26.750 0.45 8.28  ? 1   ASN A O    1
ATOM   6    O O    B ASN A 1 1  ? 14.004 -4.638 26.750 0.55 12.04 ? 1   ASN A O    1
ATOM   7    C CB   . ASN A 1 1  ? 13.493 -3.033 29.043 1.00 7.32  ? 1   ASN A CB   1
ATOM   8    C CG   . ASN A 1 1  ? 13.414 -1.559 28.658 1.00 7.93  ? 1   ASN A CG   1
ATOM   9    O OD1  . ASN A 1 1  ? 14.388 -0.782 28.842 1.00 10.99 ? 1   ASN A OD1  1
ATOM   10   N ND2  . ASN A 1 1  ? 12.315 -1.201 28.134 1.00 11.12 ? 1   ASN A ND2  1
ATOM   11   H H1   . ASN A 1 1  ? 11.748 -5.740 28.993 1.00 19.85 ? 1   ASN A H1   1
ATOM   12   H H2   . ASN A 1 1  ? 12.235 -5.016 30.155 1.00 19.85 ? 1   ASN A H2   1
ATOM   13   H H3   . ASN A 1 1  ? 13.177 -5.609 29.220 1.00 19.85 ? 1   ASN A H3   1
ATOM   14   H HA   . ASN A 1 1  ? 11.542 -3.479 28.654 1.00 8.22  ? 1   ASN A HA   1
ATOM   15   H HB2  . ASN A 1 1  ? 13.482 -3.094 30.011 1.00 8.78  ? 1   ASN A HB2  1
ATOM   16   H HB3  . ASN A 1 1  ? 14.347 -3.380 28.741 1.00 8.78  ? 1   ASN A HB3  1
ATOM   17   H HD21 . ASN A 1 1  ? 12.204 -0.386 27.881 1.00 13.35 ? 1   ASN A HD21 1
ATOM   18   H HD22 . ASN A 1 1  ? 11.681 -1.773 28.031 1.00 13.35 ? 1   ASN A HD22 1
"""

def get_pdb_inp_functions(pdb_inp):
  exclusions = [
      # Seems that nobody uses these 4 outside C++ code, could be private
      'break_indices',
      'chain_indices',
      'model_atom_counts',
      'bookkeeping_section',
      # Not very much used as well
      'extract_LINK_records',
      'extract_authors',
      'extract_connectivity',
      'extract_remark_iii_records',
      'get_link_records',
      'miscellaneous_features_section',
      'primary_structure_section',
      'record_type_counts',
      'secondary_structure_section',
      'sequence_from_SEQRES',
      'unknown_section']
  result = []
  # pdb_inp = iotbx.pdb.input(source_info=None, lines=t_pdb_str.split('\n'))
  pdb_members = dir(pdb_inp)
  for member in pdb_members:
    if not member.startswith('_') and member not in exclusions:
      result.append(member)
  return result

def find_out_extras(list1, list2):
  result = []
  for x in list1:
    if x not in list2:
      result.append(x)
  return result

def exercise():
  fail_mgs = """
  This test ensures similarity of cif_input and pdb_input interfaces.
  If it fails, inspect extra_in_pdb and extra_in_cif to figure out what
  was added to one of the class and not added to another one.
  """
  pdb_inp = iotbx.pdb.input(source_info=None, lines=t_pdb_str.split('\n'))
  cif_inp = iotbx.pdb.input(source_info=None, lines=t_cif_str.split('\n'))

  pdb_member_functions = get_pdb_inp_functions(pdb_inp)
  predicate = inspect.ismethod
  if six.PY3:
    predicate = inspect.isfunction
  cif_member_functions = [
      x[0] for x in inspect.getmembers(cif_input, predicate=predicate) \
      if not x[0].startswith('_')]
  extra_in_pdb = find_out_extras(pdb_member_functions, cif_member_functions)
  extra_in_cif = find_out_extras(cif_member_functions, pdb_member_functions)
  assert extra_in_pdb == [], " %s %s" % (extra_in_pdb, fail_mgs)
  assert extra_in_cif == [], " %s %s" % (extra_in_cif, fail_mgs)
  # check that the signatures are the same
  for cif_member in cif_member_functions:
    pdb_funct = eval('pdb_inp.%s' % cif_member)
    cif_funct = eval('cif_inp.%s' % cif_member)
    try:
      assert inspect.signature(pdb_funct) == inspect.signature(cif_funct), \
          "Function %s has different signatures: %s !=  %s" % (
              cif_member, inspect.signature(pdb_funct), inspect.signature(cif_funct))
    except ValueError:
      pass
  print("OK")
  # print extra_in_pdb
  # print "="*20
  # print extra_in_cif

if (__name__ == "__main__"):
  if sys.version_info.major > 2:
    exercise()
  else:
    print("Skipping - Python 3 required")
