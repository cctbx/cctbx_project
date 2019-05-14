from __future__ import absolute_import, division, print_function
import iotbx.pdb
from iotbx.pdb.mmcif import cif_input
import inspect

t_pdb_str = """\
ATOM      1  N   ALA     2       1.643  -2.366  -1.408  1.00
ATOM      3  CA  ALA     2       1.280  -3.608  -2.069  1.00
ATOM      6  CB  ALA     2       1.361  -4.762  -1.068  1.00
ATOM     10  C   ALA     2      -0.114  -3.466  -2.684  1.00
ATOM     11  O   ALA     2      -0.327  -3.827  -3.840  1.00
"""

def get_pdb_inp_functions():
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
  pdb_inp = iotbx.pdb.input(source_info=None, lines=t_pdb_str.split('\n'))
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
  pdb_member_functions = get_pdb_inp_functions()
  cif_member_functions = [
      x[0] for x in inspect.getmembers(cif_input, predicate=inspect.ismethod) \
      if not x[0].startswith('_')]
  extra_in_pdb = find_out_extras(pdb_member_functions, cif_member_functions)
  extra_in_cif = find_out_extras(cif_member_functions, pdb_member_functions)
  assert extra_in_pdb == [], " %s %s" % (extra_in_pdb, fail_mgs)
  assert extra_in_cif == [], " %s %s" % (extra_in_cif, fail_mgs)
  print("OK")
  # print extra_in_pdb
  # print "="*20
  # print extra_in_cif

if (__name__ == "__main__"):
  exercise()
