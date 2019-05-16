from __future__ import absolute_import, division, print_function
from mmtbx.refinement.minimization_ncs_constraints import grads_one_ncs_to_asu
from scitbx.array_family import flex
import mmtbx.ncs.ncs_utils as nu
import iotbx.ncs as ncs
import iotbx.pdb
import unittest
from six.moves import zip
from six.moves import map

__author__ = 'Youval'

class TestMinimizationFunctions(unittest.TestCase):
  """ Test NCS constraints minimization function
  exclude_selection have to be used because NAG is not included into NCS
  search procedure by default anymore"""

  def test_grads_one_ncs_to_asu(self):
    # No more NAGs in NCS selection
    # print sys._getframe().f_code.co_name
    pdb_inp = iotbx.pdb.input(lines=test_pdb_1,source_info=None)
    p = ncs.input.get_default_params()
    p.ncs_search.exclude_selection=None
    ncs_inp = ncs.input(
        hierarchy=pdb_inp.construct_hierarchy(),
        params=p.ncs_search)
    pdb_inp = iotbx.pdb.input(source_info=None, lines=test_pdb_1)
    ph = pdb_inp.construct_hierarchy()
    xrs =  pdb_inp.xray_structure_simple()
    #
    nrgl = ncs_inp.get_ncs_restraints_group_list()
    asu_length = ncs_inp.truncated_hierarchy.atoms_size()
    #
    refine_selection = nu.get_refine_selection(number_of_atoms=asu_length)
    extended_ncs_selection = nrgl.get_extended_ncs_selection(
        refine_selection=refine_selection)
    #
    self.assertEqual(asu_length, ph.atoms_size())
    self.assertEqual(asu_length, 18)
    #
    xrs_one_ncs_copy = xrs.select(extended_ncs_selection)
    master_grad = xrs_one_ncs_copy.extract_u_iso_or_u_equiv()
    #
    g = grads_one_ncs_to_asu(
      ncs_restraints_group_list=nrgl,
      total_asu_length=asu_length,
      extended_ncs_selection=extended_ncs_selection,
      master_grad=master_grad)
    #
    self.assertEqual(g.size(),18)
    masters = [[0,1,12],[2, 3, 13]]
    copies = [[[4, 5, 14],[8, 9, 16]],[[6, 7, 15],[10, 11, 17]]]
    for m,cs in zip(masters,copies):
      ml = list(g.select(flex.size_t(m)))
      for c in cs:
        cl = list(g.select(flex.size_t(c)))
        self.assertEqual(ml,cl)

#
# same chain id for HETATM at the end is valid, see
# http://deposit.rcsb.org/format-faq-v1.html
# Q.  How are three-letter residue names and chain identifiers for residue
#    modifications and ligands assigned?
#
# A.  There are four common cases: covalently bound ligands <...>
# Covalently bound ligands:
# Covalently bound ligands are assigned the chain identifier of the
# polymer chain to which the ligand is bound. The bonding between the ligand
# and the residue is specified in PDB LINK records. The ligand coordinates
# appear as HETATM records following either the TER record for the bound chain
# or after the TER record for the last polymer chain. The ligand is assigned
# a unique residue number within chain to which it is bound. The residue that
# binds the ligand retains its standard name in both coordinate and SEQRES
# records. Additional PDB records MODRES/HET/HETNAM/FORMUL and CONECT are
# provided to describe the ligand.


#
test_pdb_1 = '''\
ATOM      1  N   ASP A   1      27.619  71.759 115.947  1.00138.20           N
ATOM      2  CA  ASP A   1      27.294  72.259 114.616  1.00139.50           C
TER
ATOM   2602  N   PHE B   2      20.530  61.391  90.683  1.00101.18           N
ATOM   2603  CA  PHE B   2      19.441  60.481  91.023  1.00108.06           C
TER
ATOM   3904  N   ASP C   1      15.311  44.604 115.942  1.00136.32           N
ATOM   3905  CA  ASP C   1      15.015  44.052 114.626  1.00138.25           C
TER
ATOM   6505  N   PHE D   2      27.586  43.459  90.762  1.00 99.57           N
ATOM   6506  CA  PHE D   2      28.906  42.913  91.061  1.00109.41           C
TER
ATOM   7799  N   ASP E   1      45.077  47.388 115.919  1.00137.67           N
ATOM   7800  CA  ASP E   1      45.625  47.420 114.569  1.00138.21           C
TER
ATOM  10400  N   PHE F   2      39.568  58.805  90.625  1.00104.40           N
ATOM  10401  CA  PHE F   2      39.333  60.197  90.995  1.00112.24           C
TER
HETATM11702  C1  NAG A 401      17.881  76.973  81.221  1.00 53.50           C
HETATM11828  C1  NAG B 201      35.771  81.610 114.567  1.00 92.78           C
HETATM11842  C1  NAG C 401      15.291  33.478  81.084  1.00 54.99           C
HETATM11968  C1  NAG D 201       2.514  47.082 114.568  1.00 93.06           C
HETATM11982  C1  NAG E 401      54.251  53.233  81.004  1.00 53.68           C
HETATM12108  C1  NAG F 201      49.267  35.290 114.417  1.00 93.40           C
END
'''

def run_selected_tests():
  """  Run selected tests

  1) List in "tests" the names of the particular test you want to run
  2) Comment out unittest.main()
  3) Un-comment unittest.TextTestRunner().run(run_selected_tests())
  """
  tests = ['test_split_groups_to_spec']
  suite = unittest.TestSuite(list(map(TestSimpleAlignment,tests)))
  return suite

if __name__=='__main__':
  # use for individual tests
  # unittest.TextTestRunner().run(run_selected_tests())

  # Use to run all tests
  unittest.main(verbosity=0)
