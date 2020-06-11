from __future__ import absolute_import, division, print_function
import mmtbx.model
import iotbx.pdb
from cctbx.array_family import flex
from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal, show_diff
from iotbx.regression.ncs.tst_ncs import pdb_str_5
from mmtbx.regression.ncs.tst_ncs_search import test_pdb_6


def exercise_set_sites_cart_ncs():
  """
  No extra atoms
  """
  inp = iotbx.pdb.input(lines=pdb_str_5, source_info=None)
  model = mmtbx.model.manager(model_input=inp)
  model.search_for_ncs()
  nrgl = model.get_ncs_groups()
  nrgl._show()
  print ('n model atoms:', model.get_number_of_atoms())
  print ('n master atoms:', model.get_master_selection().count(True))
  print ('n master atoms:', model.get_master_selection().iselection().size())
  print ('n master atoms:', model.get_master_hierarchy().atoms_size())

  assert model.get_master_selection().count(True) ==\
      model.get_master_selection().iselection().size() ==\
      model.get_master_hierarchy().atoms_size()

  h = model.get_hierarchy()
  # print('h sites:', list(h.atoms().extract_xyz()))

  # Warning: here here mh is not deep-copy, therefore when we change atom coords
  # they are changing in model.get_hierarchy() as well
  mh = model.get_master_hierarchy()
  new_sites_cart = flex.vec3_double([(1.0, 1.0, 1.0)]*42)
  mh.atoms().set_xyz(new_sites_cart)
  model.set_sites_cart_from_hierarchy(multiply_ncs=True)
  h = model.get_hierarchy()
  new_xyz=h.atoms().extract_xyz()
  # print('h sites:', list(new_xyz))

  # checking if setting went as supposed:
  assert approx_equal(new_xyz[0], (1.0, 1.0, 1.0), eps=1e-4)
  assert approx_equal(new_xyz[42], (-0.6420293330506949, 1.2600792663765976, 7.999997229451341), eps=1e-4)
  assert approx_equal(new_xyz[84], (-1.396802536808536, -0.22123285934616377, 1.000000038694047), eps=1e-4)
  for i in range(3):
    for j in range(42):
      assert approx_equal(new_xyz[42*i+j], new_xyz[42*i], eps=1e-4)

def exercise_set_sites_cart_ncs_with_extra_atoms():
  inp = iotbx.pdb.input(lines=test_pdb_6, source_info=None)
  model = mmtbx.model.manager(model_input=inp)

  model.search_for_ncs()
  nrgl = model.get_ncs_groups()
  nrgl._show(brief=False)
  print ('n model atoms:', model.get_number_of_atoms())
  print ('n master atoms:', model.get_master_selection().count(True))
  print ('n master atoms:', model.get_master_selection().iselection().size())
  print ('n master isel:', list(model.get_master_selection().iselection()))
  print ('n master atoms:', model.get_master_hierarchy().atoms_size())
  print ('n master hierarchy:\n', model.get_master_hierarchy().as_pdb_string())


  # Note that "HETATM   32  C2  NDG H" and "HETATM   35  C2  NDG L"
  # don't belong to any master/copy
  assert not show_diff(model.get_master_hierarchy().as_pdb_string(),"""\
ATOM      1  N   ASP H   5      91.286 -31.834  73.572  1.00 77.83           N
ATOM      2  CA  ASP H   5      90.511 -32.072  72.317  1.00 78.04           C
ATOM      3  C   ASP H   5      90.136 -30.762  71.617  1.00 77.70           C
ATOM      4  O   ASP H   5      89.553 -29.857  72.225  1.00 77.56           O
ATOM      5  N   THR H   6      91.286 -31.834  73.572  1.00 77.83           N
ATOM      6  CA  THR H   6      90.511 -32.072  72.317  1.00 78.04           C
TER
ATOM      7  N   GLY I 501      91.286 -31.834  73.572  1.00 77.83           N
ATOM      8  CA  GLY I 501      90.511 -32.072  72.317  1.00 78.04           C
ATOM      9  C   GLY I 501      90.136 -30.762  71.617  1.00 77.70           C
ATOM     10  O   GLY I 501      89.553 -29.857  72.225  1.00 77.56           O
TER
HETATM   31  C1  NDG H 640      91.286 -31.834  73.572  1.00 77.83           C
HETATM   32  C2  NDG H 640      91.286 -31.834  73.572  1.00 77.83           C
HETATM   35  C2  NDG L 646      61.028 -14.273  81.262  1.00 69.80           C
""")

  mh = model.get_master_hierarchy()
  # Note that atoms outside NCS are getting 2.0 as xyz
  new_sites_cart = flex.vec3_double([
      (1.0, 1.0, 1.0),
      (1.0, 1.0, 1.0),
      (1.0, 1.0, 1.0),
      (1.0, 1.0, 1.0),
      (1.0, 1.0, 1.0),
      (1.0, 1.0, 1.0)]
      + [(3.0, 3.0, 3.0)]*4 +
      [(1.0, 1.0, 1.0),    # <--- Note this atom belongs the first NCS group
      (2.0, 2.0, 2.0),
      (2.0, 2.0, 2.0)]
      )
  mh.atoms().set_xyz(new_sites_cart)
  # print('='*80)
  # print (mh.as_pdb_string())
  model.set_sites_cart_from_hierarchy(multiply_ncs=True)
  h = model.get_hierarchy()
  new_xyz=h.atoms().extract_xyz()
  # print(model.model_as_pdb())
  # print (list(new_xyz))

  assert approx_equal(new_xyz[31], (2.0, 2.0, 2.0), eps=1e-4)
  assert approx_equal(new_xyz[34], (2.0, 2.0, 2.0), eps=1e-4)

  assert approx_equal(new_xyz[0], (1.0, 1.0, 1.0), eps=1e-4)
  assert approx_equal(new_xyz[5], (1.0, 1.0, 1.0), eps=1e-4)
  assert approx_equal(new_xyz[6], (3.0, 3.0, 3.0), eps=1e-4)
  assert approx_equal(new_xyz[9], (3.0, 3.0, 3.0), eps=1e-4)


  assert nrgl.check_for_max_rmsd(sites_cart=new_xyz, chain_max_rmsd=0.0)

  for i in [[0, 1, 2, 3, 4, 5, 30], [10, 11, 12, 13, 14, 15, 32],
      [20, 21, 22, 23, 24, 25, 33]]:
    for j in i:
      assert approx_equal(new_xyz[j], new_xyz[i[0]], eps=1e-4)
  for i in [[6, 7, 8, 9], [16, 17, 18, 19], [26, 27, 28, 29]]:
    for j in i:
      assert approx_equal(new_xyz[j], new_xyz[i[0]], eps=1e-4)

def exercise_set_sites_cart_no_ncs():
  inp = iotbx.pdb.input(lines=pdb_str_5, source_info=None)
  model = mmtbx.model.manager(model_input=inp, expand_with_mtrix=False)
  model.search_for_ncs()
  nrgl = model.get_ncs_groups()
  nrgl._show()
  print ('n model atoms:', model.get_number_of_atoms())
  print ('n master atoms:', model.get_master_selection().count(True))
  print ('n master atoms:', model.get_master_selection().iselection().size())
  print ('n master atoms:', model.get_master_hierarchy().atoms_size())

  assert model.get_master_selection().count(True) ==\
      model.get_master_selection().iselection().size() ==\
      model.get_master_hierarchy().atoms_size()

  h = model.get_hierarchy()
  # Warning: here here mh is not deep-copy, therefore when we change atom coords
  # they are changing in model.get_hierarchy() as well
  mh = model.get_master_hierarchy()
  new_sites_cart = flex.vec3_double([(1.0, 1.0, 1.0)]*42)
  mh.atoms().set_xyz(new_sites_cart)
  model.set_sites_cart_from_hierarchy(multiply_ncs=True)
  h = model.get_hierarchy()
  new_xyz=h.atoms().extract_xyz()
  # print('h sites:', list(new_xyz))
  # checking if setting went as supposed:
  for j in range(42):
    assert approx_equal(new_xyz[j], (1.0, 1.0, 1.0), eps=1e-4)

def run():
  exercise_set_sites_cart_ncs()
  exercise_set_sites_cart_ncs_with_extra_atoms()
  exercise_set_sites_cart_no_ncs()
  print(format_cpu_times())

if (__name__ == "__main__"):
  run()
