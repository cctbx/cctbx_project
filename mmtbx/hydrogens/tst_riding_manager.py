from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from libtbx.utils import null_out
from libtbx import group_args

#-------------------------------------------------
# Test for deep_copy method in riding H manager
#-------------------------------------------------
def exercise1():
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  model.process(make_restraints=True)
  pdb_hierarchy = model.get_hierarchy()
  sites_cart = model.get_sites_cart()
  atoms = pdb_hierarchy.atoms()

  model.setup_riding_h_manager()
  riding_h_manager = model.get_riding_h_manager()

  h_parameterization = riding_h_manager.h_parameterization
  number_h = model.get_hd_selection().count(True)
  number_h_para = len(h_parameterization) - h_parameterization.count(None)

  assert (number_h_para == number_h)

  rc = h_parameterization[20]

  new_manager = riding_h_manager.deep_copy()
  new_h_parameterization = new_manager.h_parameterization
  rc_new = new_h_parameterization[20]

  assert (rc_new.htype == rc.htype)
  assert (rc_new.ih    == rc.ih)
  assert (rc_new.a0    == rc.a0)
  assert (rc_new.a1    == rc.a1)
  assert (rc_new.a2    == rc.a2)
  assert (rc_new.a3    == rc.a3)
  assert (rc_new.a     == rc.a)
  assert (rc_new.b     == rc.b)
  assert (rc_new.h     == rc.h)
  assert (rc_new.n     == rc.n)
  assert (rc_new.disth == rc.disth)

  rc_new.htype = 'unk'
  rc_new.ih    = 0
  rc_new.a0    = 1
  rc_new.a1    = 2
  rc_new.a2    = 3
  rc_new.a3    = 4

  assert (rc_new.htype != rc.htype)
  assert (rc_new.ih != rc.ih)
  assert (rc_new.a0 != rc.a0)
  assert (rc_new.a1 != rc.a1)
  assert (rc_new.a2 != rc.a2)
  assert (rc_new.a3 != rc.a3)


def apply_selection(riding_h_manager, selection, answers):
  entries_original_h_para = 21
  entries_cpp_original    = 9
  h_parameterization = riding_h_manager.h_parameterization

  pdb_hierarchy = riding_h_manager.pdb_hierarchy
  n_atoms = pdb_hierarchy.atoms_size()

  selected_manager = riding_h_manager.select(selection)
  new_h_parameterization = selected_manager.h_parameterization

  # Check that length of h_para corresponds to known answer and to number
  # of atoms in hierarchy
  assert (len(h_parameterization) == entries_original_h_para)
  assert (len(h_parameterization) == n_atoms)
  assert (len(new_h_parameterization) == answers.entries_selected_h_para)
  assert (
    len(new_h_parameterization) == selected_manager.pdb_hierarchy.atoms_size())
  assert (len(riding_h_manager.parameterization_cpp) == entries_cpp_original)
  #print('number of entries in new_h_para_selected', len(new_h_parameterization))
  #print('para c++ new', len(selected_manager.parameterization_cpp))
  assert (
    len(selected_manager.parameterization_cpp) == answers.entries_cpp_selected)


  assert (
    selected_manager.hd_selection.count(True) == answers.h_in_selected_hierarchy)

  # check if calculated H position is close to input position (input H is ideal)
  # This will (in most cases) produce large deviations if reindexing failed
  atoms = selected_manager.pdb_hierarchy.atoms()
  sites_cart = atoms.extract_xyz()
  diagnostics = selected_manager.diagnostics(
    sites_cart         = sites_cart,
    threshold          = 0.05)
  h_distances = diagnostics.h_distances
  for ih in h_distances:
    assert (h_distances[ih] < 0.0001), \
      'distance too large: %s  atom: %s distance: %s ' \
      % (new_h_parameterization[ih].htype, atoms[ih].id_str(), h_distances[ih])

#-------------------------------------------------
# Test for select method in riding H manager
#-------------------------------------------------
def exercise2():
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  model.process(make_restraints=True)
  pdb_hierarchy = model.get_hierarchy()

  model.setup_riding_h_manager()
  riding_h_manager = model.get_riding_h_manager()

  h_parameterization = riding_h_manager.h_parameterization
  number_h = model.get_hd_selection().count(True)
  number_h_para = len(h_parameterization) - h_parameterization.count(None)

  # Check that all H atoms are parameterized in original manager
  assert (number_h_para == number_h)

  # Test several selections and compare with known answers
  #
  #print('selection1 = not (name CE1 or name CB)')
  selection1 = pdb_hierarchy.atom_selection_cache().\
      selection("not (name CE1 or name CB)")
  answers1 = group_args(
    entries_selected_h_para = 19,
    entries_cpp_selected    = 3,
    h_in_selected_hierarchy = 9)
  apply_selection(
    riding_h_manager = riding_h_manager,
    selection        = selection1,
    answers          = answers1)
  #
  #print('selection2 = not (name CD1 or name HD2 or name C)')
  selection2 = pdb_hierarchy.atom_selection_cache().\
      selection("not (name CD1 or name HD2 or name C)")
  answers2 = group_args(
    entries_selected_h_para = 18,
    entries_cpp_selected    = 4,
    h_in_selected_hierarchy = 8)
  apply_selection(
    riding_h_manager = riding_h_manager,
    selection        = selection2,
    answers          = answers2)
  #
  #print('selection3 = not (name N or name HB2 or name CE2)')
  selection3 = pdb_hierarchy.atom_selection_cache().\
      selection("not (name N or name HB2 or name CE2)")
  answers3 = group_args(
    entries_selected_h_para = 18,
    entries_cpp_selected    = 4,
    h_in_selected_hierarchy = 8)
  apply_selection(
    riding_h_manager = riding_h_manager,
    selection        = selection3,
    answers          = answers3)

# Ideal amino acid TYR
pdb_str = """\
CRYST1   17.392   14.073   15.364  90.00  90.00  90.00 P 1
SCALE1      0.057498  0.000000  0.000000        0.00000
SCALE2      0.000000  0.071058  0.000000        0.00000
SCALE3      0.000000  0.000000  0.065087        0.00000
ATOM      1  N   TYR D   8      10.136   5.241   8.802  1.00  0.00           N
ATOM      2  CA  TYR D   8      10.978   6.381   8.437  1.00  0.00           C
ATOM      3  C   TYR D   8      12.298   6.240   9.189  1.00  0.00           C
ATOM      4  O   TYR D   8      12.392   6.606  10.364  1.00  0.00           O
ATOM      5  CB  TYR D   8      10.292   7.703   8.763  1.00  0.00           C
ATOM      6  CG  TYR D   8       9.066   7.978   7.922  1.00  0.00           C
ATOM      7  CD1 TYR D   8       9.182   8.482   6.634  1.00  0.00           C
ATOM      8  CD2 TYR D   8       7.791   7.735   8.417  1.00  0.00           C
ATOM      9  CE1 TYR D   8       8.065   8.735   5.861  1.00  0.00           C
ATOM     10  CE2 TYR D   8       6.668   7.985   7.653  1.00  0.00           C
ATOM     11  CZ  TYR D   8       6.810   8.485   6.376  1.00  0.00           C
ATOM     12  OH  TYR D   8       5.694   8.736   5.611  1.00  0.00           O
ATOM     13  H   TYR D   8      10.218   5.000   9.623  1.00  0.00           H
ATOM     14  HA  TYR D   8      11.163   6.357   7.486  1.00  0.00           H
ATOM     15  HB2 TYR D   8      10.017   7.691   9.693  1.00  0.00           H
ATOM     16  HB3 TYR D   8      10.921   8.426   8.615  1.00  0.00           H
ATOM     17  HD1 TYR D   8      10.027   8.652   6.284  1.00  0.00           H
ATOM     18  HD2 TYR D   8       7.693   7.397   9.278  1.00  0.00           H
ATOM     19  HE1 TYR D   8       8.158   9.073   5.000  1.00  0.00           H
ATOM     20  HE2 TYR D   8       5.821   7.817   7.997  1.00  0.00           H
ATOM     21  HH  TYR D   8       5.000   8.541   6.042  1.00  0.00           H
TER
END
"""

if (__name__ == "__main__"):
  t0 = time.time()
  exercise1()
  exercise2()
  print("OK. Time: %8.3f"%(time.time()-t0))
