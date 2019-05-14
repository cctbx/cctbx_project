
"""
Regression tests for mmtbx.scaling.absences
"""

from __future__ import absolute_import, division, print_function
from mmtbx.scaling.absences import *
from cctbx.development import random_structure
from cctbx import sgtbx
from scitbx.array_family import flex
import libtbx.load_env
import os.path as op
import random

def exercise_abs_list():
  tmp = absences()
  assert tmp.check( "2_1 (c)", (0,0,1),True ) == (True, False)
  assert tmp.check( "2_1 (c)", (0,0,4),True ) == (True, True)
  assert tmp.check( "4_1 (a)", (4,0,0),True ) == (True, True)
  assert tmp.check( "3_1 (c)", (0,0,3),True ) == (True, True)

  tmp = sgi_iterator(chiral = True,
    crystal_system = None,
    intensity_symmetry = sgtbx.space_group_info(
      "P222").group().build_derived_reflection_intensity_group(False))
  sg_list = []
  abs_list = []
  for sg in tmp.list():
    sg_list.append( str(sg) )
    if str(sg)== "P 21 21 21":
      for s in sg.group():
        abs_list.append( conditions_for_operator( s ).absence_type() )
  assert "2_1 (a)" in abs_list
  assert "2_1 (b)" in abs_list
  assert "2_1 (c)" in abs_list

  assert "P 2 2 2" in sg_list
  assert "P 21 2 2" in sg_list
  assert "P 2 21 2" in sg_list
  assert "P 2 2 21" in sg_list
  assert "P 21 21 2" in sg_list
  assert "P 21 2 21" in sg_list
  assert "P 2 21 21" in sg_list
  assert "P 21 21 21" in sg_list

def exercise_analyze_absences():
  """
  Looping over common space groups in macromolecular structures, generate
  synthetic data converted to the derived intensity group with absent
  reflections present (I=0), and run the analysis to check that the actual
  space group is the best-scoring.
  """
  # space groups found in structures with data in the PDB as of Sep. 2011,
  # minus a few unusual ones
  symbols = ['C 2 2 21', 'P 31 1 2', 'P 62 2 2', 'P 2 3', 'P 3', 'P 1', 'P 63',
    'C 1 2 1', 'F 2 2 2', 'A 1 2 1', 'B 1 1 2', 'P 61', 'P 61 2 2', 'P 6 2 2',
    'P 2 21 21', 'P 21 21 21', 'R 3 2 :R', 'P 4 3 2', 'P 1 21 1', 'P 32 2 1',
    'P 42 2 2', 'P 43 21 2', 'P 21 2 21', 'I 2 3', 'P 43 3 2', 'I 41 3 2',
    'P 3 1 2', 'P 2 2 2', 'R 3 :H', 'P 21 21 2', 'P 65 2 2',
    'P 21 21 2', 'F 2 3', 'I 41 2 2', 'P 65', 'P 1 1 21', 'P 2 2 21',
    'P 4 2 2', 'I 2 2 2', 'I 4 3 2', 'R 3 :R', 'P 4', 'P 42', 'P 64 2 2',
    'I 1 2 1', 'P 64', 'C 1 2 1', 'I 41', 'P 63 2 2', 'P 3 2 1', 'P 41 2 2',
    'P 62', 'P 1', 'P 32 1 2', 'F 4 3 2', 'P 31 2 1', 'I 4 2 2', 'P 21 3',
    'P 43 2 2', 'C 2 2 2', 'P 4 21 2', 'B 2 21 2', 'P 41 21 2', 'P 31', 'P 41',
    'P 41 3 2', 'P 32', 'P 6', 'P 1 2 1', 'I 21 3', 'F 41 3 2',
    'P 42 3 2', 'P 42 21 2', 'P 43', 'I 21 21 21', 'R 3 2 :H', 'I 4']
  random.seed(987654321)
  flex.set_random_seed(987654321) #12345)
  for symbol in  ["P 42 3 2"] : #symbols :
    xrs = random_structure.xray_structure(
      space_group_symbol=symbol,
      elements=['O']*200,
      n_scatterers=200)
    fc = abs(xrs.structure_factors(d_min=1.5).f_calc())
    fc.set_observation_type_xray_amplitude()
    i_calc = fc.f_as_f_sq()
    space_group_info = i_calc.crystal_symmetry().space_group_info()
    n_absent = i_calc.sys_absent_flags().data().count(True)
    assert (n_absent == 0)
    ig = space_group_info.group().build_derived_reflection_intensity_group(
      False)
    complete_set = i_calc.complete_set()
    missing_set = complete_set.lone_set(i_calc)
    assert missing_set.size() == 0 # should be complete so far
    i_calc_orig = i_calc.deep_copy().customized_copy(
      sigmas=flex.double(i_calc.size(), 0.01))
    symm = i_calc.crystal_symmetry().customized_copy(
      space_group_info=ig.info())
    i_calc = i_calc.customized_copy(crystal_symmetry=symm)
    #i_calc.show_summary()
    complete_set = i_calc.complete_set()
    missing_set = complete_set.lone_set(i_calc) # these are the absences
    data = flex.double(missing_set.size(), 0.01)
    sys_abs = missing_set.array(data=data)
    i_calc = i_calc.concatenate(sys_abs)
    sigmas = flex.double(i_calc.size(), 0.01)
    i_calc = i_calc.customized_copy(sigmas=sigmas)
    i_calc = i_calc.set_observation_type_xray_intensity()
    # i_calc is now complete with respect to the derived reflection intensity
    # group, with I/sigma=1 for all systematic absences
    psgc = protein_space_group_choices(miller_array=i_calc,
      original_data=i_calc_orig)
    table = psgc.suggest_likely_candidates()
    min_score = min([ row[-1] for row in table ])
    expected_failures = ["I 4", "I 4 2 2", "P 42 3 2"]
    # The actual space group won't always be the first on the list (if the
    # absence rules are the same for more than one related space group), but
    # it should always have the lowest score.  (I think...)
    for row in table :
      group = row[0]
      score = row[-1]
      n_abs_viol = int(row[4])
      n_pres_viol = int(row[5])
      if (group == str(space_group_info)):
        assert (n_abs_viol == 0), row
        if (n_pres_viol != 0):
          print(group, n_pres_viol)
          psgc.show()
        # XXX Currently this fails for the following space groups:
        #     I 4
        #     I 4 2 2
        #     P 42 3 2
        assert ((score == min_score) or
                (str(space_group_info) in expected_failures))

# XXX depends on phenix_regression
def exercise_2():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/x_001_HD_ricardo_leal.pdb",
    test=op.isfile)

if __name__ == "__main__":
  exercise_abs_list()
  exercise_analyze_absences()
  print("OK")
