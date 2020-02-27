"""Tests for superpose.py"""
from __future__ import absolute_import, division, print_function
import random
import os
import unittest

import iotbx.pdb
import libtbx.load_env
import libtbx.test_utils

import mmtbx.superpose as superpose
from six.moves import range

if (1):
  # Fix random seed to avoid rare failures
  random.seed(0)

def get_filename(filename):
  return libtbx.env.find_in_repositories(relative_path="phenix_regression/pdb/%s"%filename, test=os.path.isfile)

def random_pair(left, right):
  assert left < right
  x = random.randrange(left, right)
  y = random.randrange(left, right)
  if(y < x):
    tmp = x
    x = y
    y = tmp
  return x,y

# Don't test yet.
# class SuperposeTestExamples(unittest.TestCase):
class SuperposeExamples(object):
  def test_example_presets(self):
    # Do a basic Ca fitting (default)
    fixed = superpose.SuperposePDB(filename="1CBS.pdb", desc='fixed', preset='ca')
    moving = superpose.SuperposePDB(filename="1HMT.pdb", desc='moving', preset='ca')
    rmsd, lsq = moving.superpose(fixed)

    # ... or will all atoms
    moving.select_update(preset="all")
    rmsd, lsq = moving.superpose(fixed)

    # ... or with a specific selection string.
    moving.select_update(selection="""pepnames and (name ca) chain A""")
    rmsd, lsq = moving.superpose(fixed)

    # ... or use the Select-o-matic.
    moving.selectomatic(fixed)
    rmsd, lsq = moving.superpose(fixed)

    # Write output.
    # moving.output(lsq, filename="1HMT-fit-1CBS.pdb")

  def test_example_na(self):
    # Fit nucleic acid -- note the default settings will realize it's NA.
    fixed = superpose.SuperposePDB(filename="na1.pdb", desc='fixed')
    moving = superpose.SuperposePDB(filename="na2.pdb", desc='moving')
    rmsd, lsq = moving.superpose(fixed)
    # moving.output(lsq, filename="na-fit.pdb")

  def test_example_multiple(self):
    # Open multiple models for fitting and align chain A Ca's.
    fixed = superpose.SuperposePDB(filename="1OAN.pdb", desc='fixed', preset='ca')
    models = superpose.SuperposePDB.open_models(filename="1S6N.pdb", desc='moving', preset='ca')
    for model in models:
      rmsd, lsq = model.superpose(fixed)

  def test_example_pairwise(self):
    result = {}
    moving = superpose.SuperposePDB(filename="1AON.pdb")
    fixed = superpose.SuperposePDB(filename="1AON.pdb")
    for moving_chain in list(moving.pdb.chains())[:7]:
      for fixed_chain in list(fixed.pdb.chains())[:7]:
        fixed.select_update(preset='ca', chain=fixed_chain.id)
        moving.select_update(preset='ca', chain=moving_chain.id)
        result[(moving_chain.id, fixed_chain.id)] = moving.superpose(fixed)

    keys = set([i[0] for i in sorted(result.keys())])
    print("\t"+"\t".join(keys))
    for key1 in keys:
      rmsds = [result[(key1,key2)][0] for key2 in keys]
      print("%s:\t"%key1, "\t".join(["%0.2f"%i for i in rmsds]))

  def test_example_sieve(self):
    pass
    # Use an alternate fitting method.
    # fixed = SuperposePDBSieve(filename="1CBS.pdb", desc='fixed', preset='ca')
    # moving = superpose.SuperposePDBSieve(filename="1HMT.pdb", desc='moving', preset='ca')
    # rmsd, lsq = moving.superpose(fixed)
    # moving.output(lsq, filename="1HMT-fit-1CBS-sieve.pdb")

class SuperposeTest(unittest.TestCase):
    # OK!
    def test_alignment_used_1(self, filename='fab_a_cut_1.pdb', test='test_alignment_used_1', tolerance=0.003):
      filename = get_filename(filename)
      filename_shifted = "%s.shifted.pdb"%test
      cmd = """phenix.pdbtools %(filename)s suffix=none output.prefix=%(filename_shifted)s remove="%(gap_selection)s" %(shift)s"""%{
          'filename': filename,
          'filename_shifted': filename_shifted.replace(".pdb",""),
          'gap_selection': "chain A and (resid 5:30 or resid 32:50)",
          'shift': 'sites.rotate="10 20 30" sites.translate="10 20 30"'
      }
      libtbx.test_utils.run_command(command=cmd, result_file_names=[filename_shifted])
      moving = superpose.SuperposePDB(filename_shifted)
      target = superpose.SuperposePDB(filename)
      rmsd, lsq = moving.superpose(target)
      self.assertLess(rmsd, tolerance)

    # OK!
    def test_alignment_used_2(self, filename='fab_a_cut.pdb', test='test_alignment_used_2', tolerance=0.003):
      filename = get_filename(filename)
      for attempt in range(1):
        filename_shifted = '%s.%s.shifted.pdb'%(test, attempt)
        gap_selection = "chain A and (resid %s:%s or resid %s:%s or resid %s:%s)"%(tuple(
          random_pair(0, 15)+
          random_pair(16, 35)+
          random_pair(36, 50)))
        cmd = """phenix.pdbtools %(filename)s suffix=none output.prefix=%(filename_shifted)s remove="%(gap_selection)s" %(shift)s"""%{
            'filename': filename,
            'filename_shifted': filename_shifted.replace(".pdb",""),
            'gap_selection': gap_selection,
            'shift': 'sites.rotate="10 20 30" sites.translate="10 20 30"'
        }
        libtbx.test_utils.run_command(command=cmd, result_file_names=[filename_shifted])
        moving = superpose.SuperposePDB(filename_shifted)
        target = superpose.SuperposePDB(filename)
        rmsd, lsq = moving.superpose(target)
        self.assertLess(rmsd, tolerance)
      return

    # TODO: Not quite OK
    def test_comprehensive(self, filename='enk_gm.pdb', test='test_comprehensive'):
      filename = get_filename(filename)
      tests = [
        ["", "None", "None"], # OK
        ["""sites.rotate="1 2 3" sites.translate="4 5 6" """, "None", "None"], # OK
        ["""sites.rotate="1 2 3" sites.translate="4 5 6" """, "chain A or resname HOH", "not (chain A or resname HOH)"], # OK
        ["""sites.rotate="1 2 3" sites.translate="4 5 6" """, "chain A or chain C", "not (chain A or chain C)"], # OK
        # ["""sites.rotate="1 2 3" sites.translate="4 5 6" """, "resname HOH", "not (resname HOH)"], # FAILED
        ["""sites.rotate="1 2 3" sites.translate="4 5 6" """, "name CA", "not (name CA)"] # OK
      ]
      count = 0
      for shift, shift_selection, remove_selection in tests:
        print("========= shift: %s -- shift_selection: %s -- remove_selection: %s"%(shift, shift_selection, remove_selection))
        filename_shifted = '%s.%s.shifted.pdb'%(test, count)
        filename_fitted = '%s.%s.fitted.pdb'%(test, count)
        cmd = """phenix.pdbtools %(filename)s suffix=none output.prefix=%(filename_shifted)s selection="%(shift_selection)s" remove="%(remove_selection)s" %(shift)s"""%{
            'filename': filename,
            'filename_shifted': filename_shifted.replace(".pdb",""),
            'shift': shift,
            'shift_selection': shift_selection,
            'remove_selection': remove_selection
        }
        libtbx.test_utils.run_command(command=cmd)

        moving = superpose.SuperposePDB(filename_shifted, selection=shift_selection)
        target = superpose.SuperposePDB(filename, selection=shift_selection)
        rmsd, lsq = moving.superpose(target)
        moving.output(lsq, filename=filename_fitted)

        xs_shifted = iotbx.pdb.input(file_name=filename_shifted).xray_structure_simple()
        xs_fitted = iotbx.pdb.input(file_name=filename_fitted).xray_structure_simple()
        check = xs_fitted.mean_distance(other=xs_shifted)
        if shift:
          self.assertGreater(check, 5.0)
        else:
          self.assertAlmostEqual(check, 0)
        count += 1
      return

    # OK!
    def test_assert_anisou_warning_message(self, filename='phe_a.pdb', test='test_assert_anisou_warning_message'):
      filename = get_filename(filename)
      filename_output = '%s.noaniso.pdb'%test
      moving = superpose.SuperposePDB(filename)
      target = superpose.SuperposePDB(filename)
      rmsd, lsq = moving.superpose(target)
      moving.output(lsq, filename=filename_output)
      pi = iotbx.pdb.input(file_name=filename_output)
      ph = pi.construct_hierarchy()
      pa = ph.atoms()
      uij = pa.extract_uij().as_double()
      print(uij.all_eq(-1.0))
      assert uij.all_eq(-1.0)

    # Mostly OK!
    def test_erik_vogan_case(self, test='test_erik_vogan_case'):
      filename1 = get_filename('rebuild03_ckpt04.pdb')
      filename2 = get_filename('2C30.pdb')
      filename_output = '%s.fitted.pdb'%test

      moving = superpose.SuperposePDB(filename2)
      target = superpose.SuperposePDB(filename1)
      rmsd, lsq = moving.superpose(target)
      moving.output(lsq, filename=filename_output)

      # ["Number of atoms for LS fitting =  36",
      # "RMSD (all matching atoms) (start): 60.830 (number of atoms: 99)",
      # "RMSD (all matching atoms) (final): 0.992 (number of atoms: 99)"]:
      self.assertLess(rmsd, 1.0)

    # Mostly OK!
    def test_altlocs_1(self, test='test_altlocs_1', tolerance=0.001):
      filename1 = get_filename("fab_1.pdb")
      filename2 = get_filename("fab_2.pdb")
      filename_output = '%s.fitted.pdb'%test

      moving = superpose.SuperposePDB(filename2)
      target = superpose.SuperposePDB(filename1)
      rmsd, lsq = moving.superpose(target)
      moving.output(lsq, filename=filename_output)

      # ["Number of atoms for LS fitting =  6",
      # "RMSD between fixed and moving atoms (start): 12.830",
      # "RMSD between fixed and moving atoms (final): 0.000",
      # "RMSD (all matching atoms) (start): 13.164 (number of atoms: 15)",
      # "RMSD (all matching atoms) (final): 0.001 (number of atoms: 15)"]:
      self.assertLess(rmsd, tolerance)

    # OK!
    def test_same_models_many_differently_ordered_chains_1(self, test='test_same_models_many_differently_ordered_chains_1', tolerance=0.03):
      filename1 = get_filename("a1.pdb")
      filename2 = get_filename("b1.pdb")
      filename_output = '%s.fitted.pdb'%test

      moving = superpose.SuperposePDB(filename2)
      target = superpose.SuperposePDB(filename1)
      rmsd, lsq = moving.superpose(target)
      rmsd_output = moving.output(lsq, filename=filename_output)

      # ["Number of atoms for LS fitting =  33",
      # "RMSD between fixed and moving atoms (start): 153.279",
      # "RMSD between fixed and moving atoms (final): 0.025"]:
      self.assertLess(rmsd, tolerance)

    # OK!
    def test_same_models_many_differently_ordered_chains_2(self, test='test_same_models_many_differently_ordered_chains_2', tolerance=0.13):
      filename1 = get_filename("a2.pdb")
      filename2 = get_filename("b2.pdb")
      filename_output = '%s.fitted.pdb'%test

      moving = superpose.SuperposePDB(filename2)
      target = superpose.SuperposePDB(filename1)
      rmsd, lsq = moving.superpose(target)
      rmsd_output = moving.output(lsq, filename=filename_output)

      # ["Number of atoms for LS fitting =  1284",
      # "RMSD (all matching atoms) (start): 159.830 (number of atoms: 26184)",
      # "RMSD (all matching atoms) (final): 67.995 (number of atoms: 26184)"]
      self.assertLess(rmsd, tolerance)

    # OK!
    def test_rna_dna(self, test='test_rna_dna', tolerance=0.7):
      filename1 = get_filename('superpose_pdbs_a.pdb')
      filename2 = get_filename('superpose_pdbs_b.pdb')
      filename_output = '%s.fitted.pdb'%test

      moving = superpose.SuperposePDB(filename2)
      target = superpose.SuperposePDB(filename1)
      rmsd, lsq = moving.superpose(target)
      rmsd_output = moving.output(lsq, filename=filename_output)

      # ["Number of atoms for LS fitting =  8142",
      # "RMSD (all matching atoms) (start): 405.009 (number of atoms: 8142)",
      # "RMSD (all matching atoms) (final): 0.607 (number of atoms: 8142)"]:
      self.assertLess(rmsd, tolerance)

if __name__ == '__main__':
    unittest.main(verbosity=0)
