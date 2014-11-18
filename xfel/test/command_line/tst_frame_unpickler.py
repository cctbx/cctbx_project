from __future__ import division
from dials.util.options import Importer, flatten_reflections, flatten_experiments, OptionParser
from xfel.command_line.frame_unpickler import *
from libtbx import easy_pickle
import libtbx.load_env
import os

class test_frame_unpickler(object):
  def __init__(self):
    self.xfel_regression_dir = libtbx.env.find_in_repositories(
      relative_path="xfel_regression",
      test=os.path.isdir)
    self.test_integration_pickles = os.path.join(self.xfel_regression_dir, "merging_test_data", "integration_pickles")
    self.test_experiment_lists = os.path.join(self.xfel_regression_dir, "merging_test_data", "experiment_lists")
    self.test_reflection_tables = os.path.join(self.xfel_regression_dir, "merging_test_data", "reflection_tables")

  def unpickle(self):
    # load contents of integration pickle
    if self.test_integration_pickles is None:
      print "Skipping test_frame_unpickler: test directory integration_pickles not found"
      return

    pickle = os.path.join(self.test_integration_pickles, "int-s00-2011-12-02T21:05Z17.167_00000.pickle")
    if not os.path.exists(pickle):
      print "Skipping test_frame_unpickler: test integration pickle not found"
      return

    self.generated = construct_reflection_table_and_experiment_list(pickle, [pickle], 0.11)

    # write files
    from libtbx.test_utils import open_tmp_directory
    tmp = open_tmp_directory(suffix="test_frame_unpickler")

    self.generated.assemble_experiments()
    self.generated.experiments_to_json(tmp)

    self.generated.assemble_reflections()
    self.generated.reflections_to_pickle(tmp)

    # load back into memory
    refl = os.path.join(tmp, "idx-s00-20111202210517167_integrated.pickle")
    if not os.path.exists(refl):
      print "Reflection table not sucessfully written to pickle file"
      return

    expt = os.path.join(tmp, "idx-s00-20111202210517167_experiments.json")
    if not os.path.exists(expt):
      print "Experiment list not successfully written to json"
      return

    importer = Importer([refl, expt], read_experiments=True, read_reflections=True, check_format=False)
    if importer.unhandled:
      print "Unable to load generated reflection table and experiment list"
      return

    self.gen_expt = flatten_experiments(importer.experiments)[0]
    self.gen_refl = flatten_reflections(importer.reflections)[0]


  def load_comparison_refl_json(self):
    if self.test_experiment_lists is None:
      print "Skipping test_frame_unpickler: test directory experiment_lists not found"
      return

    if self.test_reflection_tables is None:
      print "Skipping test_frame_unpickler: test directory reflection_tables not found"
      return

    refl = os.path.join(self.test_reflection_tables, "idx-s00-20111202210517167_integrated.pickle")
    if not os.path.exists(refl):
      print "Skipping test_frame_unpickler: reflection table not found"
      return

    expt = os.path.join(self.test_experiment_lists, "idx-s00-20111202210517167_experiments.json")
    if not os.path.exists(expt):
      print "Skipping test_frame_unpickler: experiment list not found"
      return

    importer = Importer([refl, expt], read_experiments=True, read_reflections=True, check_format=False)
    if importer.unhandled:
      print "Skipping test_frame_unpickler: unable to load reflection table and experiment list"
      return

    self.compare_expt = flatten_experiments(importer.experiments)[0]
    self.compare_refl = flatten_reflections(importer.reflections)[0]

  def compare(self):

    # beam
    b1, b2 = self.gen_expt.beam, self.compare_expt.beam
    assert b1.get_direction() == b2.get_direction()
    assert b1.get_divergence() == b2.get_divergence()
    assert b1.get_polarization_fraction() == b2.get_polarization_fraction()
    assert b1.get_polarization_normal() == b2.get_polarization_normal()
    assert b1.get_s0() == b2.get_s0()
    assert b1.get_sigma_divergence() == b2.get_sigma_divergence()
    assert b1.get_unit_s0() == b2.get_unit_s0()
    assert b1.get_wavelength() == b2.get_wavelength()

    # crystal
    c1, c2 = self.gen_expt.crystal, self.compare_expt.crystal
    assert c1.get_A() == c2.get_A()
    assert c1.get_B() == c2.get_B()
    assert c1.get_U() == c2.get_U()
    assert c1.get_mosaicity() == c2.get_mosaicity()
    assert c1.get_real_space_vectors() == c2.get_real_space_vectors()
    assert c1.get_space_group().info().symbol_and_number() == c2.get_space_group().info().symbol_and_number()
    assert c1.get_unit_cell().parameters() == c2.get_unit_cell().parameters()
    assert c1.get_unit_cell().fractionalization_matrix() == c2.get_unit_cell().fractionalization_matrix()
    assert c1.get_unit_cell().orthogonalization_matrix() == c2.get_unit_cell().orthogonalization_matrix()

    # detector
    d1, d2 = self.gen_expt.detector, self.compare_expt.detector
    assert d1.to_dict() == d2.to_dict()

    # goniometer
    g1, g2 = self.gen_expt.goniometer, self.compare_expt.goniometer
    assert g1 == g2 == None

    # imageset
    i1, i2 = self.gen_expt.imageset, self.compare_expt.imageset
    assert i1.paths() == i2.paths()
    assert i1.is_valid() == i2.is_valid() == True
    assert i1.indices() == i2.indices()
    assert i1.get_scan() == i2.get_scan() == None
    assert i1.get_goniometer() == i2.get_goniometer() == None
    assert type(i1.reader()) == type(i2.reader())

    # scan
    s1, s2 = self.gen_expt.scan, self.compare_expt.scan
    assert s1 == s2 == None
  
    # reflection table
    r1, r2 = self.gen_refl, self.compare_refl
    assert r1['background.mean'].all_eq(r2['background.mean'])
    assert r1['background.mse'].all_eq(r2['background.mse'])
    assert r1['bbox'].as_int().all_eq(r2['bbox'].as_int())
    assert r1['correlation.ideal.profile'].all_eq(r2['correlation.ideal.profile'])
    assert r1['entering'].all_eq(r2['entering'])
    assert r1['flags'].all_eq(r2['flags'])
    assert r1['id'].all_eq(r2['id'])
    assert r1['intensity.prf.value'].all_eq(r2['intensity.prf.value'])
    assert r1['intensity.prf.variance'].all_eq(r2['intensity.prf.variance'])
    assert r1['intensity.sum.value'].all_eq(r2['intensity.sum.value'])
    assert r1['intensity.sum.variance'].all_eq(r2['intensity.sum.variance'])
    assert r1['miller_index'].all_eq(r2['miller_index'])
    assert r1['n_background'].all_eq(r2['n_background'])
    assert r1['n_foreground'].all_eq(r2['n_foreground'])
    assert r1['panel'].all_eq(r2['panel'])
    assert r1['s1'].as_double().all_eq(r2['s1'].as_double())
    assert r1['xyzcal.mm'].as_double().all_eq(r2['xyzcal.mm'].as_double())
    assert r1['xyzcal.px'].as_double().all_eq(r2['xyzcal.px'].as_double())
    assert r1['xyzobs.px.value'].as_double().all_eq(r2['xyzobs.px.value'].as_double())
    assert r1['xyzobs.px.variance'].as_double().all_eq(r2['xyzobs.px.variance'].as_double())
    assert r1['zeta'].all_eq(r2['zeta'])

    # if no assertions failed, return successful result
    print 'OK'
    return

  def run_all(self):
    self.unpickle()
    self.load_comparison_refl_json()
    self.compare()

if __name__ == '__main__':
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()

  tester = test_frame_unpickler()
  tester.run_all()

