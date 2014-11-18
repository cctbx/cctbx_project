from __future__ import division
import libtbx.load_env
import os

class test_cxi_index(object):
  def __init__(self):

    self.dials_regression = libtbx.env.find_in_repositories(
        relative_path="dials_regression",
        test=os.path.isdir)

    self.labelit_regression = libtbx.env.find_in_repositories(
          relative_path="labelit_regression",
          test=os.path.isdir)

  def test_run_one_index(self):
    if self.dials_regression is None:
      print "Skipping test_run_one_index: dials_regression not present"
      return

    if self.labelit_regression is None:
      print "Skipping test_run_one_index: labelit_regression not present"
      return

    image_pickle = os.path.join(self.dials_regression, "image_examples", "LCLS_CXI",
                                "shot-s04-20111204004533388.pickle")

    if not os.path.exists(image_pickle):
      print "Skipping test_run_one_index: image pickle %s not present"%image_pickle
      return

    input_phil = os.path.join(self.labelit_regression, "xfel", "L498-thermolysin-27.phil")
    if not os.path.exists(input_phil):
      print "Skipping test_run_one_index: input phil %s not present"%input_phil
      return

    from libtbx.test_utils import open_tmp_directory
    tmp_dir = open_tmp_directory(suffix="test_cxi_index")

    int_pickle_path = os.path.join(tmp_dir, "tmp_int.pickle")
    print "Integration result will be found at", int_pickle_path

    from xfel.cxi.display_spots import run_one_index
    result = run_one_index(image_pickle, *["target=%s"%input_phil,
                                           "indexing.completeness_pickle=%s"%int_pickle_path,
                                           '--nodisplay'], **({'display':False}))

    # test the output
    from libtbx import easy_pickle
    from libtbx.test_utils import approx_equal
    data = easy_pickle.load(int_pickle_path)

    assert len(data) == 18
    required_keys = ['mapped_predictions','distance','ybeam','current_orientation','ML_half_mosaicity_deg','current_cb_op_to_primitive',
                     'effective_tiling','residual','sa_parameters','model_partialities','ML_domain_size_ang','mosaicity','observations',
                     'wavelength','xbeam','pointgroup','max_signal','correction_vectors']
    for key in required_keys:
      assert key in data

    obs = data['observations'][0]
    unit_cell = obs.unit_cell()

    expected_parameters = 91.81856980895002, 91.81856980895002, 130.56616061328216, 90.0, 90.0, 120.0
    derived_parameters = unit_cell.parameters()
    print "Expected unit cell parameters", expected_parameters
    print "Derived unit cell parameters", derived_parameters

    epsilon = 1.5
    for p_expected, p_derived in zip(expected_parameters, derived_parameters):
      assert approx_equal(p_expected, p_derived, epsilon)

    print "OK"

  def run_all(self):
    self.test_run_one_index()


if __name__ == '__main__':
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()

  tester = test_cxi_index()
  tester.run_all()

