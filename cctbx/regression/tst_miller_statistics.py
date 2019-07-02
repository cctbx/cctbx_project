
from __future__ import absolute_import, division, print_function
from iotbx import file_reader
from cctbx.array_family import flex
from cctbx import crystal
from cctbx import xray
import libtbx.load_env
from libtbx.test_utils import approx_equal
from six.moves import cStringIO as StringIO
import os.path

def exercise():
  hkl_file = libtbx.env.find_in_repositories(
    relative_path="phenix_examples/p9-sad/p9.sca",
    test=os.path.isfile)
  if (hkl_file is None):
    print("phenix_examples/p9-sad/p9.sca not available, skipping test")
    return
  # XXX very hard to avoid cross-imports here if we want to test on actual
  # real data...
  hkl_in = file_reader.any_file(hkl_file)
  i_obs = hkl_in.file_server.miller_arrays[0]
  i_obs.setup_binner(n_bins=50)
  meas = i_obs.measurability()
  assert approx_equal(meas, 0.091044)
  meas_binned = i_obs.measurability(use_binning=True)
  bijvoet_ratios = i_obs.bijvoet_ratios(measurable_only=True)
  assert (bijvoet_ratios.size() == 1713)
  assert approx_equal(flex.mean(bijvoet_ratios), 0.476718)
  bijvoet_ratios = i_obs.bijvoet_ratios(measurable_only=False)
  assert (bijvoet_ratios.size() == 18815)
  assert approx_equal(flex.mean(bijvoet_ratios), 0.2727699)
  wilson_mean = i_obs.wilson_plot()
  assert approx_equal(wilson_mean, 18481.1252)
  i_obs.setup_binner(n_bins=10)
  wilson_binned = i_obs.wilson_plot(use_binning=True)
  out = StringIO()
  wilson_binned.show(f=out)
  assert (out.getvalue() == """\
unused:         - 28.4910 [   0/7   ]
bin  1: 28.4910 -  3.7555 [4161/4167] 79994.6693
bin  2:  3.7555 -  2.9819 [4160/4164] 47274.9530
bin  3:  2.9819 -  2.6052 [4175/4177] 17682.1962
bin  4:  2.6052 -  2.3672 [4133/4139] 10444.3793
bin  5:  2.3672 -  2.1976 [4115/4133] 9079.7295
bin  6:  2.1976 -  2.0680 [4180/4210] 6529.7430
bin  7:  2.0680 -  1.9645 [4076/4123] 5050.6608
bin  8:  1.9645 -  1.8790 [4109/4156] 3617.7585
bin  9:  1.8790 -  1.8067 [4115/4157] 2096.0724
bin 10:  1.8067 -  1.7443 [3957/4151] 1472.1231
unused:  1.7443 -         [   0/0   ]
""")
  # TODO lots more stuff - eventually most of the functions used in Xtriage
  # should be tested here

def exercise_unmerged():
  quartz_structure = xray.structure(
    special_position_settings=crystal.special_position_settings(
      crystal_symmetry=crystal.symmetry(
        unit_cell=(5.01,5.01,5.47,90,90,120),
        space_group_symbol="P6222")),
    scatterers=flex.xray_scatterer([
      xray.scatterer(
        label="Si",
        site=(1/2.,1/2.,1/3.),
        u=0.2),
      xray.scatterer(
        label="O",
        site=(0.197,-0.197,0.83333),
        u=0.1)]))
  quartz_structure.set_inelastic_form_factors(
    photon=1.54,
    table="sasaki")
  fc = abs(quartz_structure.structure_factors(d_min=1.0).f_calc())
  symm = fc.crystal_symmetry()
  icalc = fc.expand_to_p1().f_as_f_sq().set_observation_type_xray_intensity()
  # generate 'unmerged' data
  i_obs = icalc.customized_copy(crystal_symmetry=symm)
  # now make up sigmas and some (hopefully realistic) error
  flex.set_random_seed(12345)
  n_refl = i_obs.size()
  sigmas = flex.random_double(n_refl) * flex.mean(fc.data())
  sigmas = icalc.customized_copy(data=sigmas).apply_debye_waller_factors(
    u_iso=0.15)
  err = (flex.double(n_refl, 0.5) - flex.random_double(n_refl)) * 2
  i_obs = i_obs.customized_copy(
    sigmas=sigmas.data(),
    data=i_obs.data() + err)
  # check for unmerged acentrics
  assert i_obs.is_unmerged_intensity_array()
  i_obs_centric = i_obs.select(i_obs.centric_flags().data())
  i_obs_acentric = i_obs.select(~(i_obs.centric_flags().data()))
  i_mrg_acentric = i_obs_acentric.merge_equivalents().array()
  i_mixed = i_mrg_acentric.concatenate(i_obs_centric)
  assert not i_mixed.is_unmerged_intensity_array()
  # XXX These results of these functions are heavily dependent on the
  # behavior of the random number generator, which is not consistent across
  # platforms - therefore we can only check for very approximate values.
  # Exact numerical results are tested with real data (stored elsewhere).
  # CC1/2, etc.
  assert approx_equal(i_obs.cc_one_half(), 0.9999, eps=0.001)
  assert approx_equal(i_obs.cc_one_half_sigma_tau(), 0.9999, eps=0.001)
  assert i_obs.resolution_filter(d_max=1.2).cc_one_half() > 0
  assert i_obs.cc_anom() > 0.1
  r_ano = i_obs.r_anom()
  assert approx_equal(r_ano, 0.080756, eps=0.0001)
  # merging stats
  i_mrg = i_obs.merge_equivalents()
  assert i_mrg.r_merge() < 0.1
  assert i_mrg.r_meas() < 0.1
  assert i_mrg.r_pim() < 0.05

if (__name__ == "__main__"):
  exercise()
  exercise_unmerged()
  print("OK")
