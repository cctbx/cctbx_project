from cctbx import miller
from cctbx import maptbx
from cctbx import sgtbx
from cctbx import xray
from cctbx import euclidean_model_matching as emma
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.development import structure_factor_utils
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal, show_diff
import pickle
from cStringIO import StringIO
import random
import sys

def check_peaks(structure, peak_sites, max_min_dist):
  for scatterer in structure.scatterers():
    site_symmetry = structure.site_symmetry(scatterer.site)
    equiv_sites = sgtbx.sym_equiv_sites(site_symmetry)
    min_dist = None
    for peak_site in peak_sites:
      dist_info = sgtbx.min_sym_equiv_distance_info(equiv_sites, peak_site)
      if (min_dist is None):
        min_dist = dist_info.dist()
      else:
        min_dist = min(min_dist, dist_info.dist())
    assert min_dist <= max_min_dist, (min_dist, max_min_dist)

def run_test(space_group_info, n_elements=5, d_min=1.5,
             grid_resolution_factor=1./3, max_prime=5, verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["Si"]*n_elements,
    volume_per_atom=200,
    min_distance=3.,
    general_positions_only=False)
  miller_set_f_obs = miller.build_set(
    crystal_symmetry=structure,
    anomalous_flag=(random.random()>0.5),
    d_min=d_min)
  f_obs = miller_set_f_obs.structure_factors_from_scatterers(
    xray_structure=structure,
    algorithm="direct").f_calc()
  structure_factor_utils.check_phase_restrictions(f_obs, verbose=verbose)
  if (0 or verbose):
    f_obs.show_summary()
  if (0 or verbose):
    f_obs.show_array()
  fft_map = f_obs.fft_map(
    resolution_factor=grid_resolution_factor,
    symmetry_flags=maptbx.use_space_group_symmetry)
  p = pickle.dumps(fft_map)
  l = pickle.loads(p)
  s1 = StringIO()
  fft_map.statistics().show_summary(f=s1)
  s2 = StringIO()
  l.statistics().show_summary(f=s2)
  assert not show_diff(s2.getvalue(), s1.getvalue())
  #
  if (not f_obs.anomalous_flag()):
    maptbx_fft_map = maptbx.fft_to_real_map_unpadded(
      space_group=fft_map.space_group(),
      n_real=fft_map.n_real(),
      miller_indices=f_obs.indices(),
      data=f_obs.data())
    fft_map_unpadded = fft_map.real_map_unpadded()
    assert approx_equal(
      flex.linear_correlation(
        fft_map_unpadded.as_1d(), maptbx_fft_map.as_1d()).coefficient(), 1)
    assert approx_equal(
      flex.max(flex.abs(maptbx_fft_map - fft_map_unpadded)), 0)
  #
  fft_map.apply_sigma_scaling()
  real_map = maptbx.copy(
    fft_map.real_map(),
    flex.grid(fft_map.real_map().focus()))
  grid_tags = maptbx.grid_tags(real_map.focus())
  grid_tags.build(
    fft_map.space_group_info().type(),
    fft_map.symmetry_flags())
  assert grid_tags.n_grid_misses() == 0
  assert grid_tags.verify(real_map)
  rms = []
  for interpolate in (False,True):
    peak_list = maptbx.peak_list(
      data=real_map,
      tags=grid_tags.tag_array(),
      peak_search_level=1,
      max_peaks=2*n_elements,
      interpolate=interpolate)
    assert peak_list.gridding() == real_map.focus()
    check_peaks(structure, peak_list.sites(), d_min * grid_resolution_factor)
    crystal_gridding_tags = fft_map.tags()
    cluster_analysis = maptbx.peak_cluster_analysis(
      peak_list=peak_list,
      special_position_settings=structure,
      general_positions_only=False,
      effective_resolution=d_min,
      min_cross_distance=2,
      max_clusters=n_elements).all()
    check_peaks(
      structure,
      cluster_analysis.sites(),
      cluster_analysis.min_cross_distance() + d_min * grid_resolution_factor)
    structure_from_peaks = xray.structure(structure)
    for site in cluster_analysis.sites():
      structure_from_peaks.add_scatterer(
        xray.scatterer(label="site", scattering_type="", site=site))
    emma_matches = emma.model_matches(
      structure.as_emma_model(),
      structure_from_peaks.as_emma_model(),
      tolerance=d_min*2)
    rms.append(emma_matches.refined_matches[0].rms)
    assert len(emma_matches.refined_matches[0].pairs) == n_elements
  if (0 or verbose):
    print "emma rms grid, interpolated: %.2f %.2f" % tuple(rms)
  assert rms[0] >= rms[1]

def exercise_average_bijvoet_mates(
      space_group_info,
      n_elements=6,
      d_min=5,
      verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=(("O","N","C")*(n_elements//3+1))[:n_elements],
    volume_per_atom=1000,
    min_distance=5,
    general_positions_only=True,
    random_f_double_prime=True)
  if (0 or verbose):
    structure.show_summary().show_scatterers()
  f_calc = structure.structure_factors(
    d_min=d_min, algorithm="direct").f_calc()
  if (0 or verbose):
    f_calc.show_summary()
  assert f_calc.anomalous_flag()
  assert abs(f_calc).anomalous_signal() > 0
  #
  fc_merged = f_calc.average_bijvoet_mates()
  fc_merged_naive = f_calc \
    .as_non_anomalous_array() \
    .merge_equivalents().array().common_set(fc_merged)
  assert fc_merged_naive.indices().all_eq(fc_merged.indices())
  assert flex.max(flex.abs(fc_merged_naive.select_acentric().data()
                          -fc_merged.select_acentric().data())) < 1.e-6
  #
  map_c = f_calc.fft_map().real_map_unpadded()
  map_r = fc_merged.fft_map().real_map_unpadded()
  assert map_c.focus() == map_r.focus()
  lc = flex.linear_correlation(map_c.as_1d(), map_r.as_1d())
  if (0 or verbose):
    print lc.coefficient()
  assert lc.coefficient() > 1-1.e-6
  assert flex.max(flex.abs(map_c.as_1d() - map_r.as_1d())) < 1.e-6
  #
  fc_ma = fc_merged.generate_bijvoet_mates().adopt_set(f_calc)
  fc_mam = fc_ma.average_bijvoet_mates().adopt_set(fc_merged)
  assert flex.max(flex.abs(fc_mam.data()-fc_merged.data())) < 1.e-6

def run_call_back(flags, space_group_info):
  run_test(space_group_info, verbose=flags.Verbose)
  if (not space_group_info.group().is_centric()):
    exercise_average_bijvoet_mates(space_group_info, verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
