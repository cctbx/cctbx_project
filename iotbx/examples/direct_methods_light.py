"""Example of how to use direct methods"""
from __future__ import absolute_import, division, print_function
from six.moves import range
from six.moves import zip
import sys

#
# Required input files:
#   http://journals.iucr.org/c/issues/2001/05/00/vj1132/vj1132Isup2.hkl
#   http://journals.iucr.org/c/issues/2001/05/00/vj1132/vj1132sup1.cif
#
# Command to run this example:
#   iotbx.python direct_methods_light.py vj1132Isup2.hkl vj1132sup1.cif
#

def run(args):
  reflection_file_name = args[0]
  import iotbx.cif
  miller_arrays = iotbx.cif.reader(
    file_path=reflection_file_name).as_miller_arrays()
  for miller_array in miller_arrays:
    s = str(miller_array.info())
    if '_meas' in s:
      if miller_array.is_xray_intensity_array():
        break
      elif miller_array.is_xray_amplitude_array():
        break
  if not ('_meas' in str(miller_array.info())
          and (miller_array.is_xray_amplitude_array()
               or miller_array.is_xray_intensity_array())):
    print("Sorry: CIF does not contain an appropriate miller array")
    return
  miller_array.show_comprehensive_summary()
  print()

  if (miller_array.is_xray_intensity_array()):
    print("Converting intensities to amplitudes.")
    miller_array = miller_array.as_amplitude_array()
    print()

  miller_array.setup_binner(auto_binning=True)
  miller_array.binner().show_summary()
  print()

  all_e_values = miller_array.quasi_normalize_structure_factors().sort(
    by_value="data")
  large_e_values = all_e_values.select(all_e_values.data() > 1.2)
  print("number of large_e_values:", large_e_values.size())
  print()

  from cctbx import dmtbx
  triplets = dmtbx.triplet_generator(large_e_values)
  from cctbx.array_family import flex
  print("triplets per reflection: min,max,mean: %d, %d, %.2f" % (
    flex.min(triplets.n_relations()),
    flex.max(triplets.n_relations()),
    flex.mean(triplets.n_relations().as_double())))
  print("total number of triplets:", flex.sum(triplets.n_relations()))
  print()

  input_phases = large_e_values \
    .random_phases_compatible_with_phase_restrictions()
  tangent_formula_phases = input_phases.data()
  for i in range(10):
    tangent_formula_phases = triplets.apply_tangent_formula(
      amplitudes=large_e_values.data(),
      phases_rad=tangent_formula_phases,
      selection_fixed=None,
      use_fixed_only=False,
      reuse_results=True)

  e_map_coeff = large_e_values.phase_transfer(
    phase_source=tangent_formula_phases)
  from cctbx import maptbx
  e_map = e_map_coeff.fft_map(
    symmetry_flags=maptbx.use_space_group_symmetry)
  e_map.apply_sigma_scaling()
  e_map.statistics().show_summary(prefix="e_map ")
  print()

  peak_search = e_map.peak_search(parameters=maptbx.peak_search_parameters(
    min_distance_sym_equiv=1.2))
  peaks = peak_search.all(max_clusters=10)
  print("e_map peak list")
  print("       fractional coordinates       peak height")
  for site,height in zip(peaks.sites(), peaks.heights()):
    print("  (%9.6f, %9.6f, %9.6f)" % site, "%10.3f" % height)
  print()

  if (len(args) > 1):
    coordinate_file_name = args[1]
    from cctbx import xray
    # xray_structure = xray.structure.from_cif(file_path=coordinate_file_name,
        # data_block_name="I")
    xray_structure = iotbx.cif.reader(
      file_path=coordinate_file_name).build_crystal_structures(
        data_block_name="I")
    xray_structure.show_summary().show_scatterers()
    print()

    f_calc = abs(miller_array.structure_factors_from_scatterers(
      xray_structure=xray_structure,
      algorithm="direct").f_calc())
    correlation = flex.linear_correlation(f_calc.data(), miller_array.data())
    assert correlation.is_well_defined()
    print("correlation of f_obs and f_calc: %.4f" % correlation.coefficient())
    print()

    reference_model = xray_structure.as_emma_model()
    assert reference_model.unit_cell().is_similar_to(e_map.unit_cell())
    assert reference_model.space_group() == e_map.space_group()
    from cctbx import euclidean_model_matching as emma
    peak_model = emma.model(special_position_settings=reference_model)
    for i,site in enumerate(peaks.sites()):
      peak_model.add_position(emma.position(label="peak%02d" % i, site=site))
    matches = emma.model_matches(
      model1=reference_model,
      model2=peak_model,
      tolerance=1.,
      models_are_diffraction_index_equivalent=True)
    for match in matches.refined_matches:
      match.show()

if (__name__ == "__main__"):
  run(sys.argv[1:])
