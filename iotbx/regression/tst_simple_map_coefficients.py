
from __future__ import absolute_import, division, print_function

def exercise():
  from iotbx import file_reader
  from cctbx.array_family import flex
  from cctbx import crystal
  from cctbx import xray
  from libtbx import easy_run
  file_base = "tmp_iotbx_simple_map_coefficients"
  quartz_structure = xray.structure(
    special_position_settings=crystal.special_position_settings(
      crystal_symmetry=crystal.symmetry(
        unit_cell=(5.01,5.01,5.47,90,90,120),
        space_group_symbol="P6222")),
    scatterers=flex.xray_scatterer([
      xray.scatterer(
        label="Si",
        site=(1/2.,1/2.,1/3.),
        u=0.2,
        fp=0.255, # from Sasaki table @ 1.54 Angstrom
        fdp=0.325),
      xray.scatterer(
        label="O",
        site=(0.197,-0.197,0.83333),
        u=0)]))
  fc = quartz_structure.structure_factors(d_min=1.0).f_calc()
  I = fc.amplitudes().f_as_f_sq()
  PHI = fc.average_bijvoet_mates().phases(deg=True).set_observation_type(None)
  FOM = PHI.customized_copy(
    data=flex.double(PHI.data().size(), 1.0)).set_observation_type(None)
  mtz = I.as_mtz_dataset(column_root_label="I")
  mtz.add_miller_array(PHI,
    column_root_label="PHI",
    column_types="P")
  mtz.add_miller_array(FOM,
    column_root_label="FOM",
    column_types="W")
  mtz_file = file_base + ".mtz"
  mtz.mtz_object().write(mtz_file)
  map_file = "%s_map.mtz" % file_base
  # Fo map
  args = [
    mtz_file,
    "output_file=%s" % map_file,
  ]
  cmd = "iotbx.simple_map_coefficients"
  result = easy_run.fully_buffered("%s %s" % (cmd, " ".join(args))
    ).raise_if_errors()
  assert ("Applying weights in FOM" in result.stdout_lines)
  map_in = file_reader.any_file(map_file).assert_file_type("hkl")
  map_coeffs = map_in.file_server.miller_arrays[0]
  real_map = map_coeffs.fft_map().apply_sigma_scaling().real_map_unpadded()
  site = quartz_structure.sites_frac()[0]
  assert (real_map.eight_point_interpolation(site) > 3)
  # anomalous map
  args.append("map_type=anom")
  result = easy_run.fully_buffered("%s %s" % (cmd, " ".join(args))
    ).raise_if_errors()
  assert (not "Applying weights in FOM" in result.stdout_lines)
  map_in = file_reader.any_file(map_file).assert_file_type("hkl")
  map_coeffs = map_in.file_server.miller_arrays[0]
  real_map = map_coeffs.fft_map().apply_sigma_scaling().real_map_unpadded()
  assert (real_map.eight_point_interpolation(site) > 3)
  # no FOM
  mtz = I.as_mtz_dataset(column_root_label="I")
  mtz.add_miller_array(PHI,
    column_root_label="PHI",
    column_types="P")
  mtz.mtz_object().write(mtz_file)
  args = [
    mtz_file,
    "output_file=%s" % map_file,
  ]
  result = easy_run.fully_buffered("%s %s" % (cmd, " ".join(args))
    ).raise_if_errors()
  assert (not "Applying weights in FOM" in result.stdout_lines)
  map_in = file_reader.any_file(map_file).assert_file_type("hkl")
  map_coeffs = map_in.file_server.miller_arrays[0]
  real_map = map_coeffs.fft_map().apply_sigma_scaling().real_map_unpadded()
  site = quartz_structure.sites_frac()[0]
  assert (real_map.eight_point_interpolation(site) > 3)
  # now check for error when use_weights=True
  args.append("use_weights=True")
  result = easy_run.fully_buffered("%s %s" % (cmd, " ".join(args)))
  assert (result.return_code == 1)
  print("OK")

if (__name__ == "__main__"):
  exercise()
