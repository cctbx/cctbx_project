def run(args):
  import cctbx.omz.dev
  import iotbx.cif
  import cctbx.xray
  #
  master_phil = cctbx.omz.dev.get_master_phil(
    iteration_limit=100,
    additional_phil_string="""\
remove_hydrogen = True
  .type = bool
f_obs_is_f_calc = True
  .type = bool
reset_u_iso = None
  .type = float
""")
  argument_interpreter = master_phil.command_line_argument_interpreter()
  phil_objects = []
  remaining_args = []
  for arg in args:
    if (arg.find("=") >= 0):
      phil_objects.append(argument_interpreter.process(arg=arg))
    else:
      remaining_args.append(arg)
  work_phil = master_phil.fetch(sources=phil_objects)
  work_phil.show()
  params = work_phil.extract()
  print
  #
  assert len(remaining_args) == 2, "refl_cif.hkl model.cif"
  refl_cif = remaining_args[0]
  model_cif = remaining_args[1]
  #
  miller_arrays = iotbx.cif.reader(file_path=refl_cif).as_miller_arrays()
  f_obs = None
  for ma in miller_arrays:
    s = str(ma.info())
    if (s.find("_meas") >= 0):
      if (ma.is_xray_amplitude_array()):
        f_obs = ma
        break
      if (ma.is_xray_intensity_array()):
        f_obs = ma.f_sq_as_f()
        break
  else:
    raise RuntimeError("Missing amplitude array.")
  f_obs.show_comprehensive_summary()
  print
  #
  structure_ideal = cctbx.xray.structure.from_cif(file_path=model_cif).values()[0]
  structure_ideal.show_summary()
  print
  #
  assert structure_ideal.is_similar_symmetry(f_obs)
  #
  structure_ideal.convert_to_isotropic()
  #
  if (params.remove_hydrogen):
    sel = structure_ideal.hd_selection()
    print "Removing hydrogen atoms:", sel.count(True)
    structure_ideal = structure_ideal.select(selection=~sel)
    print
  #
  if (params.f_obs_is_f_calc):
    f_obs = f_obs.structure_factors_from_scatterers(
      xray_structure=structure_ideal,
      algorithm="direct",
      cos_sin_table=False).f_calc().amplitudes()
  #
  structure_shake = structure_ideal.deep_copy_scatterers()
  structure_shake.shake_sites_in_place(rms_difference=params.shake_sites_rmsd)
  if (params.reset_u_iso is None):
    structure_shake.shake_adp(spread=params.shake_adp_spread)
  else:
    structure_shake.set_u_iso(value=params.reset_u_iso)
  #
  cctbx.omz.dev.run_refinement(
    structure_ideal=structure_ideal,
    structure_shake=structure_shake,
    params=params,
    f_obs=f_obs)

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
