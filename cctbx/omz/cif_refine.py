def report_fraction_of_negative_observations_if_any(id_code, obs):
  d = obs.data()
  n_neg = (d < 0).count(True)
  n_all = d.size()
  if (n_neg != 0 and n_all != 0):
    print obs.info()
    print "fraction_of_negative_obs: %.6g (%d of %d)" % (
      n_neg / n_all, n_neg, n_all)
    neg = d.select(d < 0)
    pos = d.select(d >= 0)
    from cctbx.array_family import flex
    def hist(data):
      from cStringIO import StringIO
      sio = StringIO()
      flex.histogram(data=data, n_slots=10) \
        .show(f=sio, prefix="  ", format_cutoffs="%8.2f")
      return sio.getvalue().splitlines()
    lines_neg = hist(-neg)
    lines_pos = hist(pos)
    pair_fmt = "%-35s | %s"
    print pair_fmt % (
      "Histogram of negative observations:",
        "positive observations:")
    for pair in zip(lines_neg, lines_pos):
      print pair_fmt % pair

class extract_from_cif_files(object):

  __slots__ = [
    "c_obs",
    "xray_structure",
    "non_hydrogen_selection",
    "non_hydrogen_iselection",
    "edge_list"]

  def init_slots_with_none(O):
    for slot in O.__slots__:
      setattr(O, slot, None)

  def __init__(O, report_id, refl_file, refl_cif, model_file, model_cif):
    O.init_slots_with_none()
    O.process(report_id, refl_file, refl_cif, model_file, model_cif)

  def process(O, report_id, refl_file, refl_cif, model_file, model_cif):
    from_coordinate_files = []
    from_reflection_files = []
    import iotbx.cif.builders
    def get_cs(cif, buffer):
      for cif_block in cif.model().values():
        cs = iotbx.cif.builders.crystal_symmetry_builder(
          cif_block=cif_block).crystal_symmetry
        buffer.append(cs)
    get_cs(refl_cif, from_reflection_files)
    get_cs(model_cif, from_coordinate_files)
    import cctbx.crystal
    combined_cs = cctbx.crystal.select_crystal_symmetry(
      from_coordinate_files=from_coordinate_files,
      from_reflection_files=from_reflection_files)
    if (combined_cs.unit_cell() is None):
      raise RuntimeError("Unit cell not found in both cif and hkl files.")
    if (combined_cs.space_group_info() is None):
      raise RuntimeError("Space group not found in both cif and hkl file.")
    #
    miller_arrays = refl_cif.as_miller_arrays()
    meas_a = []
    meas_i = []
    for ma in miller_arrays:
      s = str(ma.info())
      if (s.find("_meas") >= 0):
        if (ma.is_xray_intensity_array()):
          meas_i.append(ma)
        elif (ma.is_xray_amplitude_array()):
          meas_a.append(ma)
        else:
          continue
        report_fraction_of_negative_observations_if_any(report_id, ma)
    if (len(meas_i) != 0):
      O.c_obs = meas_i[0]
    elif (len(meas_a) != 0):
      O.c_obs = meas_a[0]
    else:
      raise RuntimeError("Missing diffraction data.")
    O.c_obs = O.c_obs.customized_copy(
      crystal_symmetry=combined_cs)
    print "."*79
    O.c_obs.show_comprehensive_summary()
    print "."*79
    #
    O.xray_structure = cctbx.xray.structure.from_cif(
      file_path=model_file).values()[0]
    O.xray_structure = O.xray_structure.customized_copy(
      crystal_symmetry=combined_cs)
    O.xray_structure.show_summary().show_scatterers()
    print "."*79
    #
    O.non_hydrogen_selection = ~O.xray_structure.hd_selection()
    O.non_hydrogen_iselection = O.non_hydrogen_selection.iselection()
    #
    O.edge_list = O.process_geom_bond(model_cif)
    if (O.edge_list is not None):
      print "len(edge_list):", len(O.edge_list), report_id

  def process_geom_bond(O, model_cif):
    xs = O.xray_structure
    scs = xs.scatterers()
    i_seq_by_lbl = dict(zip(scs.extract_labels(), range(scs.size())))
    if (len(i_seq_by_lbl) != scs.size()):
      return None
    edge_set = set()
    for cif_block in model_cif.model().values():
      lbl_lists = [cif_block.get("_geom_bond_atom_site_label_"+s)
        for s in "12"]
      if (lbl_lists.count(None) != 0):
        return None
      if (len(lbl_lists[0]) != len(lbl_lists[1])):
        return None
      for lbl_pair in zip(*lbl_lists):
        i_seqs = tuple(sorted([i_seq_by_lbl.get(lbl) for lbl in lbl_pair]))
        if (i_seqs.count(None) != 0):
          return None
        if (i_seqs in edge_set):
          return None
        edge_set.add(i_seqs)
    return sorted(edge_set)

def run(args):
  import cctbx.omz.dev
  import iotbx.cif
  import cctbx.xray
  #
  import random
  random.seed(0)
  from cctbx.array_family import flex
  flex.set_random_seed(0)
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
  refl_file = remaining_args[0]
  model_file = remaining_args[1]
  refl_cif = iotbx.cif.reader(file_path=refl_file)
  model_cif = iotbx.cif.reader(file_path=model_file)
  #
  mgr = extract_from_cif_files(
    report_id="cif_refine",
    refl_file=refl_file,
    refl_cif=refl_cif,
    model_file=model_file,
    model_cif=model_cif)
  #
  structure_ideal = mgr.xray_structure.deep_copy_scatterers()
  structure_ideal.convert_to_isotropic()
  if (params.remove_hydrogen):
    sel = mgr.non_hydrogen_selection
    print "Removing hydrogen atoms:", sel.count(False)
    structure_ideal = structure_ideal.select(selection=sel)
    print
  #
  c_obs = mgr.c_obs
  if (params.f_obs_is_f_calc):
    c_obs = c_obs.structure_factors_from_scatterers(
      xray_structure=structure_ideal,
      algorithm="direct",
      cos_sin_table=False).f_calc().amplitudes() \
        .set_observation_type_xray_amplitude()
  #
  assert c_obs.is_xray_intensity_array() or c_obs.is_xray_amplitude_array()
  if (c_obs.is_xray_intensity_array()):
    i_obs = c_obs
    f_obs = c_obs.f_sq_as_f(algorithm="xtal_3_7")
  else:
    f_obs = c_obs
    i_obs = c_obs.f_as_f_sq(algorithm="shelxl")
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
    f_obs=f_obs,
    i_obs=i_obs)

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
