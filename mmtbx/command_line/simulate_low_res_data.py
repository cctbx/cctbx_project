
import iotbx.phil
import libtbx.load_env
import libtbx.phil.command_line
from libtbx.str_utils import make_header
from libtbx import easy_pickle
import random
import math
import os
import sys

master_phil = iotbx.phil.parse("""
simulate_data {
  pdb_file = None
    .type = path
  hkl_file = None
    .type = path
  data_label = None
    .type = str
  r_free_label = None
    .type = str
  d_min = 3.0
    .type = float
  r_free_flags_fraction = 0.05
    .type = float
  output_file = low_res.mtz
    .type = path
  crystal_symmetry {
    space_group = None
      .type = space_group
    unit_cell = None
      .type = unit_cell
  }
  modify_pdb {
    remove_waters = True
      .type = bool
    convert_to_isotropic = True
      .type = bool
    set_mean_b_iso = None
      .type = float
    set_wilson_b = False
      .type = bool
  }
  truncate {
    add_b_iso = None
      .type = float
    add_b_aniso = 0 0 0 0 0 0
      .type = floats(size=6)
    add_random_error_percent = None
      .type = float
    remove_cone_around_axis = None
      .type = float
    axis_of_rotation = None
      .type = floats(size=3)
  }
}
fake_data {
  include scope mmtbx.command_line.fmodel.fmodel_from_xray_structure_params_str
}
""", process_includes=True)

def run (args) :
  print "WARNING: this is an experimental program - definitely NOT bug-free."
  print "         Use at your own risk!"
  if (len(args) == 0) :
    print "# full parameters:"
    master_phil.show()
    return
  from iotbx import file_reader
  interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil,
    home_scope="simulate_data")
  pdb_in = None
  pdb_hierarchy = None
  hkl_in = None
  user_phil = []
  for arg in args :
    if os.path.isfile(arg) :
      f = file_reader.any_file(arg)
      if (f.file_type == "pdb") :
        pdb_in = f.file_object
        user_phil.append(interpreter.process(arg="pdb_file=%s" % f.file_name))
      elif (f.file_type == "hkl") :
        hkl_in = f.file_object
        user_phil.append(interpreter.process(arg="hkl_file=%s" % f.file_name))
      elif (f.file_type == "phil") :
        user_phil.append(f.file_object)
    else :
      try :
        arg_phil = interpreter.process(arg=arg)
      except RuntimeError :
        print "ignoring uninterpretable argument '%s'" % arg
      else :
        user_phil.append(arg_phil)
  working_phil = master_phil.fetch(sources=user_phil)
  make_header("Working parameters", out=sys.stdout)
  working_phil.show(prefix="  ")
  params_ = working_phil.extract()
  params = params_.simulate_data
  if (params.pdb_file is None) and (params.hkl_file is None) :
    raise Sorry("No PDB file specified.")
  if (pdb_in is None) and (params.pdb_file is not None) :
    f = file_reader.any_file(params.pdb_file, force_type="pdb")
    f.assert_file_type("pdb")
    pdb_in = f.file_object
  if (hkl_in is None) and (params.hkl_file is not None) :
    f = file_reader.any_file(params.hkl_File, force_type="hkl")
    f.assert_file_type("hkl")
    hkl_in = f.file_object
  if (pdb_in is not None) :
    pdb_hierarchy = pdb_in.construct_hierarchy()
  if (hkl_in is not None) :
    make_header("Extracting experimental data", out=sys.stdout)
    F, r_free = from_hkl(hkl_in, params)
  elif (pdb_in is not None) :
    make_header("Generating fake data with phenix.fmodel", out=sys.stdout)
    F, r_free = from_pdb(pdb_in, pdb_hierarchy, params_)
  make_header("Applying low-resolution filtering", out=sys.stdout)
  print "  Final resolution: %.2f A" % params.d_min
  n_residues, n_bases = None, None
  if (pdb_in is not None) :
    n_residues, n_bases = get_counts(pdb_hierarchy)
  #if (params.auto_adjust) :
  #  if (pdb_in is None) :
  #    raise Sorry("You must supply a PDB file when auto_adjust=True.")
  F_out = truncate_data(F, params)
  print "  Completeness after processing: %.2f%%" % (F.completeness() * 100.)
  if (pdb_in is not None) :
    iso_scale, aniso_scale = wilson_scaling(F_out, n_residues, n_bases)
    print ""
    print "  Scaling statistics for output data:"
    print "    overall isotropic B-factor:   %6.2f" % iso_scale.b_wilson
    print "    overall anisotropic B-factor: %6.2f, %6.2f, %6.2f" % (
      aniso_scale.b_cart[0], aniso_scale.b_cart[3], aniso_scale.b_cart[4])
    print "                                  %14.2f, %6.2f" % (
      aniso_scale.b_cart[1], aniso_scale.b_cart[5])
    print "                                  %22.2f" % aniso_scale.b_cart[2]
    print ""
  if (F_out.sigmas() is not None) :
    mtz_dataset = F_out.as_mtz_dataset(
      column_root_label="F",
      column_types="FQ")
  else :
    mtz_dataset = F_out.as_mtz_dataset(
      column_root_label="F",
      column_types="F")
  if (r_free is not None) :
    r_free = r_free.common_set(F_out)
    mtz_dataset.add_miller_array(
      miller_array=r_free,
      column_root_label="FreeR_flag",
      column_types="I")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name=params.output_file)
  make_header("Writing output file", out=sys.stdout)
  print "  Wrote %s" % params.output_file

def from_pdb (pdb_in, pdb_hierarchy, params_) :
  params = params_.simulate_data
  pdb_sg, pdb_uc = None, None
  pdb_symm = pdb_in.crystal_symmetry()
  if (pdb_symm is not None) :
    pdb_sg = pdb_symm.space_group_info()
    pdb_uc = pdb_symm.unit_cell()
  apply_sg = pdb_sg
  apply_uc = pdb_uc
  if (params.crystal_symmetry.space_group is not None) :
    apply_sg = params.crystal_symmetry.space_group
  else :
    params.crystal_symmetry.space_group = pdb_sg
  if (params.crystal_symmetry.unit_cell is not None) :
    apply_uc = params.crystal_symmetry.unit_cell
  else :
    params.crystal_symmetry.unit_cell = pdb_uc
  if (apply_sg is None) or (apply_uc is None) :
    raise Sorry("Incomplete symmetry information - please specify a space "+
      "group and unit cell for this structure.")
  from cctbx import crystal, adptbx
  from scitbx.array_family import flex
  apply_symm = crystal.symmetry(
    unit_cell=apply_uc,
    space_group_info=apply_sg)
  if (params.modify_pdb.remove_waters) :
    print "  Removing solvent atoms..."
    for model in pdb_hierarchy.models() :
      for chain in model.chains() :
        for residue_group in chain.residue_groups() :
          for atom_group in residue_group.atom_groups() :
            if (atom_group.resname in ["HOH", "WAT"]) :
              residue_group.remove_atom_group(atom_group=atom_group)
        if (len(chain.atoms()) == 0) :
          model.remove_chain(chain=chain)
  xray_structure = pdb_in.xray_structure_simple(
    crystal_symmetry=apply_symm)
  sctr_keys = xray_structure.scattering_type_registry().type_count_dict().keys()
  if (not (("H" in sctr_keys) or ("D" in sctr_keys))) :
    print "  WARNING: this model does not contain hydrogen atoms!"
    print "           strongly recommend running phenix.ready_set or "
    print "           equivalent to ensure realistic simulated data."
    print ""
  if (params.modify_pdb.convert_to_isotropic) :
    xray_structure.convert_to_isotropic()
  set_b = None
  if (params.modify_pdb.set_mean_b_iso is not None) :
    assert (not params.modify_pdb.set_wilson_b)
    print "  Scaling B-factors to have mean of %.2f" % \
      params.modify_pdb.set_mean_b_iso
    assert (params.modify_pdb.set_mean_b_iso > 0)
    set_b = params.modify_pdb.set_mean_b_iso
  elif (params.modify_pdb.set_wilson_b) :
    print "  Scaling B-factors to match mean Wilson B for this resolution"
    set_b = get_mean_statistic_for_resolution(
      d_min=params.d_min,
      stat_type="wilson_b")
    print ""
  if (set_b is not None) :
    u_iso = xray_structure.extract_u_iso_or_u_equiv()
    u_mean = flex.mean(u_iso)
    b_mean = adptbx.u_as_b(u_mean)
    scale = set_b / b_mean
    xray_structure.set_u_iso(values=u_iso * scale)
  import mmtbx.command_line.fmodel
  from mmtbx import utils
  fmodel_params = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params.extract()
  fmodel_params.high_resolution = params.d_min
  fmodel_params.fmodel = params_.fake_data.fmodel
  if (fmodel_params.fmodel.b_sol == 0) :
    print "  b_sol is zero - will use mean value for d_min +/- 0.2A"
    print "   (this is not strongly correlated with resolution, but it's good "
    print "    to use a real value instead of leaving it set to 0)"
    fmodel_params.fmodel.b_sol = get_mean_statistic_for_resolution(
      d_min=params.d_min,
      stat_type="b_sol")
    print ""
  if (fmodel_params.fmodel.k_sol == 0) :
    print "  k_sol is zero - will use mean value for d_min +/- 0.2A"
    print "   (this is not strongly correlated with resolution, but it's good "
    print "    to use a real value instead of leaving it set to 0)"
    fmodel_params.fmodel.k_sol = get_mean_statistic_for_resolution(
      d_min=params.d_min,
      stat_type="k_sol")
    print ""
  fmodel_params.structure_factors_accuracy = params_.fake_data.structure_factors_accuracy
  fmodel_params.mask = params_.fake_data.mask
  fmodel_params.r_free_flags_fraction = params.r_free_flags_fraction
  fmodel_params.add_sigmas = False
  fmodel_params.output.type = "real"
  fmodel_ = utils.fmodel_from_xray_structure(
    xray_structure=xray_structure,
    params=fmodel_params)
  f_model = fmodel_.f_model
  r_free_flags = fmodel_.r_free_flags
  return (f_model, r_free_flags)

def get_counts (hierarchy) :
  n_residues = 0
  n_bases = 0
  for chain in hierarchy.models()[0].chains() :
    main_conf = chain.conformers()[0]
    if (main_conf.is_protein()) :
      n_residues += len(chain.residue_groups())
    elif (main_conf.is_na()) :
      n_bases += len(chain.residue_groups())
  return n_residues, n_bases

def wilson_scaling (F, n_residues, n_bases) :
  from mmtbx.scaling import absolute_scaling
  iso_scale = absolute_scaling.ml_iso_absolute_scaling(
    miller_array=F,
    n_residues=n_residues,
    n_bases=n_bases)
  aniso_scale = absolute_scaling.ml_aniso_absolute_scaling(
    miller_array=F,
    n_residues=n_residues,
    n_bases=n_bases)
  return iso_scale, aniso_scale

def from_hkl (hkl_in, params) :
  miller_arrays = hkl_in.as_miller_arrays()
  f_obs = None
  r_free = None
  for array in miller_arrays :
    if (params.data_label is not None) :
      if (array.info().label_string() == params.data_label) :
        f_obs = array.map_to_asu()
    elif (array.is_xray_amplitude_array()) :
      f_obs = array.map_to_asu()
      params.data_label = array.info().label_string()
      print "  F-obs: %s" % params.data_label
    elif (array.is_xray_intensity_array()) :
      f_obs = array.f_sq_as_f().map_to_asu()
      params.data_label = array.info().label_string()
      print "  I-obs: %s" % params.data_label
    if (params.r_free_label is not None) :
      if (array.info().label_string() == params.r_free_label) :
        r_free = array.map_to_asu()
        print "  R-free: %s" % array.info().label_string()
    elif (array.info().label_string() in ["FreeR_flag", "FREE"]) :
      r_free = array.map_to_asu()
      print "  R-free: %s" % array.info().label_string()
  f_obs = f_obs.resolution_filter(d_min=params.d_min)
  if (r_free is not None) :
    r_free = r_free.common_set(f_obs)
  return f_obs, r_free

def truncate_data (F, params) :
  from scitbx.array_family import flex
  if (params.truncate.add_b_iso is not None) :
    print "  Applying isotropic B-factor of %.2f A^2" % (
      params.truncate.add_b_iso)
    F = add_b_iso(f_obs=F, b_iso=params.truncate.add_b_iso)
  if (params.truncate.add_b_aniso != [0,0,0,0,0,0]) :
    print "  Adding anisotropy..."
    F = add_b_aniso(f_obs=F, b_cart=params.truncate.add_b_aniso)
  if (params.truncate.add_random_error_percent is not None) :
    print "  Adding random error as percent of amplitude..."
    F = add_random_error(f_obs=F,
      error_percent=params.truncate.add_random_error_percent)
  if (params.truncate.remove_cone_around_axis is not None) :
    print "  Removing cone of data around axis of rotation..."
    print "    radius = %.1f degrees" % params.truncate.remove_cone_around_axis
    assert (params.truncate.axis_of_rotation is not None)
    remove_cone_around_axis(
      f_obs=F,
      axis=params.truncate.axis_of_rotation,
      cone_radius=params.truncate.remove_cone_around_axis)
  return F

def add_b_iso (f_obs, b_iso) :
  from scitbx.array_family import flex
  assert (b_iso > 0)
  d_min_sq = flex.pow2(f_obs.d_spacings().data())
  scale = flex.exp(- b_iso / d_min_sq / 4.0)
  return f_obs.customized_copy(data=f_obs.data() * scale)

def add_b_aniso (f_obs, b_cart) :
  from cctbx import adptbx
  from scitbx.array_family import flex
  u_star = adptbx.u_cart_as_u_star(
    f_obs.unit_cell(), adptbx.b_as_u(b_cart))
  scale = flex.double()
  for i, hkl in enumerate(f_obs.indices()) :
    scale.append(f_aniso_one_h(hkl, u_star))
  return f_obs.customized_copy(data=f_obs.data() * scale)

pi_sq = math.pi**2
# see mmtbx/f_model/f_model.h
def f_aniso_one_h (h, u_star) :
  arg = -2.0 * pi_sq * (
    u_star[0] * h[0] * h[0] +
    u_star[1] * h[1] * h[1] +
    u_star[2] * h[2] * h[2] +
    u_star[3] * h[0] * h[1] * 2.0 +
    u_star[4] * h[0] * h[2] * 2.0 +
    u_star[5] * h[1] * h[2] * 2.0)
  if (arg > 40.0) : arg = 40.0 # ???
  return math.exp(arg)

def add_random_error (f_obs, error_percent) :
  from scitbx.array_family import flex
  assert (100 > error_percent > 0)
  data = f_obs.data()
  fr = F.data() * error_percent / 100.
  ri = flex.double()
  for trial in xrange(data.size()):
    r = random.randint(0,1)
    if(r == 0): r = -1
    ri.append(r)
  data = data + ri*fr
  return f_obs.customized_copy(data=data)

def remove_cone_around_axis (f_obs, axis, cone_radius) :
  assert (cone_radius < 90) and (cone_radius >= 0)
  from scitbx.matrix import rec
  indices = f_obs.indices()
  n_hkl = indices.size()
  data = f_obs.data()
  sigmas = f_obs.sigmas()
  v1 = rec(axis, (3,1))
  angle_high = 180 - cone_radius
  i = 0
  while (i < len(indices)) :
    hkl = indices[i]
    v2 = rec(hkl, (3,1))
    angle = math.degrees(v1.angle(v2))
    if (angle <= cone_radius) or (angle >= angle_high) :
      del indices[i]
      del data[i]
      if (sigmas is not None) :
        del sigmas[i]
    else :
      i += 1
  delta_n_hkl = n_hkl - indices.size()
  print "    removed %d reflections (out of %d)" % (delta_n_hkl, n_hkl)

stat_names = {
  'wilson_b' : 'Wilson B-factor',
  'k_sol' : 'Bulk solvent scale factor (k_sol)',
  'b_sol' : 'Bulk solvent B-factor (b_sol)',
}

def get_mean_statistic_for_resolution (d_min, stat_type, range=0.2) :
  from scitbx.array_family import flex
  pkl_file = libtbx.env.find_in_repositories(
    relative_path = "chem_data/polygon_data/all_mvd.pickle",
    test = os.path.isfile)
  db = easy_pickle.load(pkl_file)
  all_d_min = db['high_resolution']
  stat_values = db[stat_type]
  values_for_range = flex.double()
  for (d_, v_) in zip(all_d_min, stat_values) :
    try :
      d = float(d_)
      v = float(v_)
    except ValueError : continue
    else :
      if (d > (d_min - range)) and (d < (d_min + range)) :
        values_for_range.append(v)
  h = flex.histogram(values_for_range, n_slots=10)
  print "  %s for d_min = %.3f - %.3f A" % (stat_names[stat_type], d_min-range,
    d_min+range)
  min = flex.min(values_for_range)
  max = flex.max(values_for_range)
  mean = flex.mean(values_for_range)
  print "    count: %d" % values_for_range.size()
  print "    min: %.2f" % min
  print "    max: %.2f" % max
  print "    mean: %.2f" % mean
  print "    histogram of values:"
  h.show(prefix="      ")
  return mean

if (__name__ == "__main__") :
  run(sys.argv[1:])
