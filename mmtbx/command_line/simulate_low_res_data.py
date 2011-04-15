
import iotbx.phil
import libtbx.load_env
import libtbx.phil.command_line
from libtbx.str_utils import make_header
from libtbx.math_utils import ifloor
from libtbx import easy_pickle
import random
import math
import os
import sys

master_phil = iotbx.phil.parse("""
simulate_data {
  pdb_file = None
    .type = path
    .help = Model file.  If no reflections file is given, data will be \
      generated starting from F(model).  Otherwise it will only be used \
      to report scaling statistics.
  hkl_file = None
    .type = path
    .help = Data file.  If defined, the extracted amplitudes will be \
      truncated.  Any R-free flags will be propagated to the output file.
  data_label = None
    .type = str
    .help = Label for amplitudes or intensities in hkl_file (optional).
  d_min = 3.0
    .type = float
    .help = Resolution cutoff of output data.
  random_seed = 90125714
    .type = int
  output_file = None
    .type = path
  r_free_flags {
    file_name = None
      .type = path
      .help = File containing R-free flags.  Can be left blank if hkl_file is \
        defined and contains R-free flags.
    label = None
      .type = str
      .help = Label for R-free flags in r_free_file (or hkl_file).
    fraction = 0.05
      .type = float
      .help = Percent of reflections to flag for R-free (ignored if already \
        available).
  }
  crystal_symmetry {
    space_group = None
      .type = space_group
      .help = Space group of output data.  Ignored if reflections are used \
        as input.
    unit_cell = None
      .type = unit_cell
      .help = Unit cell of output data.  Ignored if reflections are used as \
        input.
  }
  modify_pdb {
    remove_waters = True
      .type = bool
      .help = Strip waters from input model.
    remove_alt_confs = True
      .type = bool
      .help = Strip alternate sidechains (and reset occupancy to 1).
    convert_to_isotropic = False
      .type = bool
      .help = Convert all atoms to isotropic before calculating F(model).
    set_mean_b_iso = None
      .type = float
      .help = Scale atomic B-factors to have this mean value.
    set_wilson_b = False
      .type = bool
      .help = Scale atomic B-factors to have a mean equal to the mean \
        Wilson B-factor for this resolution (+/- 0.2A)
  }
  truncate
    .help = Options for data truncation at lower resolution.
  {
    add_b_iso = None
      .type = float
      .help = Isotropic B-factor to be added to data.
    add_b_aniso = 0 0 0 0 0 0
      .type = floats(size=6)
      .help = Anisotropic B-factor to be added to data.  Severely anisotropic \
        data might have an anisotropic B of 80,80,200,0,0,0.
    add_random_error_percent = None
      .type = float
      .help = Adds random noise as a percentage of amplitude, evenly across \
        all resolutions.  This is probably inferior to the sigma-based \
        noise generation.
    remove_cone_around_axis = None
      .type = float
      .help = Radius in degrees of cone of missing data around the axis of \
        rotation (data collection pathology).  If specified, axis_of_rotation \
        must be defined.
    axis_of_rotation = None
      .type = floats(size=3)
      .help = Axis of rotation of crystal during data collection.  Only used \
        if remove_cone_around_axis is defined.  (Example value: 0,1,0)
  }
  generate_noise {
    add_noise = False
      .type = bool
    noise_profile_file = None
      .type = path
    profile_model_file = None
      .type = path
    profile_data_label = None
      .type = str
    scale_noise = 1.0
      .type = float
    n_resolution_bins = 20
      .type = int
    n_intensity_bins = 20
      .type = int
  }
  fake_data_from_fmodel
    .help = Options for generating model-based reflections using phenix.fmodel.
  {
    include scope mmtbx.command_line.fmodel.fmodel_from_xray_structure_params
  }
}
""", process_includes=True)

def run (args) :
  print ""
  print "mmtbx.simulate_low_res_data"
  print "  For generation of realistic data (model-based, or using real "
  print "  high-resolution data) for methods development."
  print ""
  print "************************* WARNING: *************************"
  print " this is an experimental program - definitely NOT bug-free."
  print "                  Use at your own risk!"
  print ""
  print "  Usage:"
  print "   mmtbx.simulate_low_res_data model.pdb [options...]"
  print "     (generate data from a PDB file)"
  print ""
  print "   mmtbx.simulate_low_res_data highres.mtz [model.pdb] [options...]"
  print "     (truncate high-resolution data)"
  print ""
  print "   mmtbx.simulate_low_res_data --help"
  print "     (print full parameters with additional info)"
  print ""
  if (len(args) == 0) or ("--help" in args) :
    print "# full parameters:"
    if ("--help" in args) :
      master_phil.show(attributes_level=1)
    else :
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
  if (params.generate_noise.add_noise) and (params.hkl_file is None) :
    if (params.generate_noise.noise_profile_file is None) :
      raise Sorry("noise_profile_file required when add_noise=True and "
        "hkl_file is undefined.")
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
    F, r_free = from_pdb(pdb_in, pdb_hierarchy, params)
  if (params.r_free_flags.file_name is not None) :
    rfree_in = file_reader.any_file(params.r_free_flags.file_name)
    rfree_in.assert_file_type("hkl")
    hkl_server = rfree_in.file_server
    r_free_raw, flag_value = hkl_server.get_r_free_flags(
      file_name=None,
      label=params.r_free_flags.label,
      test_flag_value=None,
      parameter_scope="simulate_data.r_free_flags",
      disable_suitability_test=False)
    r_free = r_free_raw.customized_copy(data=r_free_raw.data() == flag_value)
    r_free = r_free.common_set(F)
    if (F.data().size() != r_free.data().size()) :
      raise Sorry(("The specified R-free flags in %s are incomplete.  Please "+
        "generate a complete set to the desired resolution limit if you "+
        "want to use these flags with synthetic data.") %
          params.r_free_flags.file_name)
  make_header("Applying low-resolution filtering", out=sys.stdout)
  print "  Final resolution: %.2f A" % params.d_min
  n_residues, n_bases = None, None
  if (pdb_in is not None) :
    n_residues, n_bases = get_counts(pdb_hierarchy)
  #if (params.auto_adjust) :
  #  if (pdb_in is None) :
  #    raise Sorry("You must supply a PDB file when auto_adjust=True.")
  F_out = truncate_data(F, params)
  if (params.generate_noise.add_noise) :
    make_header("Adding noise using sigma profile", out=sys.stdout)
    if (F_out.sigmas() is None) :
      if (pdb_in is not None) :
        iso_scale, aniso_scale = wilson_scaling(F_out, n_residues, n_bases)
      i_obs = create_sigmas(
        f_obs=F_out,
        params=params.generate_noise,
        wilson_b=iso_scale.b_wilson,
        return_as_amplitudes=False)
    apply_sigma_noise(i_obs)
    F_out = i_obs.f_sq_as_f()
  make_header("Done processing", out=sys.stdout)
  print "  Completeness after processing: %.2f%%" % (F.completeness() * 100.)
  if (pdb_in is not None) :
    iso_scale, aniso_scale = wilson_scaling(F_out, n_residues, n_bases)
    print ""
    print "  Scaling statistics for output data:"
    show_b_factor_info(iso_scale, aniso_scale)
    print ""
    wilson_b = iso_scale.b_wilson
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
  if (params.output_file is None) :
    if (params.hkl_file is not None) :
      base_name = os.path.splitext(os.path.basename(params.hkl_file))[0]
    else :
      base_name = os.path.splitext(os.path.basename(params.pdb_file))[0]
    params.output_file = base_name + "_low_res.mtz"
  mtz_object.write(file_name=params.output_file)
  print "  Wrote %s" % params.output_file

def from_pdb (pdb_in, pdb_hierarchy, params) :
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
          if (len(residue_group.atom_groups()) == 0) :
            chain.remove_residue_group(residue_group=residue_group)
        if (len(chain.atoms()) == 0) :
          model.remove_chain(chain=chain)
  if (params.modify_pdb.remove_alt_confs) :
    print "  Removing all alternate conformations and resetting occupancies..."
    for model in pdb_hierarchy.models() :
      for chain in model.chains() :
        for residue_group in chain.residue_groups() :
          atom_groups = residue_group.atom_groups()
          assert (len(atom_groups) > 0)
          if (len(atom_groups) == 1) : continue
          for atom_group in atom_groups[1:] :
            residue_group.remove_atom_group(atom_group=atom_group)
          atom_groups[0].altloc = ""
    atoms = pdb_hierarchy.atoms()
    new_occ = flex.double(atoms.size(), 0.0)
    atoms.set_occ(new_occ)
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
  fake_data = params.fake_data_from_fmodel
  fmodel_params.fmodel = fake_data.fmodel
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
  fmodel_params.structure_factors_accuracy = fake_data.structure_factors_accuracy
  fmodel_params.mask = fake_data.mask
  fmodel_params.r_free_flags_fraction = params.r_free_flags.fraction
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
    if (params.r_free_flags.label is not None) :
      if (array.info().label_string() == params.r_free_flags.label) :
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
  if (f_obs.sigmas() is not None) :
    return f_obs.customized_copy(
      data=f_obs.data() * scale,
      sigmas=f_obs.sigmas() * scale)
  else :
    return f_obs.customized_copy(data=f_obs.data() * scale)

def add_b_aniso (f_obs, b_cart) :
  from cctbx import adptbx
  from scitbx.array_family import flex
  u_star = adptbx.u_cart_as_u_star(
    f_obs.unit_cell(), adptbx.b_as_u(b_cart))
  scale = flex.double()
  for i, hkl in enumerate(f_obs.indices()) :
    scale.append(f_aniso_one_h(hkl, u_star))
  if (f_obs.sigmas() is not None) :
    return f_obs.customized_copy(
      data=f_obs.data() * scale,
      sigmas=f_obs.sigmas() * scale)
  else :
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

def show_b_factor_info (iso_scale, aniso_scale) :
  print "    overall isotropic B-factor:   %6.2f" % iso_scale.b_wilson
  print "    overall anisotropic B-factor: %6.2f, %6.2f, %6.2f" % (
    aniso_scale.b_cart[0], aniso_scale.b_cart[3], aniso_scale.b_cart[4])
  print "                                  %14.2f, %6.2f" % (
    aniso_scale.b_cart[1], aniso_scale.b_cart[5])
  print "                                  %22.2f" % aniso_scale.b_cart[2]

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

def create_sigmas (f_obs, params, wilson_b=None, return_as_amplitudes=False) :
  assert (f_obs.sigmas() is None)
  from scitbx.array_family import flex
  i_obs = f_obs.f_as_f_sq()
  i_norm = i_obs.data() / flex.max(i_obs.data())
  profiler = profile_sigma_generator(
    mtz_file=params.noise_profile_file,
    pdb_file=params.profile_model_file,
    wilson_b=wilson_b,
    data_label=params.profile_data_label,
    n_resolution_bins=params.n_resolution_bins,
    n_intensity_bins=params.n_intensity_bins)
  sigmas = flex.double(i_norm.size(), 0.0)
  i_obs.setup_binner(n_bins=params.n_resolution_bins)
  for j_bin in i_obs.binner().range_used() :
    bin_sel = i_obs.binner().selection(j_bin)
    shell_profile = profiler.get_noise_profile_for_shell(j_bin)
    for k in bin_sel.iselection() :
      i_over_sigma = shell_profile.get_i_over_sigma(i_norm[k])
      sigmas[k] = i_obs.data()[k] / i_over_sigma
  i_new = i_obs.customized_copy(sigmas=sigmas)
  if (return_as_amplitudes) :
    return i_new.f_as_f_sq()
  else :
    return i_new

def apply_sigma_noise (i_obs) :
  assert (i_obs.is_xray_intensity_array())
  assert (i_obs.sigmas() is not None)
  data = i_obs.data()
  for i, sigma in enumerate(i_obs.sigmas()) :
    data[i] = random.gauss(data[i], sigma)

class profile_sigma_generator (object) :
  def __init__ (self,
                mtz_file,
                pdb_file,
                wilson_b=None,
                data_label=None,
                n_resolution_bins=20,
                n_intensity_bins=20) :
    if (wilson_b is None) or (pdb_file is None) :
      print "  WARNING: missing desired Wilson B-factor and/or PDB file"
      print "           for noise profile data.  Without this information"
      print "           the intensity falloff with resolution will probably"
      print "           not be the same for your synthetic data and the"
      print "           data used to generate sigmas."
    self._resolution_bins = []
    from iotbx.file_reader import any_file
    from scitbx.array_family import flex
    f = any_file(mtz_file, force_type="hkl")
    f.assert_file_type("hkl")
    miller_arrays = f.file_server.miller_arrays
    f_obs = None
    i_obs = None
    for array in miller_arrays :
      if (array.info().label_string() == data_label) or (data_label is None) :
        if (array.is_xray_amplitude_array()) and (f_obs is None) :
          f_obs = array
        elif (array.is_xray_intensity_array()) and (i_obs is None) :
          i_obs = array
    if (i_obs is None) :
      assert (f_obs is not None) and (f_obs.sigmas() is not None)
      i_obs = f_obs.f_as_f_sq()
    assert (i_obs.sigmas() is not None)
    if (wilson_b is not None) and (pdb_file is not None) :
      print "  Correcting reference data intensity falloff..."
      f_obs = i_obs.f_sq_as_f()
      pdb_hierarchy = any_file(pdb_file).file_object.construct_hierarchy()
      n_residues, n_bases = get_counts(pdb_hierarchy)
      iso_scale, aniso_scale = wilson_scaling(
        F=f_obs,
        n_residues=n_residues,
        n_bases=n_bases)
      # TODO anisotropic?
      print "  Scaling statistics for unmodified reference data:"
      show_b_factor_info(iso_scale, aniso_scale)
      delta_b = wilson_b - iso_scale.b_wilson
      f_obs = add_b_iso(f_obs, delta_b)
      i_obs = f_obs.f_as_f_sq()
    i_max = flex.max(i_obs.data())
    i_norm = i_obs.customized_copy(
      data=i_obs.data() / i_max,
      sigmas=i_obs.sigmas() / i_max)
    i_norm.setup_binner(n_bins=20)
    i_over_sigma = i_obs.data() / i_obs.sigmas()
    for i_bin in i_norm.binner().range_used() :
      sel = i_norm.binner().selection(i_bin)
      i_shell = i_norm.select(sel)
      sn_shell = i_over_sigma.select(sel)
      noise_bins = shell_intensity_bins(
        i_norm=i_shell,
        i_over_sigma=sn_shell,
        n_bins=n_intensity_bins)
      self._resolution_bins.append(noise_bins)

  def get_noise_profile_for_shell (self, i_bin) :
    return self._resolution_bins[i_bin-1]

class shell_intensity_bins (object) :
  def __init__ (self, i_norm, i_over_sigma, n_bins=20) :
    self._binner = bin_by_intensity(i_norm.data(), n_bins)
    self._sn_bins = []
    self._i_bins = []
    for n in range(n_bins) :
      bin_sel = self._binner.get_bin_selection(n)
      sn_bin = i_over_sigma.select(bin_sel)
      self._sn_bins.append(sn_bin)
      i_bin = i_norm.select(bin_sel)
      self._i_bins.append(i_bin)

  def get_i_over_sigma (self, I) :
    assert (I <= 1)
    if (I == 0) : return 1 # FIXME add French-Wilson before getting here
    k = self._binner.get_bin(I)
    sn_profile = self._sn_bins[k]
    for i_ref in self._i_bins[k] :
      if (i_ref >= I) :
        return sn_profile[k]
    return sn_profile[-1]

class bin_by_intensity (object) :
  def __init__ (self, data, n_bins=20) :
    from scitbx.array_family import flex
    self._bin_limits = []
    self._bins = flex.int(data.size(), -1)
    sorted_data = sorted(data)
    bin_size = data.size() // n_bins
    for n in range(n_bins) :
      bin_limit = bin_size * (n+1)
      if (bin_limit >= data.size()) :
        bin_limit = data.size() - 1
      i_max = sorted_data[bin_limit]
      self._bin_limits.append(i_max)
    for k, I in enumerate(data) :
      bin = self.get_bin(I)
      self._bins[k] = bin

  def get_bin (self, I) :
    assert (I >= 0)
    for k, limit in enumerate(self._bin_limits) :
      if (I < limit) :
        return k
    return len(self._bin_limits) - 1

  def get_bin_selection (self, bin) :
    from scitbx.array_family import flex
    sel = flex.bool(self._bins.size(), False)
    for k, n in enumerate(self._bins) :
      if (n == bin) :
        sel[k] = True
    return sel

if (__name__ == "__main__") :
  run(sys.argv[1:])
