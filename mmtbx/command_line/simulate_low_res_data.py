import iotbx.phil
import libtbx.load_env
from libtbx.str_utils import make_header
from libtbx.math_utils import ifloor
from libtbx.utils import Sorry
from libtbx import easy_pickle
from libtbx import adopt_init_args
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
  write_modified_pdb = True
    .type = bool
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
    missing_flags = *extend discard
      .type = choice(multi=False)
      .help = Handling of reflections for which R-free flags are not present.
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
    remove_alt_confs = False
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
    apply_b_to_sigmas = True
      .type = bool
      .help = If True, the B-factor scaling will be done on experimental \
        sigmas (if present) as well as amplitudes.
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
      .type = ints(size=3)
      .help = Axis of rotation (as Miller indices, i.e. three integers) of \
        crystal during data collection.  Only used if remove_cone_around_axis \
        is defined.
    elliptical_truncation = False
      .type = bool
      .help = Truncate data anisotropically, using the overall anisotropic \
        B-factor plus ellipse_scale to set the ellipse dimensions.
    ellipse_scale = 1.0
      .type = float
    ellipse_target_completeness = None
      .type = float
      .help = Target completeness of data (as percent, max=100) after \
        elliptical truncation
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

def run (args, out=None) :
  if (out is None) :
    out = sys.stdout
  make_header("mmtbx.simulate_low_res_data", out=out)
  print >> out, """
  For generation of realistic data (model-based, or using real
  high-resolution data) for methods development.

*********************************** WARNING: ***********************************
 this is an experimental program - definitely NOT bug-free.
                  Use at your own risk!

  Usage:
   mmtbx.simulate_low_res_data model.pdb [options...]
     (generate data from a PDB file)

   mmtbx.simulate_low_res_data highres.mtz [model.pdb] [options...]
     (truncate high-resolution data)

   mmtbx.simulate_low_res_data --help
     (print full parameters with additional info)
"""
  if (len(args) == 0) or ("--help" in args) :
    print >> out, "# full parameters:"
    if ("--help" in args) :
      master_phil.show(attributes_level=1)
    else :
      master_phil.show()
    return
  from iotbx import file_reader
  interpreter = master_phil.command_line_argument_interpreter(
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
        print >> out, "ignoring uninterpretable argument '%s'" % arg
      else :
        user_phil.append(arg_phil)
  working_phil = master_phil.fetch(sources=user_phil)
  make_header("Working parameters", out=out)
  working_phil.show(prefix="  ")
  params_ = working_phil.extract()
  params = params_.simulate_data
  prepare_data(
    params=params,
    hkl_in=hkl_in,
    pdb_in=pdb_in,
    out=out)

class prepare_data (object) :
  def __init__ (self, params, hkl_in=None, pdb_in=None, out=sys.stdout) :
    adopt_init_args(self, locals())
    self.params = params
    self.out = out
    if (params.pdb_file is None) and (params.hkl_file is None) :
      raise Sorry("No PDB file specified.")
    if (params.generate_noise.add_noise) and (params.hkl_file is None) :
      if (params.generate_noise.noise_profile_file is None) :
        raise Sorry("noise_profile_file required when add_noise=True and "
          "hkl_file is undefined.")
    if (pdb_in is None) and (params.pdb_file is not None) :
      f = file_reader.any_file(params.pdb_file, force_type="pdb")
      f.assert_file_type("pdb")
      self.pdb_in = f.file_object
    if (self.hkl_in is None) and (params.hkl_file is not None) :
      f = file_reader.any_file(params.hkl_File, force_type="hkl")
      f.assert_file_type("hkl")
      self.hkl_in = f.file_object
    if (self.pdb_in is not None) :
      self.pdb_hierarchy = self.pdb_in.construct_hierarchy()
    if (self.hkl_in is not None) :
      make_header("Extracting experimental data", out=sys.stdout)
      f_raw, r_free = self.from_hkl(self.hkl_in)
    elif (self.pdb_in is not None) :
      make_header("Generating fake data with phenix.fmodel", out=sys.stdout)
      f_raw, r_free = self.from_pdb()
    if (params.r_free_flags.file_name is not None) :
      f_raw, r_free = self.import_r_free_flags(f_raw)
    self.r_free = r_free
    make_header("Applying low-resolution filtering", out=sys.stdout)
    print >> out, "  Target resolution: %.2f A" % params.d_min
    self.n_residues, self.n_bases = None, None
    if (self.pdb_in is not None) :
      self.n_residues, self.n_bases = get_counts(self.pdb_hierarchy)
    #if (params.auto_adjust) :
    #  if (pdb_in is None) :
    #    raise Sorry("You must supply a PDB file when auto_adjust=True.")
    self.f_out = self.truncate_data(f_raw)
    if (params.generate_noise.add_noise) :
      make_header("Adding noise using sigma profile", out=sys.stdout)
      if (self.f_out.sigmas() is None) :
        if (self.pdb_in is not None) :
          iso_scale, aniso_scale = wilson_scaling(self.f_out, self.n_residues,
            self.n_bases)
        i_obs = create_sigmas(
          f_obs=self.f_out,
          params=params.generate_noise,
          wilson_b=iso_scale.b_wilson,
          return_as_amplitudes=False)
      apply_sigma_noise(i_obs)
      self.f_out = i_obs.f_sq_as_f()
    make_header("Done processing", out=sys.stdout)
    print >> out, "  Completeness after processing: %.2f%%" % (
      self.f_out.completeness() * 100.)
    print >> out, "  Final resolution: %.2f A" % self.f_out.d_min()
    if (self.pdb_in is not None) :
      iso_scale, aniso_scale = wilson_scaling(self.f_out, self.n_residues,
        self.n_bases)
      print >> out, ""
      print >> out, "  Scaling statistics for output data:"
      show_b_factor_info(iso_scale, aniso_scale, out=out)
      print >> out, ""
    self.write_output()

  def write_output (self) :
    f_out = self.f_out
    params = self.params
    out = self.out
    if (f_out.sigmas() is not None) :
      mtz_dataset = f_out.as_mtz_dataset(
        column_root_label="F",
        column_types="FQ")
    else :
      mtz_dataset = f_out.as_mtz_dataset(
        column_root_label="F",
        column_types="F")
    if (self.r_free is not None) :
      r_free = self.r_free.common_set(f_out)
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
    print >> out, "  Wrote %s" % params.output_file
    if (self.pdb_hierarchy is not None) and (params.write_modified_pdb) :
      pdb_out = os.path.splitext(params.output_file)[0] + ".pdb"
      f = open(pdb_out, "w")
      f.write("%s\n" % "\n".join(self.pdb_in.crystallographic_section()))
      f.write(self.pdb_hierarchy.as_pdb_string())
      f.close()
      print >> out, "  Wrote modified model to %s" % pdb_out

  def from_pdb (self) :
    out = self.out
    params = self.params
    pdb_hierarchy = self.pdb_hierarchy
    pdb_sg, pdb_uc = None, None
    pdb_symm = self.pdb_in.crystal_symmetry()
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
      print >> out, "  Removing solvent atoms..."
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
      print >> out, "  Removing all alternate conformations and resetting occupancies..."
      from mmtbx import pdbtools
      pdbtools.remove_alt_confs(hierarchy=pdb_hierarchy)
    xray_structure = self.pdb_in.xray_structure_simple(
      crystal_symmetry=apply_symm)
    sctr_keys = xray_structure.scattering_type_registry().type_count_dict().keys()
    if (not (("H" in sctr_keys) or ("D" in sctr_keys))) :
      print >> out, "  WARNING: this model does not contain hydrogen atoms!"
      print >> out, "           strongly recommend running phenix.ready_set or"
      print >> out, "           equivalent to ensure realistic simulated data."
      print >> out, ""
    if (params.modify_pdb.convert_to_isotropic) :
      xray_structure.convert_to_isotropic()
    set_b = None
    if (params.modify_pdb.set_mean_b_iso is not None) :
      assert (not params.modify_pdb.set_wilson_b)
      print >> out, "  Scaling B-factors to have mean of %.2f" % \
        params.modify_pdb.set_mean_b_iso
      assert (params.modify_pdb.set_mean_b_iso > 0)
      set_b = params.modify_pdb.set_mean_b_iso
    elif (params.modify_pdb.set_wilson_b) :
      print >> out, "  Scaling B-factors to match mean Wilson B for this resolution"
      set_b = get_mean_statistic_for_resolution(
        d_min=params.d_min,
        stat_type="wilson_b")
      print >> out, ""
    if (set_b is not None) :
      u_iso = xray_structure.extract_u_iso_or_u_equiv()
      u_mean = flex.mean(u_iso)
      b_mean = adptbx.u_as_b(u_mean)
      scale = set_b / b_mean
      xray_structure.scale_adps(scale)
      pdb_hierarchy.atoms().set_adps_from_scatterers(
        scatterers=xray_structure.scatterers(),
        unit_cell=xray_structure.unit_cell())
    import mmtbx.command_line.fmodel
    from mmtbx import utils
    fmodel_params = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params.extract()
    fmodel_params.high_resolution = params.d_min
    fake_data = params.fake_data_from_fmodel
    fmodel_params.fmodel = fake_data.fmodel
    if (fmodel_params.fmodel.b_sol == 0) :
      print >> out, "  b_sol is zero - will use mean value for d_min +/- 0.2A"
      print >> out, "   (this is not strongly correlated with resolution, but"
      print >> out, "    it is preferrable to use a real value instead of leaving"
      print >> out, "    it set to 0)"
      fmodel_params.fmodel.b_sol = get_mean_statistic_for_resolution(
        d_min=params.d_min,
        stat_type="b_sol")
      print >> out, ""
    if (fmodel_params.fmodel.k_sol == 0) :
      print >> out, "  k_sol is zero - will use mean value for d_min +/- 0.2A"
      print >> out, "   (this is not strongly correlated with resolution, but"
      print >> out, "    it is preferrable to use a real value instead of leaving"
      print >> out, "     it set to 0)"
      fmodel_params.fmodel.k_sol = get_mean_statistic_for_resolution(
        d_min=params.d_min,
        stat_type="k_sol")
      print >> out, ""
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

  def from_hkl (self) :
    params = self.params
    out = self.out
    miller_arrays = self.hkl_in.as_miller_arrays()
    f_obs = None
    r_free = None
    for array in miller_arrays :
      if (params.data_label is not None) :
        if (array.info().label_string() == params.data_label) :
          f_obs = array.map_to_asu()
      elif (array.is_xray_amplitude_array()) :
        f_obs = array.map_to_asu()
        params.data_label = array.info().label_string()
        print >> out, "  F-obs: %s" % params.data_label
      elif (array.is_xray_intensity_array()) :
        f_obs = array.f_sq_as_f().map_to_asu()
        params.data_label = array.info().label_string()
        print >> out, "  I-obs: %s" % params.data_label
      if (params.r_free_flags.label is not None) :
        if (array.info().label_string() == params.r_free_flags.label) :
          r_free = array.map_to_asu()
          print >> out, "  R-free: %s" % array.info().label_string()
      elif (array.info().label_string() in ["FreeR_flag", "FREE"]) :
        r_free = array.map_to_asu()
        print >> out, "  R-free: %s" % array.info().label_string()
    f_obs = f_obs.resolution_filter(d_min=params.d_min)
    if (r_free is not None) :
      r_free = r_free.common_set(f_obs)
    return f_obs, r_free

  def import_r_free_flags (self, F) :
    params = self.params.r_free_flags
    out = self.out
    from iotbx import file_reader
    rfree_in = file_reader.any_file(params.file_name)
    rfree_in.assert_file_type("hkl")
    hkl_server = rfree_in.file_server
    r_free_raw, flag_value = hkl_server.get_r_free_flags(
      file_name=None,
      label=params.label,
      test_flag_value=None,
      parameter_scope="simulate_data.r_free_flags",
      disable_suitability_test=False)
    r_free = r_free_raw.customized_copy(data=r_free_raw.data() == flag_value)
    r_free = r_free.map_to_asu().common_set(F)
    print >> out, "  Using R-free flags from %s:%s" % (rfree_in.file_name,
      r_free_raw.info().label_string())
    if (F.data().size() != r_free.data().size()) :
      n_missing = F.data().size() - r_free.data().size()
      assert (n_missing > 0)
      if (params.missing_flags == "discard") :
        print >> out, "    discarding %d amplitudes without R-free flags" % \
          n_missing
        F = F.common_set(r_free)
      else :
        print >> out, "    generating missing R-free flags for %d reflections" %\
          n_missing
        missing_set = F.lone_set(r_free)
        missing_flags = missing_set.generate_r_free_flags(
          fraction=r_free.data().count(True) / r_free.data().size(),
          max_free=None,
          use_lattice_symmetry=True)
        r_free = r_free.concatenate(other=missing_flags)
    assert (F.data().size() == r_free.data().size())
    return F, r_free

  def truncate_data (self, F) :
    params = self.params
    out = self.out
    from scitbx.array_family import flex
    if (params.truncate.add_b_iso is not None) :
      print >> out, "  Applying isotropic B-factor of %.2f A^2" % (
        params.truncate.add_b_iso)
      F = F.apply_debye_waller_factors(
        b_iso=params.truncate.add_b_iso,
        apply_to_sigmas=params.truncate.apply_b_to_sigmas)
    if (params.truncate.add_b_aniso != [0,0,0,0,0,0]) :
      print >> out, "  Adding anisotropy..."
      F = F.apply_debye_waller_factors(
        b_cart=params.truncate.add_b_aniso,
        apply_to_sigmas=params.truncate.apply_b_to_sigmas)
    if (params.truncate.add_random_error_percent is not None) :
      print >> out, "  Adding random error as percent of amplitude..."
      F = add_random_error(f_obs=F,
        error_percent=params.truncate.add_random_error_percent)
    if (params.truncate.remove_cone_around_axis is not None) :
      print >> out, "  Removing cone of data around axis of rotation..."
      print >> out, "    radius = %.1f degrees" % \
        params.truncate.remove_cone_around_axis
      assert (params.truncate.axis_of_rotation is not None)
      n_hkl, delta_n_hkl = remove_cone_around_axis(
        array=F,
        axis=params.truncate.axis_of_rotation,
        cone_radius=params.truncate.remove_cone_around_axis)
      print >> out, "    removed %d reflections (out of %d)" %(delta_n_hkl,
        n_hkl)
    if (params.truncate.elliptical_truncation) :
      print >> out, "  Truncating the data elliptically..."
      target_completeness = params.truncate.ellipse_target_completeness
      completeness_start = F.completeness() * 100.0
      if (completeness_start < target_completeness) :
        print >> out, "    completeness is already less than target value:"
        print >> out, "       %.2f versus %.2f"
        print >> out, "    elliptical truncation will be skipped."
      else :
        print >> out, "    using overall anisotropic B as a guide."
        iso_scale, aniso_scale = wilson_scaling(
          F=F,
          n_residues=self.n_residues,
          n_bases=self.n_bases)
        n_hkl, delta_n_hkl = elliptical_truncation(
          array=F,
          b_cart=aniso_scale.b_cart,
          scale_factor=params.truncate.ellipse_scale,
          target_completeness=target_completeness)
        print >> out, "    removed %d reflections (out of %d)" %(delta_n_hkl,
          n_hkl)
    return F

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

def show_b_factor_info (iso_scale, aniso_scale, out) :
  b_cart = aniso_scale.b_cart
  print >> out, "    overall isotropic B-factor:   %6.2f" % iso_scale.b_wilson
  print >> out, "    overall anisotropic B-factor: %6.2f, %6.2f, %6.2f" % (
    b_cart[0], b_cart[3], b_cart[4])
  print >> out, "                                  %14.2f, %6.2f" % (
    b_cart[1], b_cart[5])
  print >> out, "                                  %22.2f" % b_cart[2]

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

def elliptical_truncation (array,
                           b_cart,
                           scale_factor=1.0,
                           target_completeness=None) :
  from cctbx import adptbx
  from scitbx.array_family import flex
  indices = array.indices()
  axis_index = -1
  min_b_directional = sys.maxint
  for n, b_index in enumerate(b_cart[0:3]) :
    if (b_index < min_b_directional) :
      min_b_directional = b_index
      axis_index = n
  assert (0 <= axis_index <= 3)
  max_index_along_axis = [0,0,0]
  for hkl in indices :
    if (hkl[axis_index] > max_index_along_axis[axis_index]) :
      max_index_along_axis[axis_index] = hkl[axis_index]
  assert (max_index_along_axis != [0,0,0])
  u_star = adptbx.u_cart_as_u_star(array.unit_cell(), adptbx.b_as_u(b_cart))
  scale_cutoff = adptbx.debye_waller_factor_u_star(
    h=max_index_along_axis, u_star=u_star) * scale_factor
  scale = array.debye_waller_factors(b_cart=b_cart).data()
  if (target_completeness is not None) :
    assert (target_completeness > 0) and (target_completeness <= 100)
    completeness_start = array.completeness() * 100.0
    assert (completeness_start > target_completeness)
    scale_srt = sorted(scale)
    i_max = ifloor(len(scale_srt)*(completeness_start-target_completeness)/100)
    scale_cutoff = scale_srt[i_max]
  data = array.data()
  sigmas = array.sigmas()
  n_hkl = indices.size()
  i = 0
  while (i < len(indices)) :
    if (scale[i] < scale_cutoff) :
      del indices[i]
      del data[i]
      del scale[i]
      if (sigmas is not None) :
        del sigmas[i]
    else :
      i += 1
  delta_n_hkl = n_hkl - indices.size()
  return (n_hkl, delta_n_hkl)

def remove_cone_around_axis (array, axis, cone_radius) :
  assert (cone_radius < 90) and (cone_radius >= 0)
  from scitbx.matrix import rec
  unit_cell = array.unit_cell()
  indices = array.indices()
  rc_vectors = unit_cell.reciprocal_space_vector(indices)
  v1 = rec(unit_cell.reciprocal_space_vector(tuple(axis)), (3,1))
  n_hkl = indices.size()
  data = array.data()
  sigmas = array.sigmas()
  angle_high = 180 - cone_radius
  i = 0
  while (i < len(indices)) :
    v2 = rec(rc_vectors[i], (3,1))
    #print v1, v2
    angle = math.degrees(v1.angle(v2))
    if (angle <= cone_radius) or (angle >= angle_high) :
      del rc_vectors[i]
      del indices[i]
      del data[i]
      if (sigmas is not None) :
        del sigmas[i]
    else :
      i += 1
  delta_n_hkl = n_hkl - indices.size()
  return (n_hkl, delta_n_hkl)

stat_names = {
  'wilson_b' : 'Wilson B-factor',
  'k_sol' : 'Bulk solvent scale factor (k_sol)',
  'b_sol' : 'Bulk solvent B-factor (b_sol)',
}

def get_mean_statistic_for_resolution (d_min, stat_type, range=0.2, out=None) :
  if (out is None) :
    out = sys.stdout
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
  print >> out, "  %s for d_min = %.3f - %.3f A" % (stat_names[stat_type], d_min-range,
    d_min+range)
  min = flex.min(values_for_range)
  max = flex.max(values_for_range)
  mean = flex.mean(values_for_range)
  print >> out, "    count: %d" % values_for_range.size()
  print >> out, "    min: %.2f" % min
  print >> out, "    max: %.2f" % max
  print >> out, "    mean: %.2f" % mean
  print >> out, "    histogram of values:"
  h.show(prefix="      ")
  return mean

def create_sigmas (f_obs, params, wilson_b=None, return_as_amplitudes=False) :
  assert (f_obs.sigmas() is None)
  from scitbx.array_family import flex
  i_obs = f_obs.f_as_f_sq()
  i_norm = i_obs.data() / flex.mean(i_obs.data())
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
                n_intensity_bins=20,
                out=None) :
    if (out is None) :
      out = sys.stdout
    if (wilson_b is None) or (pdb_file is None) :
      print >> out, """\
  WARNING: missing desired Wilson B-factor and/or PDB file
           for noise profile data.  Without this information
           the intensity falloff with resolution will probably
           not be the same for your synthetic data and the
           data used to generate sigmas.
"""
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
      print >> out, "  Correcting reference data intensity falloff..."
      f_obs = i_obs.f_sq_as_f()
      pdb_hierarchy = any_file(pdb_file).file_object.construct_hierarchy()
      n_residues, n_bases = get_counts(pdb_hierarchy)
      iso_scale, aniso_scale = wilson_scaling(
        F=f_obs,
        n_residues=n_residues,
        n_bases=n_bases)
      # TODO anisotropic?
      print >> out, "  Scaling statistics for unmodified reference data:"
      show_b_factor_info(iso_scale, aniso_scale, out=out)
      delta_b = wilson_b - iso_scale.b_wilson
      f_obs = f_obs.apply_debye_waller_factors(b_iso=delta_b)
      i_obs = f_obs.f_as_f_sq()
    i_mean = flex.max(i_obs.data())
    i_norm = i_obs.customized_copy(
      data=i_obs.data() / i_mean,
      sigmas=i_obs.sigmas() / i_mean)
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
      i_bin = i_norm.select(bin_sel)
      assert (sn_bin.size() == i_bin.size())
      self._sn_bins.append(sn_bin)
      self._i_bins.append(i_bin)

  def get_i_over_sigma (self, I) :
    assert (I <= 1)
    if (I == 0) : return 1 # FIXME add French-Wilson before getting here
    k = self._binner.get_bin(I)
    sn_profile = self._sn_bins[k]
    for j, i_ref in enumerate(self._i_bins[k]) :
      if (i_ref >= I) :
        return sn_profile[j]
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
