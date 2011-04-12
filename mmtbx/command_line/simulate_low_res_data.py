
import iotbx.phil
import libtbx.load_env
import libtbx.phil.command_line
from libtbx import easy_pickle
import random
import os
import sys

master_phil = iotbx.phil.parse("""
simulate_data {
  pdb_file = None
    .type = path
  hkl_file = None
    .type = path
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
    convert_to_isotropic = True
      .type = bool
    set_mean_b_iso = None
      .type = float
  }
  truncate {
    overall_b = None
      .type = float
    add_random_error_percent = None
      .type = float
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
  working_phil.show()
  params_ = working_phil.extract()
  params = params_.simulate_data
  if (params.pdb_file is None) and (params.hkl_file is None) :
    raise Sorry("No PDB file specified.")
  elif (params.pdb_file is not None) and (params.hkl_file is not None) :
    raise Sorry("Please use *either* a PDB or an MTZ file (not both).")
  if (pdb_in is None) and (params.pdb_file is not None) :
    f = file_reader.any_file(params.pdb_file, force_type="pdb")
    f.assert_file_type("pdb")
    pdb_in = f.file_object
  elif (hkl_in is None) and (params.hkl_file is not None) :
    f = file_reader.any_file(params.hkl_File, force_type="hkl")
    f.assert_file_type("hkl")
    hkl_in = f.file_object
  if (pdb_in is not None) :
    print "Generating fake data with phenix.fmodel..."
    F, r_free = from_pdb(pdb_in, params_)
  elif (hkl_in is not None) :
    print "Using experimental data..."
    F, r_free = from_hkl(hkl_in, params)
  print "Applying low-resolution filtering @ %.3f A..." % params.d_min
  F_out = truncate_data(F, params)
  mtz_dataset = F_out.as_mtz_dataset(
    column_root_label="F",
    column_types="F")
  mtz_dataset.add_miller_array(
    miller_array=r_free,
    column_root_label="FreeR_flag",
    column_types="I")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name=params.output_file)
  print "Wrote %s" % params.output_file

def from_pdb (pdb_in, params_) :
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
  xray_structure = pdb_in.xray_structure_simple(
    crystal_symmetry=apply_symm)
  if (params.modify_pdb.convert_to_isotropic) :
    xray_structure.convert_to_isotropic()
  if (params.modify_pdb.set_mean_b_iso is not None) :
    assert (params.modify_pdb.set_mean_b_iso > 0)
    u_iso = xray_structure.extract_u_iso_or_u_equiv()
    u_mean = flex.mean(u_iso)
    b_mean = adptbx.u_as_b(u_mean)
    scale = params.modify_pdb.set_mean_b_iso / b_mean
    xray_structure.set_u_iso(u_iso * scale)
  import mmtbx.command_line.fmodel
  from mmtbx import utils
  fmodel_params = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params.extract()
  fmodel_params.high_resolution = params.d_min
  fmodel_params.fmodel = params_.fake_data.fmodel
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
      print "F-obs: %s" % params.data_label
    elif (array.is_xray_intensity_array()) :
      f_obs = array.f_sq_as_f().map_to_asu()
      params.data_label = array.info().label_string()
      print "I-obs: %s" % params.data_label
    if (params.r_free_label is not None) :
      if (array.info().label_string() == params.r_free_label) :
        r_free = array.map_to_asu()
  f_obs = f_obs.resolution_filter(d_min=params.d_min)
  if (r_free is not None) :
    r_free = r_free.common_set(f_obs)
  return f_obs, r_free

def truncate_data (F, params) :
  from scitbx.array_family import flex
  if (params.truncate.overall_b is not None) :
    assert (params.truncate.overall_b > 0)
    d_min_sq = flex.pow2(F.d_spacings().data()) / 4.0
    scale = flex.exp(- params.truncate.overall_b / d_min_sq)
    F = F.customized_copy(data=F.data() * scale)
  if (params.truncate.add_random_error_percent) :
    assert (100 > params.truncate.add_random_error_percent > 0)
    data = F.data()
    fr = F.data()*params.truncate.add_random_error_percent/100.
    ri = flex.double()
    for trial in xrange(data.size()):
      r = random.randint(0,1)
      if(r == 0): r = -1
      ri.append(r)
    data = data + ri*fr
    F = F.array(data=data)
  return F

def get_mean_wilson_b_for_resolution (d_min, range=0.2) :
  from scitbx.array_family import flex
  pkl_file = libtbx.env.find_in_repositories(
    relative_path = "chem_data/polygon_data/all_mvd.pickle",
    test = os.path.isfile)
  db = easy_pickle.load(pkl_file)
  all_d_min = db['high_resolution']
  wilson_b = db['wilson_b']
  b_for_range = flex.double()
  for (d_, w_) in zip(all_d_min, wilson_b) :
    try :
      d = float(d_)
      w = float(w_)
    except ValueError : continue
    else :
      if (d > (d_min - range)) and (d < (d_min + range)) :
        b_for_range.append(w)
  h = flex.histogram(b_for_range, n_slots=10)
  print "Wilson B-factor statistics for d_min = %.3f - %.3f A" % (d_min-range,
    d_min+range)
  min = flex.min(b_for_range)
  max = flex.max(b_for_range)
  mean = flex.mean(b_for_range)
  print "  count: %d" % b_for_range.size()
  print "  min: %.2f" % min
  print "  max: %.2f" % max
  print "  mean: %.2f" % mean
  h.show(prefix="    ")
  return mean

if (__name__ == "__main__") :
  run(sys.argv[1:])
