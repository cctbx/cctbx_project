from cctbx.array_family import flex
import mmtbx.f_model
from mmtbx import utils
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx.pdb import crystal_symmetry_from_pdb
from cStringIO import StringIO
from cctbx.xray import ext
from mmtbx import utils
from scitbx.array_family import shared
from cctbx import maptbx
import iotbx.phil
from iotbx import crystal_symmetry_from_any
from cctbx import adptbx
from iotbx.option_parser import iotbx_option_parser
from libtbx.utils import Sorry
import os
from cctbx import miller


master_params_str = """\
grid_resolution_factor = 1./5
  .type = float
  .help = Map resolution factor: map grid step = d_min * resolution_factor, \
          d_min is highest resolution of the dataset.
map_smearing_b_factor = 0.0
  .type = float
  .expert_level = 2
  .help = Zero values requires very small grid step. Coarse grid will require \
          some positive value for map_smearing_b_factor.
map_1
  .help = First map to use in map CC calculation
{
  pdb_file_name = None
    .type = str
    .multiple = True
    .help = PDB file name.
  reflection_file_name = None
    .type = str
    .help = File with experimental data (most of formats: CNS, SHELX, MTZ, etc).
  data_labels = None
    .type = str
    .help = Labels for experimental data.
  flags_label = None
    .type = str
    .help = Free-R flags labels.
  test_flag_value = None
    .type = int
    .help = Reflection flag value used to define cross-validation set (free-R \
            flags).
  map_type = None
    .type = str
    .help = Electron density map type. Example xmFobs-yDFcalc (for \
            maximum-likelihood weighted map) or xFobs-yFcalc (for simple \
            unweighted map), x and y are any real numbers.
  use = None
    .type = bool
    .help = Use this model to compute local map CC.
  high_resolution = None
    .type = float
  low_resolution = None
    .type = float
}
map_2
  .help = Second map to use in map CC calculation
{
  pdb_file_name = None
    .type = str
    .multiple = True
    .help = PDB file name.
  reflection_file_name = None
    .type = str
    .help = File with experimental data (most of formats: CNS, SHELX, MTZ, etc).
  data_labels = None
    .type = str
    .help = Labels for experimental data.
  flags_label = None
    .type = str
    .help = Free-R flags labels.
  test_flag_value = None
    .type = int
    .help = Reflection flag value used to define cross-validation set (free-R \
            flags).
  map_type = None
    .type = str
    .help = Electron density map type. Example xmFobs-yDFcalc (for \
            maximum-likelihood weighted map) or xFobs-yFcalc (for simple \
            unweighted map), x and y are any real numbers.
  use = None
    .type = bool
    .help = Use this model to compute local map CC.
  high_resolution = None
    .type = float
  low_resolution = None
    .type = float
}
"""
master_params = iotbx.phil.parse(master_params_str, process_includes=False)

def sampled_density_map_obj(xray_structure, fft_map, b_base):
  assert not fft_map.anomalous_flag()
  xrs = xray_structure
  sampled_density = ext.sampled_model_density(
    unit_cell                             = xrs.unit_cell(),
    scatterers                            = xrs.scatterers(),
    scattering_type_registry              = xrs.scattering_type_registry(),
    fft_n_real                            = fft_map.real_map().focus(),
    fft_m_real                            = fft_map.real_map().all(),
    u_base                                = adptbx.b_as_u(b_base),
    wing_cutoff                           = 0.001,
    exp_table_one_over_step_size          = -100,
    force_complex                         = False,
    sampled_density_must_be_positive      = False,
    tolerance_positive_definite           = 1.e-5,
    use_u_base_as_u_extra                 = True,
    store_grid_indices_for_each_scatterer = -1)
  # model_map = sampled_density.real_map() # It is in P1 !!!
  return sampled_density

class model_to_map(object):
  def __init__(self, xray_structure, r_free_flags, f_obs, map_type,
                     resolution_factor, other_fft_map = None):
    if(r_free_flags is None):
      r_free_flags = f_obs.array(data = flex.bool(f_obs.size(), False))
    self.fmodel = mmtbx.f_model.manager(
      xray_structure = xray_structure,
      r_free_flags   = r_free_flags,
      target_name    = "ls_wunit_k1",
      f_obs          = f_obs)
    self.fmodel.remove_outliers()
    self.fmodel.update_solvent_and_scale()
    self.f_model_data = self.fmodel.f_model_scaled_with_k1().data()
    map_coeff = self.fmodel.map_coefficients(map_type = map_type, b_sharp=None)
    self.fft_map = self.fmodel.electron_density_map(
      resolution_factor = resolution_factor,
      map_type          = map_type,
      map_coefficients  = map_coeff,
      other_fft_map     = other_fft_map)
    self.fft_map.apply_sigma_scaling()
    self.map_coeff_data = map_coeff.data()

class pdb_to_xrs(object):
  def __init__(self, pdb_files, crystal_symmetry):
    processed_pdb_files_srv = utils.process_pdb_file_srv(
      crystal_symmetry = crystal_symmetry,
      log = StringIO())
    self.xray_structure = None
    self.processed_pdb_file = None
    if(pdb_files != []):
      self.processed_pdb_file, pdb_inp = \
        processed_pdb_files_srv.process_pdb_files(pdb_file_names = pdb_files)
      self.xray_structure = self.processed_pdb_file.xray_structure(
        show_summary = False)
      assert self.xray_structure is not None

def extract_crystal_symmetry(params):
  crystal_symmetries = []
  crystal_symmetry = None
  for cs_source in [params.map_1.pdb_file_name,
                    params.map_2.pdb_file_name]:
    if(cs_source is not None):
      for cs_source_ in cs_source:
        cs = None
        try:
          cs = crystal_symmetry_from_pdb.extract_from(file_name = cs_source_)
        except RuntimeError, e:
          if(str(e) == "No CRYST1 record."): pass
        if(cs is not None): crystal_symmetries.append(cs)
  for cs_source in [params.map_1.reflection_file_name,
                    params.map_2.reflection_file_name]:
    if(cs_source is not None):
      cs = crystal_symmetry_from_any.extract_from(cs_source)
      if(cs is not None): crystal_symmetries.append(cs)
  if(len(crystal_symmetries)==0):
    raise Sorry("No crystal symmetry is found.")
  elif(len(crystal_symmetries)>1):
    cs0 = crystal_symmetries[0]
    for cs in crystal_symmetries[1:]:
     if(not cs0.is_similar_symmetry(cs)):
        raise Sorry("Crystal symmetry mismatch between different files.")
    crystal_symmetry = crystal_symmetries[0]
  else:
    crystal_symmetry = crystal_symmetries[0]
  assert crystal_symmetry is not None
  return crystal_symmetry

def extract_data_and_flags(params, crystal_symmetry):
  data_and_flags = None
  if(params.reflection_file_name is not None):
    reflection_file = reflection_file_reader.any_reflection_file(
      file_name = params.reflection_file_name)
    reflection_file_server = reflection_file_utils.reflection_file_server(
      crystal_symmetry = crystal_symmetry,
      force_symmetry   = True,
      reflection_files = [reflection_file])
    parameters = utils.data_and_flags.extract()
    if(params.data_labels is not None):
      parameters.labels = [params.data_labels]
    if(params.flags_label is not None):
      parameters.r_free_flags.label = params.flags_label
    if(params.test_flag_value is not None):
      parameters.r_free_flags.test_flag_value = params.test_flag_value
    if(params.high_resolution is not None):
      parameters.high_resolution = params.high_resolution
    if(params.low_resolution is not None):
      parameters.low_resolution = params.low_resolution
    data_and_flags = utils.determine_data_and_flags(
      reflection_file_server = reflection_file_server,
      parameters             = parameters,
      data_description       = "X-ray data",
      log                    = StringIO())
  return data_and_flags

def compute_map_from_model(high_resolution, low_resolution, xray_structure,
                           grid_resolution_factor, crystal_gridding = None):
  # XXX if(params.scattering_table == "neutron"):
  # XXX   xray_structure.switch_to_neutron_scattering_dictionary()
  # XXX else:
  # XXX   xray_structure.scattering_type_registry(
  # XXX     table = params.scattering_table, d_min = hr)
  f_calc = xray_structure.structure_factors(d_min = high_resolution).f_calc()
  f_calc = f_calc.resolution_filter(d_max = low_resolution)
  if(crystal_gridding is None):
    return f_calc.fft_map(
      resolution_factor = grid_resolution_factor,
      symmetry_flags    = None)
  return miller.fft_map(
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = f_calc)

def cmd_run(args, command_name):
  msg = """\
Description: Compute local map correlation coefficient (default: residue-by-residue)

How to use:

1: Run this command: phenix.real_space_correlation;
2: Copy, save into a file and edit the parameters shown between the lines *** above;
3: Run the command with this parameters file.
"""
  if(len(args) == 0):
    args = ["--help"]
    print "*"*79
    master_params.show()
    print "*"*79
  command_line = (iotbx_option_parser(
    usage="\n%s parameters.txt" % command_name,
    description=msg)
    ).process(args=args)
  if(len(args) != 1):
    raise Sorry("One command line argument expected: parameters file.")
  if(len(args) == 1):
    arg = args[0]
    if(not os.path.isfile(arg)):
      raise Sorry("%s is not a file."%arg)
    parsed_params = None
    try: parsed_params = iotbx.phil.parse(file_name=arg)
    except KeyboardInterrupt: raise
    except RuntimeError, e:
      print e
    params, unused_definitions = master_params.fetch(
        sources = [parsed_params],
        track_unused_definitions = True)
    if(len(unused_definitions)):
      print "*"*79
      print "ERROR:",
      print "Unused parameter definitions:"
      for obj_loc in unused_definitions:
        print " ", str(obj_loc)
      print "*"*79
      raise Sorry("Fix parameters file and run again.")
      print
    run(params = params.extract())

def run(params, d_min_default=1.5, d_max_default=999.9):
  # check resolution factor
  if(params.grid_resolution_factor >= 0.5):
    raise Sorry("grid_resolution_factor must be < 0.5.")
  # check for crystal_symmetry
  crystal_symmetry = extract_crystal_symmetry(params = params)
  # read in the PDB files
  if(params.map_1.pdb_file_name+params.map_2.pdb_file_name == []):
    raise Sorry("No coordinate file given.")
  pdb_to_xrs_1 = pdb_to_xrs(pdb_files        = params.map_1.pdb_file_name,
                            crystal_symmetry = crystal_symmetry)
  xray_structure_1 = pdb_to_xrs_1.xray_structure
  pdb_to_xrs_2 = pdb_to_xrs(pdb_files        = params.map_2.pdb_file_name,
                            crystal_symmetry = crystal_symmetry)
  xray_structure_2 = pdb_to_xrs_2.xray_structure
  # assert correct combination of options
  if([params.map_1.use, params.map_2.use].count(True) == 0):
    raise Sorry("No model selected to compute local density correlation.")
  if([params.map_1.use, params.map_2.use].count(True) == 2):
    raise Sorry(
      "Select one model to compute local density correlation (two selected).")
  if(params.map_1.use):
    if(xray_structure_1 is None):
      raise Sorry("PDB file for map_1 is not provided.")
  if(params.map_2.use):
    if(xray_structure_2 is None):
      raise Sorry("PDB file for map_2 is not provided.")
  # read in F-obs and free-r flags for map_1
  data_and_flags_1 = extract_data_and_flags(
    params           = params.map_1,
    crystal_symmetry = crystal_symmetry)
  data_and_flags_2 = extract_data_and_flags(
    params           = params.map_2,
    crystal_symmetry = crystal_symmetry)
  if([data_and_flags_1, data_and_flags_2].count(None) != 2):
    # compute maps
    ### first
    if(xray_structure_1 is not None): xrs = xray_structure_1
    else: xrs = xray_structure_2
    if(data_and_flags_1 is not None): data_and_flags = data_and_flags_1
    else: data_and_flags = data_and_flags_2
    model_to_map_obj_1 = model_to_map(
      xray_structure    = xrs,
      r_free_flags      = data_and_flags.r_free_flags,
      f_obs             = data_and_flags.f_obs,
      map_type          = params.map_1.map_type,
      resolution_factor = params.grid_resolution_factor)
    fft_map_1 = model_to_map_obj_1.fft_map
    map_1 = fft_map_1.real_map()
    ### second
    if(xray_structure_2 is not None): xrs = xray_structure_2
    else: xrs = xray_structure_1
    if(data_and_flags_2 is not None): data_and_flags = data_and_flags_2
    else: data_and_flags = data_and_flags_1
    model_to_map_obj_2 = model_to_map(
      xray_structure    = xrs,
      r_free_flags      = data_and_flags.r_free_flags,
      f_obs             = data_and_flags.f_obs,
      map_type          = params.map_1.map_type,
      resolution_factor = params.grid_resolution_factor,
      other_fft_map     = fft_map_1)
    map_2 = model_to_map_obj_2.fft_map.real_map()
  # compute data if no reflection files provided for both maps
  if([data_and_flags_1, data_and_flags_2].count(None) == 2):
    if(xray_structure_1 is not None): xrs = xray_structure_1
    else: xrs = xray_structure_2
    if(params.map_1.high_resolution is not None):
      hr = params.map_1.high_resolution
    elif(params.map_2.high_resolution is not None):
      hr = params.map_2.high_resolution
    else:
      hr = d_min_default
    if(params.map_1.low_resolution is not None):
      lr = params.map_1.low_resolution
    elif(params.map_2.low_resolution is not None):
      lr = params.map_2.low_resolution
    else:
      lr = d_max_default
    fft_map_1 = compute_map_from_model(
      high_resolution        = hr,
      low_resolution         = lr,
      xray_structure         = xrs,
      grid_resolution_factor = params.grid_resolution_factor,
      crystal_gridding = None)
    map_1 = fft_map_1.real_map()
    if(xray_structure_2 is not None): xrs = xray_structure_2
    else: xrs = xray_structure_1
    map_2 = compute_map_from_model(
      high_resolution        = hr,
      low_resolution         = lr,
      xray_structure         = xrs,
      grid_resolution_factor = params.grid_resolution_factor,
      crystal_gridding       = fft_map_1).real_map()
  #
  assert map_1.focus() == map_2.focus()
  assert map_1.all() == map_2.all()
  # get sampled points
  if(params.map_1.use): xrs = xray_structure_1
  elif(params.map_2.use): xrs = xray_structure_2
  else: RuntimeError
  sampled_density = sampled_density_map_obj(
    xray_structure = xrs,
    fft_map        = fft_map_1,
    b_base         = params.map_smearing_b_factor)
  # prepare selections
  pdb_to_xrs_ = None
  if(params.map_1.use): pdb_to_xrs_ = pdb_to_xrs_1
  elif(params.map_2.use): pdb_to_xrs_ = pdb_to_xrs_2
  else: RuntimeError
  models = pdb_to_xrs_.processed_pdb_file.all_chain_proxies.\
    pdb_hierarchy.models()
  res_group_selections = []
  res_names = []
  res_ids = []
  for model in models:
    for chain in model.chains():
      for rg in chain.residue_groups():
        rg_i_seqs = []
        r_name = None
        for ag in rg.atom_groups():
          if(r_name is None): r_name = ag.resname
          for atom in ag.atoms():
            rg_i_seqs.append(atom.i_seq)
        if(len(rg_i_seqs) != 0):
          res_group_selections.append(flex.size_t(rg_i_seqs))
          res_names.append(r_name)
          res_ids.append(rg.resid())
  assert len(res_group_selections) == len(res_names)
  assert len(res_group_selections) == len(res_ids)
  # combine selections
  gifes = sampled_density.grid_indices_for_each_scatterer()
  residue_selections = shared.stl_set_unsigned()
  for i_seqs in res_group_selections:
    residue_selections.append_union_of_selected_arrays(
      arrays    = gifes,
      selection = i_seqs)
  # compute and output map CC
  result = flex.double()
  print
  print "Count, residue number, name and map CC"
  for i_count, res_sel in enumerate(residue_selections):
    corr = flex.linear_correlation(
      x = map_1.select(res_sel),
      y = map_2.select(res_sel)).coefficient()
    result.append(corr)
    print "%-5d %s %s %6.3f"%(i_count,res_ids[i_count],res_names[i_count],corr)
  print
  return result
