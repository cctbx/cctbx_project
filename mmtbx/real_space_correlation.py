from cctbx.array_family import flex
import mmtbx.utils
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx.pdb import crystal_symmetry_from_pdb
from cStringIO import StringIO
from cctbx import maptbx
import iotbx.phil
from iotbx import crystal_symmetry_from_any
from cctbx import adptbx
from libtbx.utils import Sorry
import os, math
from cctbx import miller
from mmtbx import map_tools
from libtbx import group_args
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
import sys
from cctbx import crystal
import iotbx.pdb
from iotbx.pdb import combine_unique_pdb_files

core_params_str = """\
atom_radius = None
  .type = float
  .help = Atomic radius for map CC calculation. Determined automatically if \
          if None is given
  .expert_level = 2
hydrogen_atom_radius = 1.0
  .type = float
  .help = Atomic radius for map CC calculation for H or D.
  .expert_level = 2
number_of_grid_points = 50
  .type = int
  .help = Requesteed number of grid points to be used in CC calculation. \
          Together with atom_radius it will define the grid step.
  .expert_level = 2
set_cc_to_zero_if_n_grid_points_less_than = 10
  .type = int
  .help = Return zero CC if number of grid nodes is less than defined above
  .expert_level = 2
poor_cc_threshold = 0.7
  .type = float
  .help = Ad hoc definition of poor map CC
poor_map_value_threshold = 1.0
  .type = float
  .help = Ad hoc value (less than) defining a week density (in sigma) for maps \
          involved in map CC calculation
  .expert_level = 2
crystal_symmetry
  .style = noauto
{
  unit_cell=None
    .type=unit_cell
  space_group=None
    .type=space_group
}
"""

master_params_str = """\
%s
scattering_table = *n_gaussian wk1995 it1992 neutron
  .type = choice
  .help = Scattering table for structure factors calculations
details_level = atom residue *automatic
  .type = choice(multi=False)
  .help = Level of details to show CC for
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
  use_kick_map = False
    .type = bool
    .expert_level = 2
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
  use_kick_map = False
    .type = bool
    .expert_level = 2
}
"""%core_params_str

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)


def compute_grid_step_from_atom_radius_and_number_of_grid_points(
      r, n, d_min, resolution_factor=1./3):
  assert n > 0 and r < 1000 and r > 0
  result = min(r*(4*math.pi/(3*n))**(1./3), d_min*resolution_factor)
  return result

class model_to_map(object):
  def __init__(self, xray_structure, f_obs, map_type,
                     resolution_factor, high_resolution, low_resolution,
                     use_kick_map,
                     other_fft_map = None):
    # XXX
    map_name_obj = mmtbx.map_names(map_name_string = map_type)
    if(map_name_obj.k == 0 and map_name_obj.n == -1 and not map_name_obj.ml_map):
      complete_set = f_obs.complete_set(d_min = f_obs.d_min(), d_max=None)
      assert not use_kick_map
      f_calc = complete_set.structure_factors_from_scatterers(
        xray_structure = xray_structure).f_calc()
      if(other_fft_map is not None):
        self.fft_map = miller.fft_map(
          crystal_gridding     = other_fft_map,
          fourier_coefficients = f_calc)
      else:
        self.fft_map = f_calc.fft_map(
          resolution_factor = min(0.5,resolution_factor),
          symmetry_flags    = None)
      self.fft_map.apply_sigma_scaling()
      self.map_data =  self.fft_map.real_map_unpadded()
    # XXX
    else:
      r_free_flags = f_obs.array(data = flex.bool(f_obs.size(), False))
      if(high_resolution is not None):
        f_obs = f_obs.resolution_filter(d_min = high_resolution)
        r_free_flags = r_free_flags.resolution_filter(d_min = high_resolution)
      if(low_resolution is not None):
        f_obs = f_obs.resolution_filter(d_max = low_resolution)
        r_free_flags = r_free_flags.resolution_filter(d_max = low_resolution)
      #
      bss_params = bss.master_params.extract()
      bss_params.k_sol_max = 0.6
      bss_params.k_sol_min = 0.0
      bss_params.b_sol_max = 500.0
      bss_params.b_sol_min = 0.0
      bss_params.k_sol_grid_search_max = 0.6
      bss_params.k_sol_grid_search_min = 0.0
      bss_params.b_sol_grid_search_max = 80.0
      bss_params.b_sol_grid_search_min = 20.0
      bss_params.k_sol_step = 0.3
      bss_params.b_sol_step = 20.0
      #
      self.fmodel = mmtbx.utils.fmodel_simple(
        xray_structures = [xray_structure],
        f_obs           = f_obs,
        r_free_flags    = r_free_flags,
        skip_twin_detection = True, # XXX remove once Peter supports map type strings
        bss_params      = bss_params)
      if(not use_kick_map):
        map_obj = self.fmodel.electron_density_map()
        map_coeff = map_obj.map_coefficients(map_type = map_type)
        self.fft_map = map_obj.fft_map(
          resolution_factor = min(0.5,resolution_factor),
          map_coefficients  = map_coeff,
          other_fft_map     = other_fft_map)
        self.fft_map.apply_sigma_scaling()
        self.map_data =  self.fft_map.real_map_unpadded()
      else:
        km = map_tools.kick_map(
          fmodel                        = self.fmodel,
          map_type                      = map_type,
          resolution_factor             = min(0.5,resolution_factor),
          other_fft_map                 = other_fft_map,
          real_map                      = True,
          real_map_unpadded             = False,
          symmetry_flags                = maptbx.use_space_group_symmetry,
          average_maps                  = True) # XXX use map coefficients averaging
        self.fft_map = km.fft_map
        self.map_data = km.map_data

class pdb_to_xrs(object):
  def __init__(self, pdb_files, crystal_symmetry, scattering_table):
    self.xray_structure = None
    if(pdb_files != []):
      pdb_combined = combine_unique_pdb_files(file_names = pdb_files)
      raw_recs = flex.std_string()
      for rec in pdb_combined.raw_records:
        if(rec.upper().count("CRYST1")==0):
          raw_recs.append(rec)
      raw_recs.append(iotbx.pdb.format_cryst1_record(
        crystal_symmetry = crystal_symmetry))
      pdb_inp = iotbx.pdb.input(source_info = None, lines = raw_recs)
      self.xray_structure = pdb_inp.xray_structure_simple()
      self.pdb_hierarchy = pdb_inp.construct_hierarchy()
      self.pdb_hierarchy.atoms().reset_i_seq() # VERY important to do.
      mmtbx.utils.setup_scattering_dictionaries(
        scattering_table = scattering_table,
        xray_structure = self.xray_structure,
        d_min = None)

def extract_crystal_symmetry(params):
  crystal_symmetries = []
  crystal_symmetry = None
  #
  if([params.crystal_symmetry.unit_cell,
      params.crystal_symmetry.space_group].count(None)==0):
    return crystal.symmetry(
        unit_cell=params.crystal_symmetry.unit_cell,
        space_group_info=params.crystal_symmetry.space_group)
  #
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
    parameters = mmtbx.utils.data_and_flags_master_params().extract()
    parameters.force_anomalous_flag_to_be_equal_to = False
    if(params.data_labels is not None):
      parameters.labels = [params.data_labels]
    if(params.high_resolution is not None):
      parameters.high_resolution = params.high_resolution
    if(params.low_resolution is not None):
      parameters.low_resolution = params.low_resolution
    data_and_flags = mmtbx.utils.determine_data_and_flags(
      reflection_file_server = reflection_file_server,
      parameters             = parameters,
      data_description       = "X-ray data",
      extract_r_free_flags   = False,
      log                    = StringIO())
  return data_and_flags

def compute_map_from_model(high_resolution, low_resolution, xray_structure,
                           grid_resolution_factor, crystal_gridding = None):
  f_calc = xray_structure.structure_factors(d_min = high_resolution).f_calc()
  f_calc = f_calc.resolution_filter(d_max = low_resolution)
  if(crystal_gridding is None):
    return f_calc.fft_map(
      resolution_factor = min(0.5,grid_resolution_factor),
      symmetry_flags    = None)
  return miller.fft_map(
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = f_calc)

def cmd_run(args, command_name):
  msg = """\

Tool to compute local map correlation coefficient.

How to use:
1: Run this command: phenix.real_space_correlation;
2: Copy, save into a file and edit the parameters shown between the lines *** below;
3: Run the command with this parameters file:
   phenix.real_space_correlation parameters.txt
"""
  if(len(args) == 0):
    print msg
    print "*"*79
    master_params().show()
    print "*"*79
    return
  else :
    arg = args[-1]
    if(not os.path.isfile(arg)):
      raise Sorry("%s is not a file."%arg)
    parsed_params = None
    try: parsed_params = iotbx.phil.parse(file_name=arg)
    except KeyboardInterrupt: raise
    except RuntimeError, e:
      print e
    params, unused_definitions = master_params().fetch(
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

def run(params, d_min_default=1.5, d_max_default=999.9) :
  # check for crystal_symmetry
  crystal_symmetry = extract_crystal_symmetry(params = params)
  # read in the PDB files
  if(params.map_1.pdb_file_name+params.map_2.pdb_file_name == []):
    raise Sorry("No coordinate file given.")
  pdb_to_xrs_1 = pdb_to_xrs(pdb_files        = params.map_1.pdb_file_name,
                            crystal_symmetry = crystal_symmetry,
                            scattering_table = params.scattering_table)
  xray_structure_1 = pdb_to_xrs_1.xray_structure
  pdb_to_xrs_2 = pdb_to_xrs(pdb_files        = params.map_2.pdb_file_name,
                            crystal_symmetry = crystal_symmetry,
                            scattering_table = params.scattering_table)
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
  #
  # get map CC object
  def get_map_cc_obj(map_1, params, pdb_to_xrs_1, pdb_to_xrs_2, fft_map_1,
                     atom_detail, residue_detail, atom_radius,
                     hydrogen_atom_radius):
    pdb_to_xrs_ = None
    if(params.map_1.use): pdb_to_xrs_ = pdb_to_xrs_1
    elif(params.map_2.use): pdb_to_xrs_ = pdb_to_xrs_2
    else: RuntimeError
    pdb_hierarchy = pdb_to_xrs_.pdb_hierarchy
    result = map_cc_funct(
      map_1          = map_1,
      map_1_name     = params.map_1.map_type,
      xray_structure = pdb_to_xrs_.xray_structure,
      fft_map        = fft_map_1,
      pdb_hierarchy  = pdb_hierarchy,
      atom_detail    = atom_detail,
      atom_radius    = atom_radius,
      hydrogen_atom_radius = hydrogen_atom_radius,
      residue_detail = residue_detail)
    del map_1
    return result
  #
  if([data_and_flags_1, data_and_flags_2].count(None) != 2):
    # compute maps
    if([params.map_1.map_type, params.map_2.map_type].count(None) == 2):
      raise Sorry("map_type has to be defined; example: "\
                   "map_type=2Fo-Fc or map_type=0Fo--1Fc (just Fc map).")
    if(params.map_1.map_type is None):
      params.map_1.map_type = params.map_2.map_type
    if(params.map_2.map_type is None):
      params.map_2.map_type = params.map_1.map_type
    ### first
    if(xray_structure_1 is not None): xrs = xray_structure_1
    else: xrs = xray_structure_2
    if(data_and_flags_1 is not None): data_and_flags = data_and_flags_1
    else: data_and_flags = data_and_flags_2
    hr, lr = None, None
    if(params.map_1.high_resolution is not None):
      hr = params.map_1.high_resolution
    elif(params.map_2.high_resolution is not None):
      hr = params.map_2.high_resolution
    if(params.map_1.low_resolution is not None):
      lr = params.map_1.low_resolution
    elif(params.map_2.low_resolution is not None):
      lr = params.map_2.low_resolution
    #
    atom_detail, residue_detail, atom_radius = set_details_level_and_radius(
      details_level = params.details_level,
      d_min         = data_and_flags.f_obs.d_min(),
      atom_radius   = params.atom_radius)
    grid_step = compute_grid_step_from_atom_radius_and_number_of_grid_points(
      r = atom_radius, n = params.number_of_grid_points,
      d_min = data_and_flags.f_obs.d_min())
    #
    model_to_map_obj_1 = model_to_map(
      xray_structure    = xrs,
      f_obs             = data_and_flags.f_obs,
      map_type          = params.map_1.map_type,
      resolution_factor = grid_step/data_and_flags.f_obs.d_min(),
      use_kick_map      = params.map_1.use_kick_map,
      high_resolution   = hr,
      low_resolution    = lr)
    fft_map_1 = model_to_map_obj_1.fft_map
    map_1 = fft_map_1.real_map_unpadded()
    map_1_focus = map_1.focus()
    map_1_all = map_1.all()
    # prepare map CC calculation object
    map_cc_obj = get_map_cc_obj(map_1 = map_1, params = params,
      pdb_to_xrs_1 = pdb_to_xrs_1, pdb_to_xrs_2 = pdb_to_xrs_2,
      fft_map_1 = fft_map_1, atom_detail = atom_detail,
      residue_detail = residue_detail, atom_radius = atom_radius,
      hydrogen_atom_radius = params.hydrogen_atom_radius)
    #
    ### second
    if(xray_structure_2 is not None): xrs = xray_structure_2
    else: xrs = xray_structure_1
    if(data_and_flags_2 is not None): data_and_flags = data_and_flags_2
    else: data_and_flags = data_and_flags_1
    hr, lr = None, None
    if(params.map_2.high_resolution is not None):
      hr = params.map_2.high_resolution
    elif(params.map_1.high_resolution is not None):
      hr = params.map_1.high_resolution
    if(params.map_2.low_resolution is not None):
      lr = params.map_2.low_resolution
    elif(params.map_1.low_resolution is not None):
      lr = params.map_1.low_resolution
    model_to_map_obj_2 = model_to_map(
      xray_structure    = xrs,
      f_obs             = data_and_flags.f_obs,
      map_type          = params.map_2.map_type,
      resolution_factor = grid_step/data_and_flags.f_obs.d_min(),
      other_fft_map     = fft_map_1,
      use_kick_map      = params.map_2.use_kick_map,
      high_resolution   = hr,
      low_resolution    = lr)
    map_2 = model_to_map_obj_2.fft_map.real_map_unpadded()

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
    #
    atom_detail, residue_detail, atom_radius = set_details_level_and_radius(
      details_level = params.details_level,
      d_min         = hr,
      atom_radius   = params.atom_radius)
    grid_step = compute_grid_step_from_atom_radius_and_number_of_grid_points(
      r = atom_radius, n = params.number_of_grid_points, d_min = hr)
    #
    fft_map_1 = compute_map_from_model(
      high_resolution        = hr,
      low_resolution         = lr,
      xray_structure         = xrs,
      grid_resolution_factor = grid_step/hr,
      crystal_gridding = None)
    fft_map_1.apply_sigma_scaling()
    map_1 = fft_map_1.real_map_unpadded()
    map_1_focus = map_1.focus()
    map_1_all = map_1.all()
    # prepare map CC calculation object
    map_cc_obj = get_map_cc_obj(map_1 = map_1, params = params,
      pdb_to_xrs_1 = pdb_to_xrs_1, pdb_to_xrs_2 = pdb_to_xrs_2,
      fft_map_1 = fft_map_1, atom_detail = atom_detail,
      residue_detail = residue_detail, atom_radius = atom_radius,
      hydrogen_atom_radius = params.hydrogen_atom_radius)
    #
    if(xray_structure_2 is not None): xrs = xray_structure_2
    else: xrs = xray_structure_1
    if(params.map_2.high_resolution is not None):
      hr = params.map_2.high_resolution
    elif(params.map_1.high_resolution is not None):
      hr = params.map_1.high_resolution
    else:
      hr = d_min_default
    if(params.map_2.low_resolution is not None):
      lr = params.map_2.low_resolution
    elif(params.map_1.low_resolution is not None):
      lr = params.map_1.low_resolution
    else:
      lr = d_max_default
    fft_map_2 = compute_map_from_model(
      high_resolution        = hr,
      low_resolution         = lr,
      xray_structure         = xrs,
      grid_resolution_factor = grid_step/hr,
      crystal_gridding       = fft_map_1)
    fft_map_2.apply_sigma_scaling()
    map_2 = fft_map_2.real_map_unpadded()
  #
  assert map_1_focus == map_2.focus()
  assert map_1_all == map_2.all()
  # get map cc
  result = map_cc_obj.map_cc(
    map_2                                     = map_2,
    map_2_name                                = params.map_2.map_type,
    set_cc_to_zero_if_n_grid_points_less_than = params.set_cc_to_zero_if_n_grid_points_less_than,
    poor_cc_threshold                         = params.poor_cc_threshold,
    poor_map_value_threshold                  = params.poor_map_value_threshold)
  show_result(result = result, show_hydrogens = False)
  map_cc_obj.overall_correlation_min_max_standard_deviation()

class map_cc_funct(object):

  def __init__(self, map_1,
                     xray_structure,
                     fft_map,
                     atom_radius,
                     hydrogen_atom_radius,
                     atom_detail,
                     residue_detail,
                     map_1_name = None,
                     selection = None,
                     pdb_hierarchy = None):
    self.map_1_name = map_1_name
    self.xray_structure = xray_structure
    self.selection = selection
    self.atom_detail = atom_detail
    self.residue_detail = residue_detail
    self.pdb_hierarchy = pdb_hierarchy
    self.result = []
    self.map_1_size = map_1.size()
    self.map_1_stat = maptbx.statistics(map_1)
    assert [self.atom_detail, self.residue_detail].count(True) == 1
    self.atoms_with_labels = None
    if(pdb_hierarchy is not None and self.atom_detail):
      self.atoms_with_labels = list(pdb_hierarchy.atoms_with_labels())
    scatterers = self.xray_structure.scatterers()
    if(self.selection is None):
      self.selection = flex.bool(scatterers.size(), True)
    real_map_unpadded = fft_map.real_map_unpadded()
    sites_cart = self.xray_structure.sites_cart()
    if(self.atom_detail):
      self.gifes = [None,]*scatterers.size()
      self._result = [None,]*scatterers.size()
      #
      atom_radii = flex.double(scatterers.size(), atom_radius)
      for i_seq, sc in enumerate(scatterers):
        if(self.selection[i_seq]):
          if(sc.element_symbol().strip().lower() in ["h","d"]):
            atom_radii[i_seq] = hydrogen_atom_radius
      #
      for i_seq, site_cart in enumerate(sites_cart):
        if(self.selection[i_seq]):
          sel = maptbx.grid_indices_around_sites(
            unit_cell  = self.xray_structure.unit_cell(),
            fft_n_real = real_map_unpadded.focus(),
            fft_m_real = real_map_unpadded.all(),
            sites_cart = flex.vec3_double([site_cart]),
            site_radii = flex.double([atom_radii[i_seq]]))
          self.gifes[i_seq] = sel
          m1 = map_1.select(sel)
          ed1 = map_1.eight_point_interpolation(scatterers[i_seq].site)
          a = None
          if(self.atoms_with_labels is not None):
            a = self.atoms_with_labels[i_seq]
          self._result[i_seq] = group_args(atom = a, m1 = m1, ed1 = ed1,
            xyz=site_cart)
    if(self.residue_detail):
      assert self.pdb_hierarchy is not None
      residues = self.extract_residues()
      self.gifes = [None,]*len(residues)
      self._result = [None,]*len(residues)
      for i_seq, residue in enumerate(residues):
        residue_sites_cart = sites_cart.select(residue.selection)
        if 0: print i_seq, list(residue.selection) # DEBUG
        sel = maptbx.grid_indices_around_sites(
          unit_cell  = self.xray_structure.unit_cell(),
          fft_n_real = real_map_unpadded.focus(),
          fft_m_real = real_map_unpadded.all(),
          sites_cart = residue_sites_cart,
          site_radii = flex.double(residue.selection.size(), atom_radius))
        self.gifes[i_seq] = sel
        m1 = map_1.select(sel)
        ed1 = flex.double()
        for i_seq_r in residue.selection:
          ed1.append(map_1.eight_point_interpolation(scatterers[i_seq_r].site))
        self._result[i_seq] = \
          group_args(residue = residue, m1 = m1, ed1 = flex.mean(ed1),
            xyz=residue_sites_cart.mean(), n_atoms=residue_sites_cart.size())
    del map_1

  def map_cc(self, map_2,
                   set_cc_to_zero_if_n_grid_points_less_than,
                   poor_cc_threshold,
                   poor_map_value_threshold,
                   map_2_name = None):
    assert self.map_1_size == map_2.size()
    self.map_2_stat = maptbx.statistics(map_2)
    scatterers = self.xray_structure.scatterers()
    unit_cell = self.xray_structure.unit_cell()
    if(self.atom_detail):
      for i_seq, scatterer in enumerate(scatterers):
        if(self.selection[i_seq]):
          sel = list(self.gifes[i_seq])
          m1 = self._result[i_seq].m1
          m2 = map_2.select(sel)
          assert m1.size() == m2.size()
          if(m1.size() < set_cc_to_zero_if_n_grid_points_less_than or
            scatterer.occupancy == 0.0): corr = 0.
          else: corr = flex.linear_correlation(x = m1, y = m2).coefficient()
          ed1 = self._result[i_seq].ed1
          ed2 = map_2.eight_point_interpolation(scatterer.site)
          xyz = self._result[i_seq].xyz
          poor_flag = False
          if(((ed2 < poor_map_value_threshold or ed1 < poor_map_value_threshold)
             or corr < poor_cc_threshold)):
            poor_flag = True
          a = None
          if(self._result[i_seq].atom is not None):
            assert self._result[i_seq].atom is self.atoms_with_labels[i_seq]
            a = self._result[i_seq].atom
          if(scatterer.u_iso == -1):
            b_iso = adptbx.u_as_b(adptbx.u_star_as_u_iso(unit_cell, scatterer.u_star))
          else:
            b_iso = adptbx.u_as_b(scatterer.u_iso)
          self.result.append(group_args(
            atom             = a,
            b_iso            = b_iso,
            occupancy        = scatterer.occupancy,
            xyz              = xyz,
            cc               = corr,
            map_1_val        = ed1,
            map_2_val        = ed2,
            map_1_name       = self.map_1_name,
            map_2_name       = map_2_name,
            residual_map_val = None,
            scatterer        = scatterers[i_seq],
            data_points      = m2.size(),
            i_seq            = i_seq,
            poor_flag        = poor_flag))
    elif(self.residue_detail):
      scatterers = self.xray_structure.scatterers()
      occupancies = scatterers.extract_occupancies()
      b_isos = self.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
      for i, result in enumerate(self._result):
        sel = list(self.gifes[i])
        m1 = self._result[i].m1
        m2 = map_2.select(sel)
        assert m1.size() == m2.size()
        if(m1.size() < set_cc_to_zero_if_n_grid_points_less_than): corr = 0.
        else: corr = flex.linear_correlation(x = m1, y = m2).coefficient()
        ed1 = self._result[i].ed1
        ed2_ = flex.double()
        xyz = self._result[i].xyz
        n_atoms = self._result[i].n_atoms
        for i_seq in self._result[i].residue.selection:
          ed2_.append(map_2.eight_point_interpolation(scatterers[i_seq].site))
        ed2 = flex.mean(ed2_)
        poor_flag = False
        if(((ed2 < poor_map_value_threshold or ed1 < poor_map_value_threshold)
           or corr < poor_cc_threshold)):
          poor_flag = True
        self.result.append(group_args(
          residue     = result.residue,
          occupancy   = flex.mean(occupancies.select(result.residue.selection)),
          xyz         = xyz,
          n_atoms     = n_atoms,
          b_iso       = flex.mean(b_isos.select(result.residue.selection)),
          cc          = corr,
          map_1_val   = ed1,
          map_2_val   = ed2,
          data_points = m2.size(),
          poor_flag   = poor_flag))
    else: raise RuntimeError
    del map_2
    del self._result
    return self.result

  def overall_correlation_min_max_standard_deviation(self, log = None):
    if(log is None): log = sys.stdout
    fmt = "Map %d: min, max, mean, standard deviation: %6.3f %6.3f %6.3f %6.3f"
    print >> log, fmt%(1, self.map_1_stat.min(), self.map_1_stat.max(),
      self.map_1_stat.mean(), self.map_1_stat.sigma())
    print >> log, fmt%(2, self.map_2_stat.min(), self.map_2_stat.max(),
      self.map_2_stat.mean(), self.map_2_stat.sigma())

  def extract_residues(self, combine = True):
    result = []
    for i_model, model in enumerate(self.pdb_hierarchy.models()):
      rm = []
      for chain in model.chains():
        for rg in chain.residue_groups():
          rg_i_seqs = []
          r_name = None
          for ag in rg.atom_groups():
            if(r_name is None): r_name = ag.resname
            for atom in ag.atoms():
              if(self.selection[atom.i_seq]):
                rg_i_seqs.append(atom.i_seq)
          if(len(rg_i_seqs) != 0):
            rm.append(group_args(
              selection = flex.size_t(rg_i_seqs),
              name      = r_name,
              model_id  = i_model,
              resid     = rg.resid(),
              chain_id  = chain.id))
      result.append(rm)
    #
    if(combine):
      r0 = result[0]
      for r in result[1:]:
        for i, ri in enumerate(r):
          r0[i].selection.extend(ri.selection)
          assert r0[i].name == ri.name
    else:
      r0 = result[0]
      for r in result[1:]:
        r0.extend(r)
    return r0

def set_details_level_and_radius(details_level, d_min, atom_radius):
  assert details_level in ["atom","residue","automatic"]
  if(details_level == "automatic"):
    assert atom_radius is None
    if(d_min<2.5):
      atom_detail    = True
      residue_detail = False
      atom_radius    = 1.0
    else:
      atom_detail    = False
      residue_detail = True
      atom_radius    = 1.5
  elif(details_level == "atom"):
    atom_detail    = True
    residue_detail = False
    if(atom_radius is None): atom_radius = 1.0
  elif(details_level == "residue"):
    atom_detail    = False
    residue_detail = True
    if(atom_radius is None): atom_radius = 1.5
  return atom_detail, residue_detail, atom_radius

def simple(fmodel,
           pdb_hierarchy         = None,
           map_1_name            = "Fc",
           map_2_name            = "2mFo-DFc",
           details_level         = "automatic",
           atom_radius           = None,
           hydrogen_atom_radius  = 1.0,
           number_of_grid_points = 100,
           show                  = True,
           log                   = None,
           show_hydrogens        = False,
           selection             = None,
           diff_map              = "mFo-DFc",
           set_cc_to_zero_if_n_grid_points_less_than = 50,
           poor_cc_threshold                         = 0.7,
           poor_map_value_threshold                  = 1.0):
    if(fmodel.twin):
      raise Sorry("Not available for twinned data.")
    atom_detail, residue_detail, atom_radius = set_details_level_and_radius(
      details_level = details_level,
      d_min         = fmodel.f_obs.d_min(),
      atom_radius   = atom_radius)
    map_name_obj = mmtbx.map_names(map_name_string = map_1_name)
    grid_step = compute_grid_step_from_atom_radius_and_number_of_grid_points(
      r = atom_radius, n = number_of_grid_points, d_min = fmodel.f_obs.d_min())
    resolution_factor = grid_step/fmodel.f_obs.d_min()
    if([map_name_obj.k, map_name_obj.n] == [0,-1] and not map_name_obj.ml_map):
      complete_set = fmodel.f_obs.complete_set(d_min = fmodel.f_obs.d_min(),
        d_max=None)
      f_calc = complete_set.structure_factors_from_scatterers(
        xray_structure = fmodel.xray_structure).f_calc()
      fft_map_1 = f_calc.fft_map(
        resolution_factor = min(0.5,resolution_factor),
        symmetry_flags    = maptbx.use_space_group_symmetry)
    else:
      fft_map_1 = fmodel.electron_density_map().fft_map(
        resolution_factor = min(0.5,resolution_factor),
        map_type          = map_1_name,
        symmetry_flags    = maptbx.use_space_group_symmetry)
    fft_map_1.apply_sigma_scaling()
    map_1 = fft_map_1.real_map_unpadded()
    assert fmodel.xray_structure is not None
    map_cc_obj = map_cc_funct(
      map_1          = map_1,
      map_1_name     = map_1_name,
      xray_structure = fmodel.xray_structure,
      selection      = selection,
      fft_map        = fft_map_1,
      pdb_hierarchy  = pdb_hierarchy,
      atom_detail    = atom_detail,
      atom_radius    = atom_radius,
      hydrogen_atom_radius = hydrogen_atom_radius,
      residue_detail = residue_detail)
    del map_1
    fft_map_2 = fmodel.electron_density_map().fft_map(
      other_fft_map  = fft_map_1,
      map_type       = map_2_name,
      symmetry_flags = maptbx.use_space_group_symmetry)
    fft_map_2.apply_sigma_scaling()
    map_2 = fft_map_2.real_map_unpadded()
    result = map_cc_obj.map_cc(
      map_2                                     = map_2,
      map_2_name                                = map_2_name,
      set_cc_to_zero_if_n_grid_points_less_than = set_cc_to_zero_if_n_grid_points_less_than,
      poor_cc_threshold                         = poor_cc_threshold,
      poor_map_value_threshold                  = poor_map_value_threshold)
    del map_2
    if(atom_detail and diff_map is not None):
      fft_map_3 = fmodel.electron_density_map().fft_map(
        other_fft_map  = fft_map_1,
        map_type       = diff_map,
        symmetry_flags = maptbx.use_space_group_symmetry)
      fft_map_3.apply_sigma_scaling()
      map_3 = fft_map_3.real_map_unpadded()
      for i_seq, r in enumerate(result):
        ed3 = map_3.eight_point_interpolation(r.scatterer.site)
        r.residual_map_val = ed3
      del map_3
    if(show):
      show_result(result = result, show_hydrogens = show_hydrogens, log = log)
      map_cc_obj.overall_correlation_min_max_standard_deviation(log = log)
    return result

def show_result(result, show_hydrogens = False, log = None):
  if(log is None): log = sys.stdout
  keys = result[0].__dict__.keys()
  if("atom" in keys):
    assert not "residue" in keys
    if(result[0].residual_map_val is not None):
      print >> log, "i_seq :   PDB_string      element   occ      b      CC   %s   %s  mFo-DFc  No.Points  FLAG"%(result[0].map_1_name,result[0].map_2_name)
      fmt = "%5d : %s %7s %5.2f %6.2f %7.4f %6.2f %6.2f   %6.2f     %4d  %s"
      for i_seq, r in enumerate(result):
        w_msg = ""
        if(r.poor_flag): w_msg = " <<<"
        if(abs(r.residual_map_val) > 2.5): w_msg = " <<<"
        print_line = True
        if(not show_hydrogens and
           r.scatterer.element_symbol().strip().upper() in ["H","D"]):
          print_line = False
        if(print_line):
          print >> log, fmt % (
            i_seq,
            r.atom.id_str()[4:],
            r.atom.element,
            r.occupancy,
            r.b_iso,
            r.cc,
            r.map_1_val,
            r.map_2_val,
            r.residual_map_val,
            r.data_points,
            w_msg)
    else:
      print >> log, "i_seq :   PDB_string      element   occ      b      CC   %s   %s  No.Points FLAG"%(result[0].map_1_name,result[0].map_2_name)
      fmt = "%5d : %s %7s %5.2f %6.2f %7.4f %6.2f %6.2f   %d %s"
      for i_seq, r in enumerate(result):
        w_msg = ""
        if(r.poor_flag): w_msg = " <<<"
        print_line = True
        if(not show_hydrogens and
           r.scatterer.element_symbol().strip().upper() in ["H","D"]):
          print_line = False
        if(print_line):
          print >> log, fmt % (
            i_seq,
            r.atom.id_str()[4:],
            r.atom.element,
            r.occupancy,
            r.b_iso,
            r.cc,
            r.map_1_val,
            r.map_2_val,
            r.data_points,
            w_msg)
  elif("residue" in keys):
    assert not "atom" in keys
    print >> log, "i_seq : model chain resseq resname   occ      b      CC   map1   map2  No.Points"
    fmt = "%5d : %5s %5s %6s %7s %5.2f %6.2f %7.4f %6.2f %6.2f   %d"
    for i_seq, r in enumerate(result):
      print >> log, fmt % (
        i_seq,
        r.residue.model_id,
        r.residue.chain_id,
        r.residue.resid,
        r.residue.name,
        r.occupancy,
        r.b_iso,
        r.cc,
        r.map_1_val,
        r.map_2_val,
        r.data_points)
  else: raise RuntimeError
