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
from libtbx.utils import Sorry
import os, math, time
from cctbx import miller
from mmtbx import map_tools
from libtbx import adopt_init_args
from libtbx import Auto, group_args
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss

master_params_str = """\
grid_step = 0.5
  .type = float
  .help = Defines finess of the grid at which the map is computed.
  .expert_level = 2
scattering_table = *n_gaussian wk1995 it1992 neutron
  .type = choice
  .help = Scattering table for structure factors calculations
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
"""
master_params = iotbx.phil.parse(master_params_str, process_includes=False)

class model_to_map(object):
  def __init__(self, xray_structures, f_obs, map_type,
                     resolution_factor, high_resolution, low_resolution,
                     use_kick_map,
                     other_fft_map = None):
    # XXX
    map_name_obj = mmtbx.map_names(map_name_string = map_type)
    if(map_name_obj.k == 0 and map_name_obj.n == -1 and not map_name_obj.ml_map):
      complete_set = f_obs.complete_set(d_min = f_obs.d_min(), d_max=None)
      assert len(xray_structures) == 1
      assert not use_kick_map
      f_calc = complete_set.structure_factors_from_scatterers(
        xray_structure = xray_structures[0],
        grid_resolution_factor = resolution_factor).f_calc()
      if(other_fft_map is not None):
        self.fft_map = miller.fft_map(
          crystal_gridding     = other_fft_map,
          fourier_coefficients = f_calc)
      else:
        self.fft_map = f_calc.fft_map(
          resolution_factor = resolution_factor,
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
      self.fmodel = utils.fmodel_simple(
        xray_structures = xray_structures,
        f_obs           = f_obs,
        r_free_flags    = r_free_flags,
        bss_params      = bss_params,
        twin_law        = None) # XXX support twin_law in future
      if(not use_kick_map):
        map_obj = self.fmodel.electron_density_map()
        map_coeff = map_obj.map_coefficients(map_type = map_type)
        self.fft_map = map_obj.fft_map(
          resolution_factor = resolution_factor,
          map_coefficients  = map_coeff,
          other_fft_map     = other_fft_map)
        self.fft_map.apply_sigma_scaling()
        self.map_data =  self.fft_map.real_map_unpadded()
      else:
        km = map_tools.kick_map(
          fmodel                        = self.fmodel,
          map_type                      = map_type,
          resolution_factor             = resolution_factor,
          other_fft_map                 = other_fft_map,
          real_map                      = True,
          real_map_unpadded             = False,
          symmetry_flags                = maptbx.use_space_group_symmetry,
          average_maps                  = True) # XXX use map coefficients averaging
        self.fft_map = km.fft_map
        self.map_data = km.map_data

class pdb_to_xrs(object):
  def __init__(self, pdb_files, crystal_symmetry, scattering_table):
    processed_pdb_files_srv = utils.process_pdb_file_srv(
      crystal_symmetry = crystal_symmetry,
      log = StringIO())
    self.xray_structures = None
    self.processed_pdb_file = None
    if(pdb_files != []):
      self.processed_pdb_file, pdb_inp = \
        processed_pdb_files_srv.process_pdb_files(pdb_file_names = pdb_files)
      xsfppf = mmtbx.utils.xray_structures_from_processed_pdb_file(
        processed_pdb_file = self.processed_pdb_file,
        scattering_table   = scattering_table,
        d_min              = None)
      self.xray_structures = xsfppf.xray_structures

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
    parameters.force_anomalous_flag_to_be_equal_to = False
    if(params.data_labels is not None):
      parameters.labels = [params.data_labels]
    if(params.high_resolution is not None):
      parameters.high_resolution = params.high_resolution
    if(params.low_resolution is not None):
      parameters.low_resolution = params.low_resolution
    data_and_flags = utils.determine_data_and_flags(
      reflection_file_server = reflection_file_server,
      parameters             = parameters,
      data_description       = "X-ray data",
      extract_r_free_flags   = False,
      log                    = StringIO())
  return data_and_flags

def compute_map_from_model(high_resolution, low_resolution, xray_structures,
                           grid_resolution_factor, crystal_gridding = None):
  if(len(xray_structures) > 1):
    raise Sorry("Multiple models cannot be used with this option.")
  xray_structure = xray_structures[0]
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

Description:
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
    master_params.show()
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

def run(params, d_min_default=1.5, d_max_default=999.9) :
  # check for crystal_symmetry
  crystal_symmetry = extract_crystal_symmetry(params = params)
  # read in the PDB files
  if(params.map_1.pdb_file_name+params.map_2.pdb_file_name == []):
    raise Sorry("No coordinate file given.")
  pdb_to_xrs_1 = pdb_to_xrs(pdb_files        = params.map_1.pdb_file_name,
                            crystal_symmetry = crystal_symmetry,
                            scattering_table = params.scattering_table)
  xray_structures_1 = pdb_to_xrs_1.xray_structures
  pdb_to_xrs_2 = pdb_to_xrs(pdb_files        = params.map_2.pdb_file_name,
                            crystal_symmetry = crystal_symmetry,
                            scattering_table = params.scattering_table)
  xray_structures_2 = pdb_to_xrs_2.xray_structures
  # assert correct combination of options
  if([params.map_1.use, params.map_2.use].count(True) == 0):
    raise Sorry("No model selected to compute local density correlation.")
  if([params.map_1.use, params.map_2.use].count(True) == 2):
    raise Sorry(
      "Select one model to compute local density correlation (two selected).")
  if(params.map_1.use):
    if(xray_structures_1 is None):
      raise Sorry("PDB file for map_1 is not provided.")
    if(len(xray_structures_1) > 1):
      raise Sorry("use=true option cannot be used for a PDB file with multiple models.")
      if(params.map_1.use_kick_map):
        raise Sorry("kick map cannot be used for a PDB file with multiple models.")
  if(params.map_2.use):
    if(xray_structures_2 is None):
      raise Sorry("PDB file for map_2 is not provided.")
    if(len(xray_structures_2) > 1):
      raise Sorry("use=true option cannot be used for a PDB file with multiple models.")
      if(params.map_2.use_kick_map):
        raise Sorry("kick map cannot be used for a PDB file with multiple models.")
  # read in F-obs and free-r flags for map_1
  data_and_flags_1 = extract_data_and_flags(
    params           = params.map_1,
    crystal_symmetry = crystal_symmetry)
  data_and_flags_2 = extract_data_and_flags(
    params           = params.map_2,
    crystal_symmetry = crystal_symmetry)
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
    if(xray_structures_1 is not None): xrs = xray_structures_1
    else: xrs = xray_structures_2
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
    model_to_map_obj_1 = model_to_map(
      xray_structures   = xrs,
      f_obs             = data_and_flags.f_obs,
      map_type          = params.map_1.map_type,
      resolution_factor = params.grid_step/data_and_flags.f_obs.d_min(),
      use_kick_map      = params.map_1.use_kick_map,
      high_resolution   = hr,
      low_resolution    = lr)
    fft_map_1 = model_to_map_obj_1.fft_map
    map_1 = fft_map_1.real_map_unpadded()
    ### second
    if(xray_structures_2 is not None): xrs = xray_structures_2
    else: xrs = xray_structures_1
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
      xray_structures   = xrs,
      f_obs             = data_and_flags.f_obs,
      map_type          = params.map_2.map_type,
      resolution_factor = params.grid_step/data_and_flags.f_obs.d_min(),
      other_fft_map     = fft_map_1,
      use_kick_map      = params.map_2.use_kick_map,
      high_resolution   = hr,
      low_resolution    = lr)
    map_2 = model_to_map_obj_2.fft_map.real_map_unpadded()
  # compute data if no reflection files provided for both maps
  if([data_and_flags_1, data_and_flags_2].count(None) == 2):
    if(xray_structures_1 is not None): xrs = xray_structures_1
    else: xrs = xray_structures_2
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
      xray_structures        = xrs,
      grid_resolution_factor = params.grid_step/hr,
      crystal_gridding = None)
    fft_map_1.apply_sigma_scaling()
    map_1 = fft_map_1.real_map_unpadded()
    if(xray_structures_2 is not None): xrs = xray_structures_2
    else: xrs = xray_structures_1
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
      xray_structures        = xrs,
      grid_resolution_factor = params.grid_step/hr,
      crystal_gridding       = fft_map_1)
    fft_map_2.apply_sigma_scaling()
    map_2 = fft_map_2.real_map_unpadded()
  #
  assert map_1.focus() == map_2.focus()
  assert map_1.all() == map_2.all()
  # get sampled points
  if(params.map_1.use): xrs = xray_structures_1
  elif(params.map_2.use): xrs = xray_structures_2
  else: RuntimeError
  if(len(xrs) > 1):
    raise Sorry("Multiple models cannot be used with this option.")
  # get map cc
  pdb_to_xrs_ = None
  if(params.map_1.use): pdb_to_xrs_ = pdb_to_xrs_1
  elif(params.map_2.use): pdb_to_xrs_ = pdb_to_xrs_2
  else: RuntimeError
  compute_map_cc_result = compute_map_cc(map_1 = map_1, map_2 = map_2,
    xray_structure = xrs[0], fft_map = fft_map_1)
  per_atom_result = compute_map_cc_result.atoms(
    pdb_hierarchy = pdb_to_xrs_.processed_pdb_file.all_chain_proxies.\
      pdb_hierarchy, show = True)
  compute_map_cc_result.overall_correlation_min_max_standard_deviation(
    show = True)

class compute_map_cc(object):

  def __init__(self, map_1, map_2, xray_structure, fft_map):
    self.map_1 = map_1
    self.map_2 = map_2
    assert self.map_1.size() == self.map_1.size()
    self.xray_structure = xray_structure
    self.fft_map = fft_map
    self.sampled_density = self.sampled_density_map()

  def overall_correlation_min_max_standard_deviation(self, show = False):
    corr = flex.linear_correlation(
      x = self.map_1.as_1d(),
      y = self.map_2.as_1d()).coefficient()
    st1 = maptbx.statistics(self.map_1)
    st2 = maptbx.statistics(self.map_2)
    if(show):
      print "Overall map correlation: %5.2f"%corr
      print "Map 1: min, max, mean, standard deviation: %6.3f %6.3f %6.3f %6.3f" % (
        st1.min(), st1.max(), st1.mean(), st1.sigma())
      print "Map 2: min, max, mean, standard deviation: %6.3f %6.3f %6.3f %6.3f" % (
        st2.min(), st2.max(), st2.mean(), st2.sigma())
    return group_args(correlation = corr, map_stat_1 = st1, map_stat_2 = st2)

  def sampled_density_map(self):
    # Very specific code to get the sampling points only.
    # Do not use the actual map.
    assert not self.fft_map.anomalous_flag()
    xrs = self.xray_structure.deep_copy_scatterers()
    xrs.convert_to_isotropic()
    u_iso_all = adptbx.b_as_u(30.0)
    xrs.set_u_iso(value = u_iso_all)
    real_map_unpadded = self.fft_map.real_map_unpadded()
    sampled_density = ext.sampled_model_density(
      unit_cell                             = xrs.unit_cell(),
      scatterers                            = xrs.scatterers(),
      scattering_type_registry              = xrs.scattering_type_registry(),
      fft_n_real                            = real_map_unpadded.focus(),
      fft_m_real                            = real_map_unpadded.all(),
      u_base                                = adptbx.b_as_u(0),
      wing_cutoff                           = 1.e-1,
      exp_table_one_over_step_size          = -100,
      force_complex                         = False,
      sampled_density_must_be_positive      = False,
      tolerance_positive_definite           = 1.e-5,
      use_u_base_as_u_extra                 = True,
      store_grid_indices_for_each_scatterer = -1)
    return sampled_density

  def show(self, histogram):
    h_1 = histogram
    lc_1 = histogram.data_min()
    s_1 = enumerate(histogram.slots())
    for (i_1,n_1) in s_1:
      hc_1 = h_1.data_min() + h_1.slot_width() * (i_1+1)
      print "%8.3f - %8.3f: %5d" % (lc_1,hc_1,n_1)
      lc_1 = hc_1

  def atoms(self, pdb_hierarchy,
                  show = False,
                  poor_cc_threshold = 0.7,
                  poor_map_value_threshold = 1.0,
                  ignore_points_with_map_values_less_than = 0.1,
                  set_cc_to_zero_if_n_grid_points_less_than = 50,
                  show_hydrogens = True):
    result = []
    atoms = pdb_hierarchy.atoms_with_labels()
    gifes = self.sampled_density.grid_indices_for_each_scatterer()
    scatterers = self.xray_structure.scatterers()
    for i_seq, a in enumerate(atoms):
       sel = list(gifes[i_seq])
       m1 = self.map_1.select(sel)
       m2 = self.map_2.select(sel)
       sel_flat  = m1 > ignore_points_with_map_values_less_than
       sel_flat &= m2 > ignore_points_with_map_values_less_than
       m1 = m1.select(sel_flat)
       m2 = m2.select(sel_flat)
       if(m1.size() < set_cc_to_zero_if_n_grid_points_less_than): corr = 0.
       else:
         corr = flex.linear_correlation(
           x = m1,
           y = m2).coefficient()
       ed1 = self.map_1.eight_point_interpolation(scatterers[i_seq].site)
       ed2 = self.map_2.eight_point_interpolation(scatterers[i_seq].site)
       poor_flag = False
       if(corr < poor_cc_threshold or
          ((ed2 < poor_map_value_threshold or ed1 < poor_map_value_threshold)
          and corr < poor_cc_threshold) or
          ed2 < ignore_points_with_map_values_less_than or
          ed1 < ignore_points_with_map_values_less_than):
         poor_flag = True
       result.append(group_args(
         atom      = a,
         cc        = corr,
         map_1_val = ed1,
         map_2_val = ed2,
         poor_flag = poor_flag))
    if(show):
      print "i_seq : chain resseq resname altloc name element   occ      b      CC   map1   map2  FLAG"
      fmt = "%5d : %5s %6s %7s %6s %4s %7s %5.2f %6.2f %7.4f %6.2f %6.2f %s"
      cc_values = flex.double()
      map_values = flex.double()
      for i_seq, r in enumerate(result):
        assert scatterers[i_seq].element_symbol().strip().upper() == \
          r.atom.element.strip().upper()
        w_msg = ""
        if(r.poor_flag): w_msg = " <<< WEAK DENSITY"
        print_line = True
        if(not show_hydrogens and r.atom.element.strip().upper() in ["H","D"]):
          print_line = False
        if(print_line):
          cc_values.append(r.cc)
          map_values.append(r.map_2_val)
          print fmt % (i_seq, r.atom.chain_id, r.atom.resseq, r.atom.resname,
            r.atom.altloc, r.atom.name, r.atom.element, r.atom.occ, r.atom.b,
            r.cc, r.map_1_val, r.map_2_val, w_msg)
    return result
