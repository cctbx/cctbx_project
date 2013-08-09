from __future__ import division
import mmtbx.utils
import iotbx.phil
from scitbx.array_family import flex
from libtbx.utils import Sorry, date_and_time
from libtbx import adopt_init_args
from libtbx.str_utils import show_string
from libtbx.math_utils import ifloor, iceil
import libtbx.callbacks # import dependency
import os
import sys
import random
from mmtbx import map_tools
from cctbx import miller
from cctbx import maptbx
from libtbx import group_args

map_coeff_params_base_str = """\
  map_coefficients
    .multiple = True
    .short_caption = Map coefficients
    .style = auto_align
  {
    map_type = None
      .type = str
      .style = bold renderer:draw_map_type_widget
    format = *mtz phs
      .type = choice(multi=True)
    mtz_label_amplitudes = None
      .type = str
      .short_caption = MTZ label for amplitudes
      .style = bold
    mtz_label_phases = None
      .type = str
      .short_caption = MTZ label for phases
      .style = bold
    kicked = False
      .type = bool
      .short_caption = Kicked map
    fill_missing_f_obs = False
      .type = bool
      .short_caption = Fill missing F(obs) with F(calc)
    acentrics_scale = 2.0
      .type = float
      .help = Scale terms corresponding to acentric reflections (residual maps only: k==n)
      .expert_level = 2
    centrics_pre_scale = 1.0
      .type = float
      .help = Centric reflections, k!=n and k*n != 0: \
              max(k-centrics_pre_scale,0)*Fo-max(n-centrics_pre_scale,0)*Fc
      .expert_level = 2
    sharpening = False
      .type = bool
      .help = Apply B-factor sharpening
      .short_caption = Apply B-factor sharpening
      .style = bold
    sharpening_b_factor = None
      .type = float
      .help = Optional sharpening B-factor value
      .short_caption = Sharpening B-factor value (optional)
    exclude_free_r_reflections = False
      .type = bool
      .help = Exclude free-R selected reflections from output map coefficients
      .short_caption = Exclude R-free set from map coefficients
    isotropize = True
      .type = bool
    dev
      .expert_level=3
    {
      complete_set_up_to_d_min = False
        .type = bool
      aply_same_incompleteness_to_complete_set_at = randomly low high
        .type = choice(multi=False)
    }
    %s
  }
"""

ncs_average_param_str = """
ncs_average = False
  .type = bool
  .expert_level = 2
  .help = Perform NCS averaging on map using RESOLVE (without density \
      modification).  Will be ignored if NCS is not present.
  .short_caption = NCS average
"""

# for phenix.maps
map_coeff_params_str = map_coeff_params_base_str % ""
# for phenix.refine
map_coeff_params_ncs_str = map_coeff_params_base_str % ncs_average_param_str

map_params_base_str ="""\
  map
    .short_caption = XPLOR or CCP4 map
    .multiple = True
    .style = auto_align
  {
    map_type = None
      .type = str
      .expert_level=0
      .style = bold renderer:draw_map_type_widget
    format = xplor *ccp4
      .type = choice
      .short_caption = File format
      .caption = XPLOR CCP4
      .style = bold
    file_name = None
      .type = path
      .style = bold new_file
    kicked = False
      .type = bool
      .expert_level=0
    fill_missing_f_obs = False
      .type = bool
      .expert_level=0
    grid_resolution_factor = 1/4.
      .type = float
      .expert_level=0
    scale = *sigma volume
      .type = choice(multi=False)
      .expert_level=2
    region = *selection cell
      .type = choice
      .caption = Atom_selection Unit_cell
      .short_caption=Map region
    atom_selection = None
      .type = atom_selection
      .short_caption = Atom selection
    atom_selection_buffer = 3
      .type = float
    acentrics_scale = 2.0
      .type = float
      .help = Scale terms corresponding to acentric reflections (residual maps only: k==n)
      .expert_level=2
    centrics_pre_scale = 1.0
      .type = float
      .help = Centric reflections, k!=n and k*n != 0: \
              max(k-centrics_pre_scale,0)*Fo-max(n-centrics_pre_scale,0)*Fc
      .expert_level=2
    sharpening = False
      .type = bool
      .help = Apply B-factor sharpening
      .short_caption = Apply B-factor sharpening
      .style = bold
    sharpening_b_factor = None
      .type = float
      .help = Optional sharpening B-factor value
      .short_caption = Sharpening B-factor value (optional)
    exclude_free_r_reflections = False
      .type = bool
      .help = Exclude free-R selected reflections from map calculation
    isotropize = True
      .type = bool
    %s
  }
"""

map_params_str = map_params_base_str % ""
map_params_ncs_str = map_params_base_str % ncs_average_param_str

# XXX for phenix.maps
map_and_map_coeff_params_str = """\
%s
%s
"""%(map_coeff_params_str, map_params_str)

# XXX for phenix.refine
map_and_map_coeff_params_ncs_str = """\
%s
%s
"""%(map_coeff_params_ncs_str, map_params_ncs_str)


def map_and_map_coeff_master_params():
  return iotbx.phil.parse(map_and_map_coeff_params_str, process_includes=False)

maps_including_IO_params_str = """\
maps {
  input {
    pdb_file_name = None
      .type = path
      .optional = False
      .short_caption = Model file
      .style = bold file_type:pdb input_file
    reflection_data {
      %s
      r_free_flags {
        %s
      }
    }
  }
  output {
    directory = None
      .type = path
      .short_caption = Output directory
      .help = For GUI only.
      .style = bold output_dir noauto
    prefix = None
      .type = str
      .input_size = 100
      .short_caption = Output prefix
      .style = bold noauto
    title = None
      .type = str
      .short_caption = Job title
      .input_size = 400
      .style = noauto
    fmodel_data_file_format = mtz
      .optional=True
      .type=choice
      .help=Write Fobs, Fmodel, various scales and more to MTZ file
    include_r_free_flags = False
      .type = bool
      .short_caption = Include R-free flags in output MTZ file
  }
  scattering_table = wk1995  it1992  *n_gaussian  neutron
    .type = choice
    .help = Choices of scattering table for structure factors calculations
  wavelength = None
    .type = float(value_min=0.2, value_max=10.)
    .input_size = 80
    .help = Optional X-ray wavelength (in Angstroms), which will be used to \
      set the appropriate anomalous scattering factors for the model.  This \
      will only affect the LLG map from Phaser.
  bulk_solvent_correction = True
    .type = bool
  anisotropic_scaling = True
    .type = bool
  skip_twin_detection = False
    .type = bool
    .short_caption = Skip automatic twinning detection
    .help = Skip automatic twinning detection
  omit {
    method = *simple
      .type = choice(multi=False)
    selection = None
      .type = str
      .short_caption = Omit selection
      .input_size = 400
  }
  %s
  %s
}
"""%(mmtbx.utils.data_and_flags_str_part1,
     mmtbx.utils.data_and_flags_str_part2,
     map_coeff_params_str,
     map_params_str)

# XXX for documentation
master_params = maps_including_IO_params_str

def maps_including_IO_master_params():
  return iotbx.phil.parse(maps_including_IO_params_str, process_includes=True)

def cast_map_coeff_params(map_type_obj):
  map_coeff_params_str = """\
    map_coefficients
    {
      format = *mtz phs
      mtz_label_amplitudes = %s
      mtz_label_phases = P%s
      map_type = %s
      kicked = %s
      fill_missing_f_obs = %s
    }
"""%(map_type_obj.format(), map_type_obj.format(), map_type_obj.format(),
     map_type_obj.kicked, map_type_obj.f_obs_filled)
  return iotbx.phil.parse(map_coeff_params_str, process_includes=False)

class map_coeffs_mtz_label_manager:

  def __init__(self, map_params):
    self._amplitudes = map_params.mtz_label_amplitudes
    self._phases = map_params.mtz_label_phases
    if(self._amplitudes is None): self._amplitudes = str(map_params.map_type)
    if(self._phases is None): self._phases = "PH"+str(map_params.map_type)

  def amplitudes(self):
    return self._amplitudes

  def phases(self, root_label, anomalous_sign=None):
    assert anomalous_sign is None or not anomalous_sign
    return self._phases

class write_xplor_map_file(object):

  def __init__(self, params, coeffs, atom_selection_manager=None,
               xray_structure=None):
    adopt_init_args(self, locals())
    fft_map = coeffs.fft_map(resolution_factor =
      self.params.grid_resolution_factor)
    if(self.params.scale == "volume"): fft_map.apply_volume_scaling()
    elif(self.params.scale == "sigma"): fft_map.apply_sigma_scaling()
    else: raise RuntimeError
    title_lines=["REMARK file: %s" %
      show_string(os.path.basename(self.params.file_name))]
    title_lines.append("REMARK directory: %s" %
      show_string(os.path.dirname(self.params.file_name)))
    title_lines.append("REMARK %s" % date_and_time())
    assert self.params.region in ["selection", "cell"]
    if(self.params.region == "selection" and xray_structure is not None) :
      map_iselection = None
      if atom_selection_manager is not None :
        map_iselection = self.atom_iselection()
      frac_min, frac_max = self.box_around_selection(
        iselection = map_iselection,
        buffer     = self.params.atom_selection_buffer)
      n_real = fft_map.n_real()
      gridding_first=[ifloor(f*n) for f,n in zip(frac_min,n_real)]
      gridding_last=[iceil(f*n) for f,n in zip(frac_max,n_real)]
      title_lines.append('REMARK map around selection')
      title_lines.append('REMARK   atom_selection=%s' %
        show_string(self.params.atom_selection))
      title_lines.append('REMARK   atom_selection_buffer=%.6g' %
        self.params.atom_selection_buffer)
      if(map_iselection is None):
        sel_size = self.xray_structure.scatterers().size()
      else:
        sel_size = map_iselection.size()
      title_lines.append('REMARK   number of atoms selected: %d' % sel_size)
    else:
      gridding_first = None
      gridding_last = None
      title_lines.append("REMARK map covering the unit cell")
    if params.format == "xplor" :
      fft_map.as_xplor_map(
        file_name      = self.params.file_name,
        title_lines    = title_lines,
        gridding_first = gridding_first,
        gridding_last  = gridding_last)
    else :
      fft_map.as_ccp4_map(
        file_name      = self.params.file_name,
        gridding_first = gridding_first,
        gridding_last  = gridding_last,
        labels=title_lines)

  def box_around_selection(self, iselection, buffer):
    sites_cart = self.xray_structure.sites_cart()
    if(iselection is not None):
      sites_cart = sites_cart.select(iselection)
    return self.xray_structure.unit_cell().box_frac_around_sites(
      sites_cart = sites_cart, buffer = buffer)

  def atom_iselection(self):
    if(self.params.region != "selection" or self.params.atom_selection is None):
      return None
    try:
      result = self.atom_selection_manager.selection(string =
        self.params.atom_selection).iselection()
    except KeyboardInterrupt: raise
    except Exception:
      raise Sorry('Invalid atom selection: %s' % self.params.atom_selection)
    if(result.size() == 0):
      raise Sorry('Empty atom selection: %s' % self.params.atom_selection)
    return result

def compute_f_calc(fmodel, params):
  from cctbx import miller
  coeffs_partial_set = fmodel.f_obs().structure_factors_from_scatterers(
    xray_structure = fmodel.xray_structure).f_calc()
  if(hasattr(params,"dev") and params.dev.complete_set_up_to_d_min):
    coeffs = fmodel.xray_structure.structure_factors(
      d_min = fmodel.f_obs().d_min()).f_calc()
    frac_inc = 1.*coeffs_partial_set.data().size()/coeffs.data().size()
    n_miss = coeffs.data().size() - coeffs_partial_set.data().size()
    if(params.dev.aply_same_incompleteness_to_complete_set_at == "randomly"):
      sel = flex.random_bool(coeffs.data().size(), frac_inc)
      coeffs = coeffs.select(sel)
    elif(params.dev.aply_same_incompleteness_to_complete_set_at == "low"):
      coeffs = coeffs.sort()
      coeffs = miller.set(
        crystal_symmetry = coeffs,
        indices = coeffs.indices()[n_miss+1:],
        anomalous_flag = coeffs.anomalous_flag()).array(
        data = coeffs.data()[n_miss+1:])
    elif(params.dev.aply_same_incompleteness_to_complete_set_at == "high"):
      coeffs = coeffs.sort(reverse=True)
      coeffs = miller.set(
        crystal_symmetry = coeffs,
        indices = coeffs.indices()[n_miss+1:],
        anomalous_flag = coeffs.anomalous_flag()).array(
        data = coeffs.data()[n_miss+1:])
  else:
    coeffs = coeffs_partial_set
  return coeffs

def map_coefficients_from_fmodel(
      params,
      fmodel = None,
      map_calculation_server = None,
      post_processing_callback=None,
      pdb_hierarchy=None):
  assert [fmodel, map_calculation_server].count(None) == 1
  from mmtbx import map_tools
  import mmtbx
  from cctbx import miller
  mnm = mmtbx.map_names(map_name_string = params.map_type)
  if(mnm.k==0 and abs(mnm.n)==1):
    return compute_f_calc(fmodel, params)
  if(fmodel is not None and
     fmodel.is_twin_fmodel_manager() and
     mnm.phaser_sad_llg):
    return None
  if(fmodel is not None):
    e_map_obj = fmodel.electron_density_map(update_f_part1=True)
    xrs = fmodel.xray_structure
  else:
    e_map_obj = map_calculation_server
    xrs = map_calculation_server.fmodel.xray_structure
  coeffs = None
  if(not params.kicked):
    coeffs = e_map_obj.map_coefficients(
      map_type           = params.map_type,
      acentrics_scale    = params.acentrics_scale,
      centrics_pre_scale = params.centrics_pre_scale,
      fill_missing       = params.fill_missing_f_obs,
      isotropize         = params.isotropize,
      exclude_free_r_reflections=params.exclude_free_r_reflections,
      ncs_average=getattr(params, "ncs_average", False),
      post_processing_callback=post_processing_callback,
      pdb_hierarchy=pdb_hierarchy)
    if (coeffs is None) : return None
    if(coeffs.anomalous_flag()) :
      coeffs = coeffs.average_bijvoet_mates()
    if(params.sharpening):
      from mmtbx import map_tools
      coeffs, b_sharp = map_tools.sharp_map(
        sites_frac = xrs.sites_frac(),
        map_coeffs = coeffs,
        b_sharp    = params.sharpening_b_factor)
  else:
    if(params.map_type.count("anom")==0):
      coeffs = kick(
        fmodel   = e_map_obj.fmodel,
        map_type = params.map_type).map_coefficients
  # XXX need to figure out why this happens
  if (coeffs is None) :
    raise RuntimeError(("Map coefficient generation failed (map_type=%s, "
      "kicked=%s, sharpening=%s, isotropize=%s, anomalous=%s.") %
        (params.map_type, params.kicked, params.sharpening, params.isotropize,
         fmodel.f_obs().anomalous_flag()))
  # XXX is this redundant?
  if(coeffs.anomalous_flag()) :
    coeffs = coeffs.average_bijvoet_mates()
  return coeffs

def compute_xplor_maps(
    fmodel,
    params,
    atom_selection_manager=None,
    file_name_prefix=None,
    file_name_base=None,
    post_processing_callback=None) :
  assert ((post_processing_callback is None) or
          (hasattr(post_processing_callback, "__call__")))
  output_files = []
  for mp in params:
    if(mp.map_type is not None):
      coeffs = map_coefficients_from_fmodel(fmodel = fmodel, params = mp,
        post_processing_callback=post_processing_callback)
      if (coeffs is None) :
        raise Sorry("Couldn't generate map type '%s'." % mp.map_type)
      if(mp.file_name is None):
        output_file_name = ""
        if(file_name_prefix is not None): output_file_name = file_name_prefix
        if(file_name_base is not None):
          if(len(output_file_name)>0):
            output_file_name = output_file_name + "_"+file_name_base
          else: output_file_name = output_file_name + file_name_base
        if mp.format == "xplor" :
          ext = ".xplor"
        else :
          ext = ".ccp4"
        output_file_name = output_file_name + "_" + mp.map_type + "_map" + ext
        mp.file_name = output_file_name
      write_xplor_map_file(params = mp, coeffs = coeffs,
        atom_selection_manager = atom_selection_manager,
        xray_structure = fmodel.xray_structure)
      output_files.append(mp.file_name)
  return output_files

class compute_map_coefficients(object):

  def __init__(self,
               fmodel,
               params,
               mtz_dataset = None,
               post_processing_callback=None,
               pdb_hierarchy=None,
               log=sys.stdout):
    assert ((post_processing_callback is None) or
            (hasattr(post_processing_callback, "__call__")))
    self.mtz_dataset = mtz_dataset
    coeffs = None
    # Avoid doing slow calculation several times!
    map_calculation_server = fmodel.electron_density_map(update_f_part1=True)
    self.map_coeffs = []
    for mcp in params:
      if(mcp.map_type is not None):
        if(fmodel.is_twin_fmodel_manager() and mcp.kicked) :
          raise Sorry("Kicked maps are not supported when twinning is present.  "+
            "You can disable the automatic twin law detection by setting the "+
            "parameter maps.skip_twin_detection to True (or check the "+
            "corresponding box in the Phenix GUI).")
        if(fmodel.is_twin_fmodel_manager()) and (mcp.isotropize) :
          mcp.isotropize = False
        coeffs = map_coefficients_from_fmodel(
          map_calculation_server   = map_calculation_server,
          params                   = mcp,
          post_processing_callback = post_processing_callback,
          pdb_hierarchy            = pdb_hierarchy)
        if("mtz" in mcp.format and coeffs is not None):
          lbl_mgr = map_coeffs_mtz_label_manager(map_params = mcp)
          if(self.mtz_dataset is None):
            self.mtz_dataset = coeffs.as_mtz_dataset(
              column_root_label = lbl_mgr.amplitudes(),
              label_decorator   = lbl_mgr)
          else:
            self.mtz_dataset.add_miller_array(
              miller_array      = coeffs,
              column_root_label = lbl_mgr.amplitudes(),
              label_decorator   = lbl_mgr)
          self.map_coeffs.append(coeffs)
        elif (coeffs is None) :
          if ((mcp.map_type == "anomalous") and
              (not fmodel.f_obs().anomalous_flag())) :
            # since anomalous map is included in the defaults, even if the
            # data are merged, no warning is issued here
            pass
          else :
            libtbx.warn(("Map coefficients not available for map type '%s'; "+
              "usually means you have requested an anomalous map but supplied "+
              "merged data, or indicates a twinning-related incompatibility.")%
              mcp.map_type)

  def write_mtz_file(self, file_name, mtz_history_buffer = None,
      r_free_flags=None):
    from cctbx.array_family import flex
    if(self.mtz_dataset is not None):
      if (r_free_flags is not None) :
        self.mtz_dataset.add_miller_array(r_free_flags,
          column_root_label="FreeR_flag")
      if(mtz_history_buffer is None):
        mtz_history_buffer = flex.std_string()
      mtz_history_buffer.append(date_and_time())
      mtz_history_buffer.append("> file name: %s" % os.path.basename(file_name))
      mtz_object = self.mtz_dataset.mtz_object()
      mtz_object.add_history(mtz_history_buffer)
      mtz_object.write(file_name = file_name)
      return True
    return False

class kick(object):

  def __init__(
      self,
      fmodel,
      crystal_gridding = None,
      map_type         = "2mFo-DFc",
      number_of_kicks  = 100,
      delta_fofc_outliers_threshold = 20,
      resolution_factor = 0.25,
      update_r_free_flags= True,
      use_complete_with  = True,
      use_intersection   = True,
      shake_f_calc = True,
      shake_f_mask = True,
      shake_f_model= True,
      shake_f_map  = True):
    self.map_type = map_type
    self.update_r_free_flags = update_r_free_flags
    self.delta_fofc_outliers_threshold = delta_fofc_outliers_threshold
    self.mean_phase_error = None
    self.use_complete_with = use_complete_with
    self.map_type = map_type
    self.number_of_kicks = number_of_kicks
    assert self.number_of_kicks > 0
    self.shake_f_calc  = shake_f_calc
    self.shake_f_mask  = shake_f_mask
    self.shake_f_model = shake_f_model
    self.shake_f_map   = shake_f_map
    if(crystal_gridding is None):
      crystal_gridding = fmodel.f_obs().crystal_gridding(
        d_min                   = fmodel.f_obs().d_min(),
        resolution_factor       = resolution_factor,
        grid_step               = None,
        symmetry_flags          = None,
        mandatory_factors       = None,
        max_prime               = 5,
        assert_shannon_sampling = True)
    fmodel = self.convert_to_non_anomalous(fmodel=fmodel)
    # starting map coefficients and sigma-scaled map
    self.complete_set = fmodel.electron_density_map(
      update_f_part1=False).map_coefficients(
        map_type     = self.map_type,
        isotropize   = True,
        fill_missing = True)
    fft_map = miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = self.complete_set)
    fft_map.apply_sigma_scaling()
    self.map_data_orig = fft_map.real_map_unpadded()
    if(self.use_complete_with):
      self.outlier_substitutes = self.get_outlier_substitutes(fmodel = fmodel)
    self.fmodel_dc = fmodel.deep_copy()
    self.f_calc   = self.fmodel_dc.f_calc()
    self.f_mask   = self.fmodel_dc.f_masks()[0]
    self.f_model  = self.fmodel_dc.f_model_no_scales()
    k_isotropic   = self.fmodel_dc.k_isotropic()
    k_anisotropic = self.fmodel_dc.k_anisotropic()
    self.fmodel_dc2 = mmtbx.f_model.manager(
      f_obs         = fmodel.f_obs(),
      r_free_flags  = fmodel.r_free_flags(),
      k_isotropic   = k_isotropic,
      k_anisotropic = k_anisotropic,
      f_calc        = self.f_model,
      f_mask        = self.f_mask.customized_copy(data = self.f_mask.data()*0))
    map_data = None
    for opt in [("f_calc", 5,self.shake_f_calc),
                ("f_mask", 5,self.shake_f_mask),
                ("f_model",10,self.shake_f_model),
                ("f_map",  10,self.shake_f_map)]:
      array, number_of_trials, run_status = opt
      if(not run_status): continue
      for it in xrange(number_of_trials):
        print "%7s: step %2d out of %2d" % (array, it, number_of_trials)
        complete_set = self.alter_complete_set_with_delta_fofc_outliers()
        mc = self.shake(complete_set=complete_set, array=array)
        if(self.use_complete_with):
          s = flex.random_bool(mc.size(), 0.8)
          mc = mc.select(s)
          mc = mc.complete_with(complete_set, scale=True)
          if(1): # XXX option to explore
            sel = flex.random_bool(mc.size(), 0.99)
            mc = mc.select(sel)
        fft_map = miller.fft_map(
          crystal_gridding     = crystal_gridding,
          fourier_coefficients = mc)
        fft_map.apply_sigma_scaling()
        m = fft_map.real_map_unpadded()
        if(map_data is None): map_data = m
        else:
          # critical to do both
          if(use_intersection):
            for i in [0,0.1,0.2,0.3,0.4,0.5]:
              maptbx.intersection(
                map_data_1 = m,
                map_data_2 = map_data,
                threshold  = i)
          map_data = (m+map_data)/2
    if(use_intersection):
      for i in xrange(3):
        maptbx.map_box_average(
          map_data   = map_data,
          cutoff     = 0.5,
          index_span = 1)
    self.map_data2 = map_data.deep_copy()
    if(use_intersection):
      # Very powerful in combination with NO HP modification.
      sd = map_data.sample_standard_deviation()
      map_data = map_data/sd
      maptbx.reset(data=map_data, substitute_value=0.0, less_than_threshold=0.0)
    sd = map_data.sample_standard_deviation()
    map_data = map_data/sd
    self.map_data = map_data
    self.map_coefficients = self.complete_set.structure_factors_from_map(
      map            = self.map_data,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)

  def shake(self, complete_set, array):
    assert array in ["f_mask", "f_calc", "f_model", "f_map", None]
    def get_mc(fm):
      return fm.electron_density_map(
        update_f_part1=False).map_coefficients(
          map_type     = self.map_type,
          isotropize   = True,
          fill_missing = False)
    if(self.shake_f_calc and array == "f_calc"):
      self.mean_phase_error=[0,10,30,40,60]
      f_calc_kick = self.call_run_kick_loop(map_coeffs=self.f_calc)
      fm = self.recreate_r_free_flags(fmodel = self.fmodel_dc)
      fm.update(f_calc = f_calc_kick)
      mc = get_mc(fm=fm)
    elif(self.shake_f_mask and array == "f_mask"):
      self.mean_phase_error=[0,10,30,40,60]
      f_mask_kick = self.call_run_kick_loop(map_coeffs=self.f_mask)
      fm = self.recreate_r_free_flags(fmodel = self.fmodel_dc)
      fm.update(f_mask=f_mask_kick)
      mc = get_mc(fm=fm)
    elif(self.shake_f_model and array == "f_model"):
      self.mean_phase_error=[0,5,10]
      f_model_kick = self.call_run_kick_loop(map_coeffs=self.f_model)
      fm = self.recreate_r_free_flags(fmodel = self.fmodel_dc2)
      fm.update(f_calc = f_model_kick)
      mc = get_mc(fm=fm)
    elif(self.shake_f_map and array == "f_map"):
      self.mean_phase_error=[0,5,10]
      mc = self.call_run_kick_loop(map_coeffs=complete_set)
    else:
      mc = complete_set
    return mc

  def alter_complete_set_with_delta_fofc_outliers(self):
    if(self.use_complete_with):
      rc = random.choice([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
      s = flex.random_bool(self.outlier_substitutes.outlier_ma.data().size(),rc)
      outl = self.outlier_substitutes.outlier_ma.select(s)
      complete_set2_, dummy = self.complete_set.lone_sets(outl)
      C = random.choice([0,1,2])
      if(C==0):
        complete_set = complete_set2_.complete_with(
          other=self.outlier_substitutes.outlier_substitute1, scale=True)
      elif(C==1):
        complete_set = complete_set2_.complete_with(
          other=self.outlier_substitutes.outlier_substitute2, scale=True)
      else:
        complete_set = self.complete_set
    else:
      complete_set = self.complete_set
    return complete_set

  def recreate_r_free_flags(self, fmodel):
    if(not self.update_r_free_flags): return fmodel
    rc = random.choice([0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
    r_free_flags = flex.random_bool(fmodel.f_obs().indices().size(), rc)
    fmodel._r_free_flags._data = r_free_flags
    return fmodel

  def get_outlier_substitutes(self, fmodel):
    threshold_value, delta_fofc = fmodel.top_largest_f_obs_f_model_differences(
      threshold_percent=self.delta_fofc_outliers_threshold)
    SELR = delta_fofc>threshold_value
    f_obs          = fmodel.f_obs().deep_copy()
    r_free_flags   = fmodel.r_free_flags().deep_copy()
    f_obs_outliers = f_obs.select(SELR)
    SELR_ = ~SELR
    fmodel_dc = mmtbx.f_model.manager(
      f_obs        = f_obs.select(SELR_),
      r_free_flags = r_free_flags.select(SELR_),
      xray_structure = fmodel.xray_structure)
    fmodel_dc.update_all_scales(update_f_part1_for="map", map_neg_cutoff=0)
    self.complete_set_2 = fmodel_dc.electron_density_map(
      update_f_part1 = True).map_coefficients(
        map_type     = self.map_type,
        isotropize   = True,
        fill_missing = True)
    self.complete_set, self.complete_set_2 = self.complete_set.common_sets(
      self.complete_set_2)
    fo = fmodel.f_obs().customized_copy(
      data = fmodel.f_obs().data()/(fmodel.k_isotropic()*fmodel.k_anisotropic()))
    outlier_substitute1, dummy = self.complete_set_2.common_sets(f_obs_outliers)
    outlier_substitute2, dummy = fo.common_sets(f_obs_outliers)
    outlier_substitute1, outlier_substitute2 = outlier_substitute1.common_sets(
      outlier_substitute2)
    assert outlier_substitute1.indices().all_eq(outlier_substitute2.indices())
    assert f_obs_outliers.size() == outlier_substitute1.size()
    outlier_substitute2 = outlier_substitute2.phase_transfer(
      phase_source=outlier_substitute1)
    return group_args(
      outlier_ma = f_obs_outliers,
      outlier_substitute1 = outlier_substitute1,
      outlier_substitute2 = outlier_substitute2)

  def convert_to_non_anomalous(self, fmodel):
    if(fmodel.f_obs().anomalous_flag()):
      f_obs        = fmodel.f_obs().average_bijvoet_mates()
      r_free_flags = fmodel.r_free_flags().average_bijvoet_mates()
      fmodel = mmtbx.f_model.manager(
        f_obs = f_obs,
        r_free_flags = r_free_flags,
        xray_structure = fmodel.xray_structure)
      fmodel.update_all_scales(update_f_part1_for="map")
    return fmodel

  def call_run_kick_loop(self, map_coeffs):
    map_coeff_data = self.run_kick_loop(map_coeffs = map_coeffs)
    return miller.set(
      crystal_symmetry = map_coeffs.crystal_symmetry(),
      indices          = map_coeffs.indices(),
      anomalous_flag   = False).array(data = map_coeff_data)

  def run_kick_loop(self, map_coeffs):
    map_coeff_data = None
    for kick in xrange(self.number_of_kicks):
      rc, ar = random.choice([(0.1,0.05),(0.2,0.05),
                              (0.3,0.04),(0.4,0.04),
                              (0.5,0.03),(0.6,0.03),
                              (0.7,0.02),(0.8,0.02),
                              (0.9,0.01),(1.0,0.01)])
      sel = flex.random_bool(map_coeffs.size(), rc)
      pr = random.choice(self.mean_phase_error)
      mc = map_coeffs.randomize_amplitude_and_phase(
        amplitude_error=ar, phase_error_deg=pr, selection=sel)
      if(map_coeff_data is None): map_coeff_data = mc.data()
      else:                       map_coeff_data = map_coeff_data + mc.data()
    return map_coeff_data/self.number_of_kicks

def fem(ko, crystal_gridding, mc_orig, fmodel):
  if(1): # logically seems better: why modify a maniulated (by FT) map!
    map_data = ko.map_data.deep_copy()
  else:
    fft_map = miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = ko.map_coefficients)
    fft_map.apply_sigma_scaling()
    map_data = fft_map.real_map_unpadded()
  map_dataK = map_data.deep_copy()
  if(0):
    o = maptbx.volume_scale(map = map_data, n_bins = 10000)
    fem = ko.complete_set.structure_factors_from_map(
      map            = o.map_data(),
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
  else:
    fem = maptbx.local_scale(
      map_data         = map_data,
      crystal_gridding = crystal_gridding,
      crystal_symmetry = fmodel.f_obs().crystal_symmetry(),
      miller_array     = ko.complete_set).map_coefficients
  fft_map = miller.fft_map(
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = mc_orig)
  fft_map.apply_sigma_scaling()
  map_orig = fft_map.real_map_unpadded()
  #
  fft_map = miller.fft_map(
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = fem)
  fft_map.apply_sigma_scaling()
  map_fem = fft_map.real_map_unpadded()
  cut_by_map = map_orig
  if(0): cut_by_map = ko.map_data
  maptbx.cut_by(kick=cut_by_map, fem=map_fem, cut_by_threshold=0.5)
  maptbx.cut_by(kick=map_dataK, fem=map_fem, cut_by_threshold=0.25)
  # these are alternatives to explore
  #oo = maptbx.omit(crystal_gridding=crystal_gridding, fmodel=fmodel)
  #maptbx.cut_by(kick=oo.map_result, fem=map_fem, cut_by_threshold=0.5)
  #
  #maptbx.cut_by(kick=map_orig, fem=map_fem)
  #maptbx.cut_by(kick=ko.map_data, fem=map_fem)
  #
  #maptbx.cut_by(kick=map_orig, fem=map_fem)
  #maptbx.cut_by(kick=map_dataK, fem=map_fem)
  return ko.complete_set.structure_factors_from_map(
    map            = map_fem,
    use_scale      = True,
    anomalous_flag = False,
    use_sg         = False)

class omit(object):
  def __init__(
        self,
        crystal_gridding,
        fmodel):
    b = boxes(
      n_real     = crystal_gridding.n_real(),
      unit_cell  = fmodel.f_obs().crystal_symmetry().unit_cell(),
      box_size_as_unit_cell_fraction = 0.11)
    #
    import mmtbx.f_model
    def get_map_data(fmodel, map_type):
      f_map = fmodel.electron_density_map(
        update_f_part1=False).map_coefficients(
          map_type     = map_type,
          isotropize   = True,
          fill_missing = False)
      fft_map = cctbx.miller.fft_map(
        crystal_gridding     = crystal_gridding,
        fourier_coefficients = f_map)
      fft_map.apply_volume_scaling()
      return fft_map.real_map_unpadded()
    f_map_data_orig = get_map_data(fmodel=fmodel, map_type="2mFo-DFc")
    f_model = fmodel.f_model_scaled_with_k1()
    fft_map = cctbx.miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = f_model)
    fft_map.apply_volume_scaling()
    f_model_map_data = fft_map.real_map_unpadded()
    zero_complex_ma = f_model.customized_copy(data = f_model.data()*0)
    self.map_result = flex.double(flex.grid(b.n_real))
    if 0: print "start1:",fmodel.r_work()
    fmodel_ = mmtbx.f_model.manager(
        f_obs        = fmodel.f_obs(),
        r_free_flags = fmodel.r_free_flags(),
        f_calc       = f_model,
        f_mask       = zero_complex_ma)
    if 0: print "strat2:",fmodel_.r_work()
    self.map_ave = None
    for s,e in zip(b.starts, b.ends):
      #print s,e
      box = copy(f_model_map_data, s, e)
      box.reshape(flex.grid(box.all()))
      # modify map in box
      #reset(data=box, substitute_value=0.0, less_than_threshold=0.1)
      #for i in xrange(3):
      #  map_box_average(
      #    map_data   = box,
      #    cutoff     = 1,
      #    index_span = 1)
      reset(data=box, substitute_value=0.0, less_than_threshold=9999)
      tmp = f_model_map_data.deep_copy()
      set_box(
        map_data_from = box,
        map_data_to   = tmp,
        start         = s,
        end           = e)
      f_model_ = f_model.structure_factors_from_map(
        map            = tmp,
        use_scale      = True,
        anomalous_flag = False,
        use_sg         = False)
      fmodel_ = mmtbx.f_model.manager(
        f_obs        = fmodel.f_obs(),
        r_free_flags = fmodel.r_free_flags(),
        f_calc       = f_model_,
        f_mask       = zero_complex_ma)
      print fmodel_.r_work()
      f_map_data = get_map_data(fmodel=fmodel_, map_type="mFo-DFc")
      #if(self.map_ave is None): self.map_ave = f_map_data
      #else: self.map_ave = self.map_ave + f_map_data
      box = copy(f_map_data, s, e)
      box.reshape(flex.grid(box.all()))
      #reset(data=box, substitute_value=0.0, less_than_threshold=0)
      set_box(
        map_data_from = box,
        map_data_to   = self.map_result,
        start         = s,
        end           = e)
    #
    #self.map_ave =  self.map_ave/len(b.starts)
    #self.map_coefficients = f_model.structure_factors_from_map(
    #  map            = self.map_ave,
    #  use_scale      = True,
    #  anomalous_flag = False,
    #  use_sg         = False)
   # sd = self.map_result.sample_standard_deviation()
   # self.map_result = self.map_result/sd
   #
   # for i in [0,0.1,0.2,0.3,0.4,0.5]:
   #   intersection(
   #     map_data_1 = f_map_data_orig,
   #     map_data_2 = self.map_result,
   #     threshold  = i)
    self.map_coefficients = f_model.structure_factors_from_map(
      map            = self.map_result,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
