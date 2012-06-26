import mmtbx.utils
from libtbx.utils import Sorry, date_and_time
import iotbx.phil
from libtbx import adopt_init_args
from libtbx.str_utils import show_string
from libtbx.math_utils import ifloor, iceil
import os
import sys
from scitbx.array_family import flex

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
      .short_caption = PDB file
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
      .style = bold directory noauto
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
    fmodel_data_file_format = mtz cns
      .optional=True
      .type=choice
      .help=Write Fobs, Fmodel, various scales and more to MTZ or CNS file
  }
  scattering_table = wk1995  it1992  *n_gaussian  neutron
    .type = choice
    .help = Choices of scattering table for structure factors calculations
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

def map_coefficients_from_fmodel (fmodel, params,
    post_processing_callback=None):
  from mmtbx import map_tools
  import mmtbx
  from cctbx import miller
  mnm = mmtbx.map_names(map_name_string = params.map_type)
  if(mnm.k==0 and abs(mnm.n)==1):
    return compute_f_calc(fmodel, params)
  #XXXsave_k_part, save_b_part = None, None
  #XXXif(mnm.k is not None and abs(mnm.k) == abs(mnm.n) and fmodel.k_part()!=0):
  #XXX  save_k_part = fmodel.k_part()
  #XXX  save_b_part = fmodel.b_part()
  #XXX  fmodel.update_core(k_part=0, b_part=0)
  e_map_obj = fmodel.electron_density_map()
  coeffs = None
  if(not params.kicked):
    coeffs = e_map_obj.map_coefficients(
      map_type           = params.map_type,
      acentrics_scale    = params.acentrics_scale,
      centrics_pre_scale = params.centrics_pre_scale,
      ncs_average=getattr(params, "ncs_average", False),
      post_processing_callback=post_processing_callback)
    if (coeffs is None) : return None
    if(coeffs.anomalous_flag() and not
       mmtbx.map_names(params.map_type).anomalous):
      coeffs = coeffs.average_bijvoet_mates()
    if(params.isotropize and e_map_obj.mch is not None):
      isotropize_helper = e_map_obj.fmodel.map_calculation_helper().isotropize_helper
      if(isotropize_helper is not None):
        isc = isotropize_helper.iso_scale.data()
      else: isc = 1
      coeffs = miller.set(
        crystal_symmetry=coeffs.crystal_symmetry(),
        indices = coeffs.indices(),
        anomalous_flag=False).array(data=coeffs.data()*isc)
    if(params.sharpening):
      from mmtbx import map_tools
      coeffs, b_sharp = map_tools.sharp_map(
        sites_frac = fmodel.xray_structure.sites_frac(),
        map_coeffs = coeffs,
        b_sharp    = params.sharpening_b_factor)
  else:
    if fmodel.__class__.__name__ == "twin_model_manager" :
      raise Sorry("Kicked maps are not supported when twinning is present.  "+
        "You can disable the automatic twin law detection by setting the "+
        "parameter maps.skip_twin_detection to True (or check the "+
        "corresponding box in the Phenix GUI).")
    if(params.map_type.count("anom")==0):
      coeffs = map_tools.kick_map(
        fmodel     = e_map_obj.fmodel,
        map_type   = params.map_type,
        isotropize = params.isotropize,
        sharp      = params.sharpening).map_coeffs
  # XXX need to figure out why this happens
  if (coeffs is None) :
    raise RuntimeError(("Map coefficient generation failed (map_type=%s, "
      "kicked=%s, sharpening=%s, isotropize=%s, anomalous=%s.") %
        (params.map_type, params.kicked, params.sharpening, params.isotropize,
         fmodel.f_obs().anomalous_flag()))
  # XXX is this redundant?
  if(coeffs.anomalous_flag() and not
     mmtbx.map_names(params.map_type).anomalous):
    coeffs = coeffs.average_bijvoet_mates()
  if(params.exclude_free_r_reflections):
    r_free_flags = e_map_obj.fmodel.r_free_flags()
    if (r_free_flags.anomalous_flag()) :
      r_free_flags = r_free_flags.average_bijvoet_mates()
    coeffs = coeffs.select(~r_free_flags.data())
  #XXXif(mnm.k is not None and abs(mnm.k) == abs(mnm.n) and save_k_part is not None):
  #XXX  fmodel.update_core(k_part=save_k_part, b_part=save_b_part)
  if(params.fill_missing_f_obs):
    from mmtbx import map_tools
    coeffs = map_tools.fill_missing_f_obs(coeffs, fmodel)
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
               log=sys.stdout):
    assert ((post_processing_callback is None) or
            (hasattr(post_processing_callback, "__call__")))
    self.mtz_dataset = mtz_dataset
    coeffs = None
    self.map_coeffs = []
    for mcp in params:
      if(mcp.map_type is not None):
        # XXX
        if(fmodel.__class__.__name__ == "twin_model_manager") :
          if (mcp.map_type == "anomalous") :
            print >> log, "Anomalous maps not supported when a twin law is used."
            continue
          elif (mcp.isotropize) :
            mcp.isotropize = False
        # XXX
        coeffs = map_coefficients_from_fmodel(fmodel = fmodel,
          params = mcp,
          post_processing_callback=post_processing_callback)
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

  def write_mtz_file(self, file_name, mtz_history_buffer = None):
    from cctbx.array_family import flex
    if(self.mtz_dataset is not None):
      if(mtz_history_buffer is None):
        mtz_history_buffer = flex.std_string()
      mtz_history_buffer.append(date_and_time())
      mtz_history_buffer.append("> file name: %s" % os.path.basename(file_name))
      mtz_object = self.mtz_dataset.mtz_object()
      mtz_object.add_history(mtz_history_buffer)
      mtz_object.write(file_name = file_name)
      return True
    return False
