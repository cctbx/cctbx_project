from cctbx.array_family import flex
import iotbx.phil
from mmtbx import utils as mmtbx_utils
import sys, os
from libtbx.utils import Sorry, date_and_time
from mmtbx import map_tools
import mmtbx
from cctbx import maptbx
from libtbx import adopt_init_args
from libtbx.str_utils import show_string
from libtbx.math_utils import ifloor, iceil

map_params_str ="""\
  map_coefficients
    .multiple = True
  {
    format = *mtz phs
      .type = choice(multi=True)
    mtz_label_amplitudes = None
      .type = str
      .expert_level=1
    mtz_label_phases = None
      .type = str
      .expert_level=1
    map_type = None
      .type = str
      .expert_level=1
    kicked = False
      .type = bool
      .expert_level=1
    fill_missing_f_obs = False
      .type = bool
      .expert_level=1
  }
  map
    .multiple = True
  {
    map_type = None
      .type = str
      .expert_level=1
    kicked = False
      .type = bool
      .expert_level=1
    fill_missing_f_obs = False
      .type = bool
      .expert_level=1
    grid_resolution_factor = 1/4.
      .type = float
      .expert_level=1
    region = *selection cell
      .type = choice
      .expert_level=1
      .short_caption=Map region
    atom_selection = None
      .type = str
      .expert_level=2
      .style = selection
    scale = *sigma volume
      .type = choice(multi=False)
      .expert_level = 1
    atom_selection_buffer = 3
      .type = float
      .expert_level=2
    file_name = None
      .type = str
  }
"""

def map_master_params():
  return iotbx.phil.parse(map_params_str, process_includes=False)

class map_coeffs_mtz_label_manager:

  def __init__(self, map_params):
    self._amplitudes = map_params.mtz_label_amplitudes
    self._phases = map_params.mtz_label_phases
    if (self._amplitudes is None or self._phases is None):
      raise Sorry("MTZ labels for map type %s undefined." % map_type)

  def amplitudes(self):
    return self._amplitudes

  def phases(self, root_label, anomalous_sign=None):
    assert anomalous_sign is None or not anomalous_sign
    return self._phases

class write_xplor_map_file(object):

  def __init__(self, params, coeffs, all_chain_proxies, xray_structure):
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
    if(self.params.region == "selection"):
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
    fft_map.as_xplor_map(
      file_name      = self.params.file_name,
      title_lines    = title_lines,
      gridding_first = gridding_first,
      gridding_last  = gridding_last)

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
      result = self.all_chain_proxies.selection(string =
        self.params.atom_selection).iselection()
    except KeyboardInterrupt: raise
    except Exception:
      raise Sorry('Invalid atom selection: %s' % self.params.atom_selection)
    if(result.size() == 0):
      raise Sorry('Empty atom selection: %s' % self.params.atom_selection)
    return result

class compute_maps(object):

  def __init__(self,
               fmodel,
               params,
               mtz_dataset = None,
               all_chain_proxies = None):
    adopt_init_args(self, locals())
    # map coefficients
    for mcp in params.map_coefficients:
      if(mcp.map_type is not None):
        coeffs = self.compute_map_coefficients(map_params = mcp)
        if("mtz" in mcp.format):
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
        if("phs" in mcp.format):
          raise RuntimeError("Not implemented") # XXX add later
    # xplor maps
    for mp in params.map:
      if(mp.map_type is not None):
        assert all_chain_proxies is not None
        coeffs = self.compute_map_coefficients(map_params = mp)
        write_xplor_map_file(params = mp, coeffs = coeffs,
          all_chain_proxies = self.all_chain_proxies,
          xray_structure = fmodel.xray_structure)

  def compute_map_coefficients(self, map_params):
    e_map_obj = self.fmodel.electron_density_map(
      fill_missing_f_obs = map_params.fill_missing_f_obs,
      fill_mode          = "dfmodel")
    if(not map_params.kicked):
      coeffs = e_map_obj.map_coefficients(map_type = map_params.map_type)
    else:
      coeffs = map_tools.kick_map(
        fmodel            = e_map_obj.fmodel,
        map_type          = map_params.map_type,
        real_map          = True,
        real_map_unpadded = False,
        symmetry_flags    = maptbx.use_space_group_symmetry,
        average_maps      = False).map_coeffs
    if(coeffs.anomalous_flag() and not
       mmtbx.map_names(map_params.map_type).anomalous):
      coeffs = coeffs.average_bijvoet_mates()
    return coeffs

  def write_mtz_file(self, file_name, mtz_history_buffer = None):
    if(self.mtz_dataset is not None):
      if(mtz_history_buffer is None):
        mtz_history_buffer = flex.std_string()
      mtz_history_buffer.append(date_and_time())
      mtz_history_buffer.append("> file name: %s" % os.path.basename(file_name))
      mtz_object = self.mtz_dataset.mtz_object()
      mtz_object.add_history(mtz_history_buffer)
      mtz_object.write(file_name = file_name)
