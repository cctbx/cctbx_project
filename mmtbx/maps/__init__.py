from cctbx.array_family import flex
import iotbx.phil
from mmtbx import utils as mmtbx_utils
import sys, os
from libtbx.utils import Sorry, date_and_time
from mmtbx import map_tools
import mmtbx
from cctbx import maptbx

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
    atom_selection_buffer = 3
      .type = float
      .expert_level=2
    apply_sigma_scaling = True
      .type = bool
      .expert_level = 1
    apply_volume_scaling = False
      .type = bool
      .expert_level = 2
    output_file_name = None
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

class compute_maps(object):

  def __init__(self,
               fmodel,
               params,
               mtz_dataset = None):
    # map coefficients
    self.mtz_dataset = mtz_dataset
    self.params = params
    for mcp in params.map_coefficients:
      if(mcp.map_type is not None):
        e_map_obj = fmodel.electron_density_map(
          fill_missing_f_obs = mcp.fill_missing_f_obs,
          fill_mode          = "dfmodel")
        if(not mcp.kicked):
          coeffs = e_map_obj.map_coefficients(map_type = mcp.map_type)
        else:
          km = map_tools.kick_map(
            fmodel            = e_map_obj.fmodel,
            map_type          = mcp.map_type,
            real_map          = True,
            real_map_unpadded = False,
            symmetry_flags    = maptbx.use_space_group_symmetry,
            average_maps      = False)
          coeffs = km.map_coeffs
        if(coeffs.anomalous_flag() and not
           mmtbx.map_names(mcp.map_type).anomalous):
          coeffs = coeffs.average_bijvoet_mates()
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
        raise RuntimeError("Not implemented.") # XXX add later

  def write_mtz_file(self, file_name, mtz_history_buffer = None):
    if(self.mtz_dataset is not None):
      if(mtz_history_buffer is None):
        mtz_history_buffer = flex.std_string()
      mtz_history_buffer.append(date_and_time())
      mtz_history_buffer.append("> file name: %s" % os.path.basename(file_name))
      mtz_object = self.mtz_dataset.mtz_object()
      mtz_object.add_history(mtz_history_buffer)
      mtz_object.write(file_name = file_name)
