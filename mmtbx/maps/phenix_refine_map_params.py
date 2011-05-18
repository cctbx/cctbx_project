import iotbx.phil

# XXX temporary place for phenix.refine map calculation parameters.
# XXX It will be replaced with mmtbx/maps/__init__.py.

map_params_str ="""\
  map_format = *xplor
    .optional = True
    .type = choice(multi=True)
  map_coefficients_format = *mtz phs
    .optional = True
    .type = choice(multi=True)
  suppress = None
    .type = strings
    .help = List of mtz_label_amplitudes of maps to be suppressed. \
            Intended to selectively suppress computation and \
            writing of the standard maps.
    .expert_level = 1
  map
    .multiple = True
    .short_caption=Electron density map
    .style=noauto auto_align scrolled
  {
    mtz_label_amplitudes = None
      .type = str
      .short_caption=Amplitude label
    mtz_label_phases = None
      .type = str
      .short_caption=Phase label
    likelihood_weighted = None
      .type = bool
      .expert_level=1
    obs_factor = None
      .type = float
      .short_caption=Multiply Fobs by
    calc_factor = None
      .type = float
      .short_caption=Multiply Fcalc by
    kicked = False
      .type = bool
    fill_missing_f_obs_with_weighted_f_model = True
      .type = bool
  }
  map {
    mtz_label_amplitudes = 2FOFCWT
    mtz_label_phases = PH2FOFCWT
    likelihood_weighted = True
    obs_factor = 2
    calc_factor = 1
    kicked = False
    fill_missing_f_obs_with_weighted_f_model = True
  }
  map {
    mtz_label_amplitudes = FOFCWT
    mtz_label_phases = PHFOFCWT
    likelihood_weighted = True
    obs_factor = 1
    calc_factor = 1
    kicked = False
    fill_missing_f_obs_with_weighted_f_model = False
  }
  map {
    mtz_label_amplitudes = 2FOFCWT_no_fill
    mtz_label_phases = PH2FOFCWT_no_fill
    likelihood_weighted = True
    obs_factor = 2
    calc_factor = 1
    kicked = False
    fill_missing_f_obs_with_weighted_f_model = False
  }
  anomalous_difference_map
    .short_caption=Anomalous difference map
    .style = box auto_align
  {
    mtz_label_amplitudes = ANOM
      .type = str
      .short_caption=Amplitude label
    mtz_label_phases = PHANOM
      .type = str
      .short_caption=Phase label
  }
  reverse_scale = True
    .type = bool
    .expert_level = 2
    .help = Apply scales to Fobs and remove them from Fmodel (for maps only).
  grid_resolution_factor = 1/4
    .type = float
    .expert_level=1
  region = *selection cell
    .type = choice
    .expert_level=1
    .short_caption=Map region
  atom_selection = None
    .type = atom_selection
    .expert_level=2
  atom_selection_buffer = 3
    .type = float
    .expert_level=2
  apply_sigma_scaling = True
    .type = bool
    .expert_level = 1
  apply_volume_scaling = False
    .type = bool
    .expert_level = 2
"""

map_params = iotbx.phil.parse(map_params_str, process_includes=True)
