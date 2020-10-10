
from simtbx.diffBragg.refiners.local_refiner import LocalRefiner
import numpy as np
from simtbx.diffBragg import utils
from copy import deepcopy


def local_refiner_from_parameters(refls, expt, params, miller_data=None):

  for model in ["crystal", "detector", "beam", "imageset"]:
    if not hasattr(expt, model):
      print("No %s in experiment, exiting. " % model)
      return

  n_trials = len(params.refiner.max_calls)
  if n_trials == 0:
    print("No trials available")
    return

  img_data = utils.image_data_from_expt(expt)
  if params.refiner.adu_per_photon is not None:
    img_data /= params.refiner.adu_per_photon

  background_mask = utils.load_mask(params.roi.background_mask)
  hotpix_mask = utils.load_mask(params.roi.hotpixel_mask)
  if hotpix_mask is not None:
      hotpix_mask = ~hotpix_mask

  # prepare ROI information
  panel_selection = [True]*len(refls)
  if params.roi.panels is not None:
    keeper_panels = utils.parse_panel_input_string(params.roi.panels)
    for i_pid, pid in enumerate(refls['panel']):
      panel_selection[i_pid] = pid in keeper_panels

  rois, pids, tilt_abc, selection_flags = utils.get_roi_background_and_selection_flags(
    refls, img_data, shoebox_sz=params.roi.shoebox_size,
    reject_edge_reflections=params.roi.reject_edge_reflections,
    reject_roi_with_hotpix=params.roi.reject_roi_with_hotpix,
    background_mask=background_mask, hotpix_mask=hotpix_mask,
    bg_thresh=params.roi.background_threshold,
    set_negative_bg_to_zero=params.roi.force_negative_background_to_zero,
    pad_for_background_estimation=params.roi.pad_shoebox_for_background_estimation)
  
  selection_flags = [sel1 and sel2 for sel1, sel2 in zip(selection_flags, panel_selection)]
  if not np.any(selection_flags):
      print("No spots for refinement")
      return None 

  
  # TODO consider panel selection combine here instead of select above

  nanoBragg_rois, xrel, yrel, roi_imgs = [], [], [], []
  for i_roi, (x1, x2, y1, y2) in enumerate(rois):
    nanoBragg_rois.append(((int(x1), int(x2)), (int(y1), int(y2))))
    yr, xr = np.indices((y2 - y1 + 1, x2 - x1 + 1))
    xrel.append(xr)
    yrel.append(yr)
    panel_id = pids[i_roi]
    roi_imgs.append(img_data[panel_id, y1:y2 + 1, x1:x2 + 1])

  # preprare arguments for refinement class instance
  UcellMan = utils.manager_from_crystal(expt.crystal)

  NCELLS_MASK = utils.get_ncells_mask_from_string(params.refiner.ncells_mask)
  if all(NCELLS_MASK):
    if not params.simulator.crystal.has_isotropic_ncells:
      print("WARNING: NCELLS mask specifies isotropic ncells, but simulator.crystal.has_isotropic_ncells is set to False")
    params.simulator.crystal.has_isotropic_ncells = True
  else:
    if params.simulator.crystal.has_isotropic_ncells:
      print("WARNING: NCELLS mask specifies anisotropic ncells, but simulator.crystal.has_isotropic_ncells is set to True")
    params.simulator.crystal.has_isotropic_ncells = False

  SIM = utils.simulator_from_expt_and_params(expt, params)
  if miller_data is not None:
    SIM.crystal.miller_array = miller_data.as_amplitude_array()
    SIM.update_Fhkl_tuple()
  shot_ucell_managers = {0: UcellMan}
  shot_rois = {0: rois}
  shot_nanoBragg_rois = {0: nanoBragg_rois}
  shot_roi_imgs = {0: roi_imgs}
  shot_spectra = {0: SIM.beam.spectrum}
  shot_crystal_models = {0: expt.crystal}
  shot_crystal_model_refs = {0: deepcopy(expt.crystal)}
  shot_xrel = {0: xrel}
  shot_yrel = {0: yrel}
  shot_abc_inits = {0: tilt_abc}
  shot_panel_ids = {0: pids}
  shot_originZ_init = {0:0}


  # determine number of parameters:
  if not any(NCELLS_MASK):
    n_ncells_param = 3
  elif all(NCELLS_MASK):
    n_ncells_param = 1
  else:
    n_ncells_param = 2

  nrot_params = 3
  n_unitcell_params = len(UcellMan.variables)
  n_spotscale_params = 1
  n_originZ_params = 1
  n_tilt_params = 3 * len(nanoBragg_rois)
  n_local_unknowns = nrot_params + n_unitcell_params + n_ncells_param + n_spotscale_params + n_originZ_params \
                     + n_tilt_params

  panel_group_from_id = {pid: 0 for pid in range(len(expt.detector))}
  if params.refiner.panel_group_file is not None:
    panel_group_from_id = utils.load_panel_group_file(params.refiner.panel_group_file)
  panel_groups = set(panel_group_from_id.values())
  n_panel_groups = len(panel_groups)
  panels_per_group = {group_id: [] for group_id in panel_groups}
  for pid in panel_group_from_id:
    group_id = panel_group_from_id[pid]
    panels_per_group[group_id].append(pid)
  
  n_spectra_params = 2 if params.refiner.refine_spectra is not None else 0
  n_panelRot_params = 3*n_panel_groups
  n_panelXYZ_params = 3*n_panel_groups
  n_global_params = n_spectra_params + n_panelRot_params + n_panelXYZ_params

  RUC = None
  x_init = None
  for i_macro_cyc in range(params.refiner.num_macro_cycles):
    for i_trial in range(n_trials):

      RUC = LocalRefiner(
        n_total_params=n_local_unknowns + n_global_params,
        n_local_params=n_local_unknowns,
        local_idx_start=0,
        shot_ucell_managers=shot_ucell_managers,
        shot_rois=shot_rois,
        shot_nanoBragg_rois=shot_nanoBragg_rois,
        shot_roi_imgs=shot_roi_imgs, shot_spectra=shot_spectra,
        shot_crystal_GTs=shot_crystal_model_refs, shot_crystal_models=shot_crystal_models,
        shot_xrel=shot_xrel, shot_yrel=shot_yrel, shot_abc_inits=shot_abc_inits,
        shot_asu=None,
        global_param_idx_start=n_local_unknowns,
        shot_panel_ids=shot_panel_ids,
        all_crystal_scales=None,
        global_ncells=False, global_ucell=False,
        shot_detector_distance_init=shot_originZ_init)

      if not params.refiner.only_predict_model:
        if params.refiner.refine_Bmatrix is not None:
          RUC.refine_Bmatrix = params.refiner.refine_Bmatrix[i_trial]

        if params.refiner.refine_Umatrix is not None:
          RUC.refine_Umatrix = params.refiner.refine_Umatrix[i_trial]

        if params.refiner.refine_ncells is not None:
          RUC.refine_ncells = params.refiner.refine_ncells[i_trial]

        if params.refiner.refine_bg is not None:
          RUC.refine_background_planes = params.refiner.refine_bg[i_trial]

        if params.refiner.refine_spot_scale is not None:
          RUC.refine_crystal_scale = params.refiner.refine_spot_scale[i_trial]

        if params.refiner.refine_spectra is not None:
          RUC.refine_spectra = params.refiner.refine_spectra[i_trial]

        if params.refiner.refine_detdist is not None:
          RUC.refine_detdist = params.refiner.refine_detdist[i_trial]

        if params.refiner.refine_panelZ is not None:
          RUC.refine_panelZ = params.refiner.refine_panelZ[i_trial]

        if params.refiner.refine_panelRotO is not None:
          RUC.refine_panelRotO = params.refiner.refine_panelRotO[i_trial]
        if params.refiner.refine_panelRotF is not None:
          RUC.refine_panelRotF = params.refiner.refine_panelRotF[i_trial]

        if params.refiner.refine_panelRotS is not None:
          RUC.refine_panelRotS = params.refiner.refine_panelRotS[i_trial]

        if params.refiner.refine_panelXY is not None:
          RUC.refine_panelXY = params.refiner.refine_panelXY[i_trial]

      if RUC.refine_detdist and RUC.refine_panelZ:
        raise ValueError("Cannot refine panelZ and detdist simultaneously")

      RUC.panel_group_from_id = panel_group_from_id

      RUC.panel_reference_from_id = {}
      for pid in panel_group_from_id:
        group_id = panel_group_from_id[pid]
        reference = expt.detector[panels_per_group[group_id][0]]
        RUC.panel_reference_from_id[pid] = reference.get_origin()
      refined_groups = []
      for i,pid in enumerate(pids):
        if selection_flags[i] and panel_group_from_id[pid] not in refined_groups:
          refined_groups.append(panel_group_from_id[pid])
      RUC.panel_groups_being_refined = refined_groups 

      RUC.panelRot_sigma = params.refiner.sensitivity.panelRotOFS
      RUC.panelX_sigma, RUC.panelY_sigma = params.refiner.sensitivity.panelXY
      RUC.panelZ_sigma = params.refiner.sensitivity.panelZ
      RUC.panelX_range = params.refiner.ranges.panel_X
      RUC.panelY_range = params.refiner.ranges.panel_Y
      RUC.panelZ_range = params.refiner.ranges.panel_Z
      RUC.panelRot_range = [[ang*np.pi/180 for ang in params.refiner.ranges.panel_rotO],
                            [ang*np.pi/180 for ang in params.refiner.ranges.panel_rotF],
                            [ang*np.pi/180 for ang in params.refiner.ranges.panel_rotS]]

      RUC.refine_Fcell = False
      RUC.max_calls = params.refiner.max_calls[i_trial]

      RUC.x_init = x_init
      RUC.only_pass_refined_x_to_lbfgs = False
      RUC.bg_extracted = False
      RUC.save_model = params.refiner.save_models

      if n_ncells_param == 2:
        INVERTED_NCELLS_MASK = [int(not mask_val) for mask_val in NCELLS_MASK]
        RUC.ncells_mask = INVERTED_NCELLS_MASK
      RUC.n_ncells_param = n_ncells_param
      RUC.recenter = True
      RUC.rescale_params = True
      RUC.rescale_fcell_by_resolution = params.refiner.rescale_fcell_by_resolution

      RUC.ignore_line_search_failed_step_at_lower_bound = True

      # INIT VALUES
      RUC.spot_scale_init = {0: params.refiner.init.spot_scale}  # self.spot_scale_init

      m_init = params.simulator.crystal.ncells_abc
      if n_ncells_param == 2:
        m_init = [m_init[i_ncell] for i_ncell in sorted(set(INVERTED_NCELLS_MASK))]
      RUC.m_init = {0: m_init}

      RUC.ucell_inits = {0: shot_ucell_managers[0].variables}
      if params.refiner.ranges.ucell_edge_percentage is not None:
        names = shot_ucell_managers[0].variable_names
        maxs, mins = [], []
        for i_n, n in enumerate(names):
          val = shot_ucell_managers[0].variables[i_n]
          if "Ang" in n:
            perc = params.refiner.ranges.ucell_edge_percentage*0.01
            valmin = val - val * perc
            valmax = val + val * perc
          else:
            deviation = params.refiner.ranges.ucell_angle_deviation
            valmin = val - deviation/2.
            valmax = val + deviation/2.
          mins.append(valmin)
          maxs.append(valmax)
        RUC.ucell_mins = {0: mins}
        RUC.ucell_maxs = {0: maxs}
        RUC.use_ucell_ranges = True

      # SIGMA VALUES
      RUC.rotX_sigma, RUC.rotY_sigma, RUC.rotZ_sigma = params.refiner.sensitivity.rotXYZ
      RUC.detector_distance_sigma = params.refiner.sensitivity.originZ
      RUC.ucell_sigmas = utils.unitcell_sigmas(UcellMan, params.refiner.sensitivity.unitcell)
      RUC.m_sigma = params.refiner.sensitivity.ncells_abc
      RUC.spot_scale_sigma = params.refiner.sensitivity.spot_scale
      RUC.a_sigma, RUC.b_sigma, RUC.c_sigma = params.refiner.sensitivity.tilt_abc
      RUC.fcell_sigma_scale = params.refiner.sensitivity.fcell

      RUC.print_all_missets = True 

      RUC.n_spectra_param = n_spectra_params
      RUC.spectra_coefficients_sigma = params.refiner.sensitivity.spectra_coefficients  # .01, .01
      RUC.spectra_coefficients_init = params.refiner.init.spectra_coefficients  # 0, 1
      RUC.lambda_coef_ranges = [params.refiner.ranges.spectra0, params.refiner.ranges.spectra1]
      RUC.detector_distance_range = params.refiner.ranges.originZ

      RUC.fcell_resolution_bin_Id = None
      RUC.compute_image_model_correlation = params.refiner.compute_image_model_correlation

      # plot things
      RUC.sigma_r = params.refiner.sigma_r / params.refiner.adu_per_photon
      RUC.trial_id = i_trial
      RUC.refine_rotZ = not params.refiner.fix_rotZ
      RUC.plot_images = params.refiner.plot.display
      RUC.plot_residuals = params.refiner.plot.as_residuals
      RUC.plot_stride = params.refiner.plot.iteration_stride
      RUC.setup_plots()
      RUC.log_fcells = True
      RUC.big_dump = params.refiner.big_dump
      RUC.request_diag_once = False
      RUC.S = SIM
      RUC.restart_file = params.refiner.io.restart_file
      RUC.has_pre_cached_roi_data = True
      RUC.trad_conv = True
      RUC.S.update_nanoBragg_instance('update_oversample_during_refinement', False)
      RUC.refine_gain_fac = False
      RUC.use_curvatures_threshold = params.refiner.use_curvatures_threshold
      if not params.refiner.curvatures:
        RUC.S.update_nanoBragg_instance('compute_curvatures', False)
      RUC.calc_curvatures = params.refiner.curvatures
      RUC.poisson_only = params.refiner.poissononly
      RUC.trad_conv_eps = params.refiner.tradeps
      RUC.verbose = params.refiner.verbose
      # TODO optional properties.. make this obvious
      RUC.FNAMES = []
      RUC.PROC_FNAMES = []
      RUC.PROC_IDX = []
      RUC.BBOX_IDX = []
      RUC.output_dir = params.refiner.io.output_dir
      RUC.run(setup_only=True)

      RUC.num_positive_curvatures = 0
      RUC.use_curvatures = params.refiner.start_with_curvatures
      RUC.hit_break_to_use_curvatures = False
      RUC.selection_flags = {0: selection_flags}
      RUC.record_model_predictions = params.refiner.record_xy_calc
      RUC.run(setup=False)
      if RUC.hit_break_to_use_curvatures:
        RUC.fix_params_with_negative_curvature = False
        RUC.num_positive_curvatures = 0
        RUC.use_curvatures = True
        RUC.run(setup=False)

      x_init = RUC.x
      if params.refiner.only_predict_model:
        if RUC.gnorm > 0:
          raise ValueError("Only predciting model, but the gradient is finite! This means the model changed, somethings wrong!")
      
      # TODO: CRITICAL ecause nanobragg requires a free_all we have to do it manually here to avoid MEM leak
      # RUC>S>D>free_all()
      #if x_init is not None:

  return RUC

