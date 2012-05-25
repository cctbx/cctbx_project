#Phil parameters required to process data from Stanford LCLS CXI instrument

def cxi_basic_start():
  from labelit import preferences
  from rstbx.command_line.index import special_defaults_for_new_horizons
  new_horizons_phil = preferences.RunTimePreferences()
  special_defaults_for_new_horizons( new_horizons_phil )

  common_arguments = [
          "distl.bins.verbose=True",
          "distl.minimum_spot_area=3",
          "distl.peripheral_margin=1",
          "distl.peak_intensity_maximum_factor=10000.", #avoids intensity filter
          "distl.compactness_filter=True",
  ]

  new_horizons_phil.merge_command_line(common_arguments)

  return new_horizons_phil

def cxi_versioned_extract(*args):

  # args is one or more lists of phil parameters, as would be passed
  # in through the command line; to be processed sequentially.

  working_phil = cxi_basic_start()

  for arg in args:
    working_phil.merge_command_line(arg)

  #distl_args = [a.object.as_str().strip() for a in working_phil.phil_scope.all_definitions()]

  cxi_version = working_phil.phil_scope.get("distl.detector_format_version"
                ).extract().detector_format_version

  if cxi_version in ["CXI 3.1","CXI 3.2"]:
    from spotfinder.applications.xfel import cxi_run3
    # run 3 tiles from distl.find_active_area, plus initial tile translations
    #   derived from a lysozyme powder pattern, before auxiliary adjustments from Bragg spots
    run_3_tiling_arguments = [
       "distl.detector_tiling=%s"%cxi_run3.run3_cxi_limits().as_string(),
    ] + cxi_run3.lysozyme_calibration()

    working_phil.merge_command_line(run_3_tiling_arguments)

    working_extract = working_phil.command_extractor

    corrected_auxiliary_translations = [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,
                               0,1,0,0,0,0,0,0,0,0,-1,3,-1,4,2,0,3,0,0,0,
                               0,0,0,0,1,3,-1,4,4,0,-1,1,0,4,3,-3,2,-1,1,3,
                               0,0,5,0,-1,2,0,0,1,0,0,0,-1,1,0,2,2,1,2,-1,
                               2,0,1,2,1,1,1,-6,0,0,2,3,0,2,2,3,2,-2,0,0,
                               0,0,0,0,0,0,0,0,-2,1,0,0,0,2,0,0,0,0,0,0,]
    from scitbx.array_family import flex
    total_tile_translations = flex.int(
      [int(a) for a in working_extract.distl.tile_translations]
      )+flex.int(corrected_auxiliary_translations)

    TT = list(total_tile_translations)
    working_extract.distl.tile_translations = TT

    return working_extract
  elif cxi_version in ["CXI 4.1"]:
    working_extract = working_phil.command_extractor

    corrected_auxiliary_translations = [
                               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,-1,
                               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,
                               0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                               0,0,0,0,0,0,0,0]

    from scitbx.array_family import flex
    total_tile_translations = flex.int(corrected_auxiliary_translations)

    TT = list(total_tile_translations)
    working_extract.distl.tile_translations = TT
    #working_extract.distl.quad_translations = [0,0,0,0,0,0,0,0]
    #working_extract.distl.quad_translations = [5,-5,6,-9,-5,-2,0,-8]
    working_extract.distl.quad_translations = [6,-2,7,-4,-4,1,1,-5]
    return working_extract

  elif cxi_version in ["CXI 5.1"]:
    working_extract = working_phil.command_extractor

    # The auxiliary translations are modified with respect to CXI 4.1.  If the
    # SLAC-provided metrology were to be trusted, this would be be all
    # zeros?
    corrected_auxiliary_translations = [
                               1,-1,1,0,0,-1,0,-1,2,0,0,0,3,1,2,0,2,-1,2,-1,
                               0,0,1,0,1,-2,1,-3,0,-1,0,-1,-2,1,0,1,-1,0,0,0,
                               0,2,-1,1,0,2,-1,2,0,0,0,0,0,0,1,1,1,0,0,0,
                               -1,0,-1,0,1,1,0,0,0,1,0,1,0,-1,0,0,1,-2,0,0,
                               -1,-2,-1,-1,0,0,-1,-1,0,0,0,1,-1,-1,0,0,1,-2,1,-1,
                               1,-1,0,0,0,-2,0,-3,1,-4,1,-4,-2,0,-1,-1,0,0,-1,-1,
                               -1,1,-1,0,0,1,0,1]

    from scitbx.array_family import flex
    total_tile_translations = flex.int(corrected_auxiliary_translations)

    TT = list(total_tile_translations)
    working_extract.distl.tile_translations = TT

    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y
    working_extract.distl.quad_translations = [-3,-1,-1,-5,-13,2,-7,-4]
    return working_extract

  else:
    return working_phil.command_extractor
