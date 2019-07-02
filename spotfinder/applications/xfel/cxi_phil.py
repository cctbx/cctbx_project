# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-

from __future__ import absolute_import, division, print_function
#Phil parameters required to process data from Stanford LCLS CXI instrument

def cxi_basic_start():
  try:
    from labelit import preferences
    from rstbx.command_line.index import special_defaults_for_new_horizons
  except Exception:
    # option to view images or open pickle files without any LABELIT dependency
    from rstbx.phil import preferences
    def special_defaults_for_new_horizons(phil_scope):
      # for integration, do not want 2x2 binning
      phil_scope.merge_command_line(["distl_permit_binning=False"])

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
  import copy

  # args is one or more lists of phil parameters, as would be passed
  # in through the command line; to be processed sequentially.

  working_phil = cxi_basic_start()

  #for arg in args:
  for arg in copy.deepcopy(args):
    working_phil.merge_command_line(arg)
  working_extract = working_phil.command_extractor

  if working_extract.distl.tile_translations is not None and \
     working_extract.distl.quad_translations is not None:
    return working_extract

  # Get the working tile translations for the given detector format version
  # from previous experiments
  legacy_extract = cxi_versioned_extract_detail(args)

  if working_extract.distl.tile_translations is None:
    working_extract.distl.tile_translations = legacy_extract.distl.tile_translations

  if working_extract.distl.quad_translations is None:
    working_extract.distl.quad_translations = legacy_extract.distl.quad_translations

  return working_extract

def cxi_versioned_extract_detail(args):

  # args is one or more lists of phil parameters, as would be passed
  # in through the command line; to be processed sequentially.

  working_phil = cxi_basic_start()

  for arg in args:
    working_phil.merge_command_line(arg)

  #distl_args = [a.object.as_str().strip() for a in working_phil.phil_scope.all_definitions()]

  cxi_version = working_phil.phil_scope.get("distl.detector_format_version"
                ).extract().detector_format_version

  print("cxi_versioned_extract()::cxi_version:", cxi_version)

  if cxi_version in ["CXI 3.1","CXI 3.2"]:
    working_extract = working_phil.command_extractor

    # Three sensors were disabled at this time.  The 2D-translations
    # will have (32 - 3) * 2 * 2 = 116 components.
    corrected_auxiliary_translations = [
        1,   1,   1,   1,   1,   2,   1,   2,   0,   0,   0,   0,   1,  -1,   1,  -1,
        1,  -2,   1,  -2,   0,  -2,   0,  -2,  -1,   1,  -1,   1,   1,   0,   1,   0,
        2,  -1,   2,  -1,  -3,   1,  -3,   1,  -4,   0,  -4,   0,   4,   0,   4,   0,
       -7,   4,  -7,   4,   2,  -1,   2,  -1,   1,  -2,   1,  -2,  -1,   1,  -1,   1,
       -1,   1,  -1,   1,   0,  -1,   0,  -1,   2,  -2,   2,  -2,   0,   1,   0,   1,
        0,   2,   0,   2,  -1,  -2,  -1,  -2,  -1,  -3,  -1,  -3,   1,   1,   1,   1,
        0,   1,   0,   1,   0,  -1,   0,  -1,   1,  -2,   1,  -2,  -2,  -1,  -2,  -1,
        2,   2,   2,   2]

    from scitbx.array_family import flex
    total_tile_translations = flex.int(corrected_auxiliary_translations)

    TT = list(total_tile_translations)
    working_extract.distl.tile_translations = TT
    working_extract.distl.quad_translations = [-19, 15, 21, 14, -25, -20, 20, -15]

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
    #added NKS 8/19/14
    #TT = [0]*128
    working_extract.distl.tile_translations = TT

    print("IN CXI %>! WITH ",working_extract.distl.tile_translations)

    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y
    working_extract.distl.quad_translations = [-3,-1,-1,-5,-13,2,-7,-4]
    return working_extract

  elif cxi_version in ["CXI 6.1"]:
    working_extract = working_phil.command_extractor
    corrected_auxiliary_translations = [
       2,  1,  1,  1,  1,  3, -1,  2,  3,  1,
       1,  0,  5,  2,  4,  1,  2, -1,  2,  0,
       2,  1,  2,  0, -1, -2, -1, -2, -1,  0,
      -2,  1, -1,  0,  0,  1,  1,  0,  1,  1,
      -1,  0, -1,  0, -1,  0, -1,  0,  0,  0,
      -1,  0,  1,  0,  1,  0,  0,  1,  1,  2,
      -1,  1,  0,  2, -1,  1,  0,  1,  0,  0,
       1,  0, -2,  0, -1,  1, -2, -1, -2,  1,
      -2,  0, -1,  1, -3, -1, -3,  0,  0,  1,
       0,  1,  0,  0,  1,  0,  2,  0,  2, -1,
       1,  0,  0, -1,  0,  0,  1, -2,  1, -1,
       2, -2, -1,  0,  0, -1, -2,  1,  0,  0,
      -1,  0, -1, -1,  1,  1,  1, -1]

    from scitbx.array_family import flex
    total_tile_translations = flex.int(corrected_auxiliary_translations)

    TT = list(total_tile_translations)
    working_extract.distl.tile_translations = TT

    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y
    working_extract.distl.quad_translations = [0,7,13,8,-8,0,11,-3]
    return working_extract

  elif cxi_version in ["CXI 7.1"]:
    working_extract = working_phil.command_extractor
    corrected_auxiliary_translations = [
       2,  1,  1,  1,  1,  2,  0,  2,  3,  0,
       0,  0,  4,  1,  3,  1,  1, -2,  2, -1,
      -1, -1,  1, -1, -1, -2, -1, -1, -1,  0,
      -1,  1, -2,  0,  0,  1,  0,  0,  1,  1,
      -1,  1, -2,  1, -1,  0, -1,  1, -1,  1,
      -1,  1, -2,  1, -1,  1,  1,  1,  0,  2,
      -1,  1,  0,  1,  0,  1,  0,  1,  0,  0,
       0,  1, -1,  0,  0,  1,  0, -1, -1,  1,
      -1,  0, -1,  1, -1,  0, -2,  0,  0,  1,
       0,  1, -1,  0,  1,  0,  2, -1,  2, -1,
       1,  0,  0,  0,  1, -1,  1, -2,  2, -2,
       2, -2,  0, -1,  1, -1,  2,  1,  0,  0,
       0,  0,  0, -2,  1,  1,  1,  0]

    from scitbx.array_family import flex
    total_tile_translations = flex.int(corrected_auxiliary_translations)

    TT = list(total_tile_translations)
    working_extract.distl.tile_translations = TT

    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y
    working_extract.distl.quad_translations = [2,-6,3,-6,-7,0,-1,-4]
    return working_extract

  elif cxi_version in ["CXI 7.d"]:
    working_extract = working_phil.command_extractor

    corrected_auxiliary_translations = [
       1,  1,  1,  0,  0,  0,  0,  0, -2, -1,
      -1, -1, -1, -3,  0, -3,  2,  3,  2,  1,
      -1,  4, -1,  2,  0,  1,  1,  1,  0,  2,
      -1,  2,  0,  0,  0,  0,  0,  0,  0,  0,
      -3, -2, -3, -1, -2, -1, -2, -1,  3, -3,
       3, -3,  5, -2,  4, -2,  2, -1,  2, -1,
       2, -2,  2, -1,  1,  1,  1,  0,  0,  0,
       0, -1, -1,  2, -2,  2, -2,  0, -1,  0,
      -1, -2, -1, -1, -3, -1, -3, -1,  3,  0,
       2,  0,  2,  0,  2,  0,  1,  0,  0,  0,
       0,  0, -1,  0,  0, -1,  0, -1,  1, -2,
       1, -2, -6,  0, -6,  0, -5,  1, -5,  1,
      -5, -2, -5, -3, -4, -2, -4, -2]

    from scitbx.array_family import flex
    total_tile_translations = flex.int(corrected_auxiliary_translations)

    TT = list(total_tile_translations)
    working_extract.distl.tile_translations = TT

    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y
    working_extract.distl.quad_translations = [-6,-2,9,2,-14,-8,7,-12]
    return working_extract

  elif cxi_version in ["XPP 7.1"]:
    working_extract = working_phil.command_extractor

    corrected_auxiliary_translations = [
       0,  0,  0,  0,  0,  0,  0,  0,
       0, -2,  0, -2, -1, -1, -1, -1,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0, -2, -2, -2, -2,
      -1, -1, -1, -1, -1, -1, -1, -1,
       0,  0,  0,  0, -2, -1, -2, -1,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0, -1,  0, -1,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0, -1, -1, -1, -1,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0]

    L748_corrections_post103 = [
       0,  0,  0,  0,  0,  0,  0,  0,
       0, -1,  0, -1, -1,  0, -1,  0,
       1, -1,  1, -1,  0,  0,  0,  0,
       0,  0,  0,  0,  2,  1,  2,  1,

       0,  0,  0,  0, -1,  0, -1,  0,
       3,  0,  3,  0,  0,  1,  0,  1,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  1, -1,  1, -1,

       0,  0,  0,  0,  0,  1,  0,  1,
       1,  3,  1,  3,  1,  3,  1,  3,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  1,  1,  1,  1,

       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  2,  0,  2,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0, -2,  0, -2,  0,  0,  0,  0,
    ] # for 110 mm

    from scitbx.array_family import flex
    total_tile_translations = flex.int(corrected_auxiliary_translations) - \
                              flex.int(L748_corrections_post103)
    TT = list(total_tile_translations)
    working_extract.distl.tile_translations = TT

    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y.  For the
    # XPP CSPAD, this is effectively correcting for the beam center.
    working_extract.distl.quad_translations = [-3, -21,
                                               -3, -22,
                                               -3, -21,
                                                2, -21]
    return working_extract

  elif cxi_version in ["XPP 8.1"]:
    working_extract = working_phil.command_extractor

    corrected_auxiliary_translations = [
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0]

    from scitbx.array_family import flex
    total_tile_translations = flex.int(corrected_auxiliary_translations)

    TT = list(total_tile_translations)
    working_extract.distl.tile_translations = TT

    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y.  For the
    # XPP CSPAD, this is effectively correcting for the beam center.
    working_extract.distl.quad_translations = [0, 1,
                                               0, 1,
                                               0, 1,
                                               0, 1]
    return working_extract


  elif cxi_version in ["XPP 9.1"]:
    working_extract = working_phil.command_extractor

    # metrology from trial 14, runs 40, 41, xppe0314
    corrected_auxiliary_translations = [
       1,  2,  1,  2,  1,  2,  1,  2,
       2,  1,  2,  1,  2,  2,  2,  2,
      -1, -1, -1, -1,  0,  0,  0,  0,
      -2,  1, -2,  1, -1,  2, -1,  2,

       0, -1,  0, -1,  0,  0,  0,  0,
      -2, -1, -2, -1, -2, -1, -2, -1,
      -2,  0, -2,  0, -3, -1, -3, -1,
       0, -2,  0, -2,  0, -1,  0, -1,

       1, -1,  1, -1,  1, -1,  1, -1,
       0, -3,  0, -3,  1, -2,  1, -2,
      -2, -3, -2, -3,  0,  0,  0,  0,
      -3,  0, -3,  0,  0, -1,  0, -1,

      -1,  0, -1,  0,  0,  0,  0,  0,
      -3,  1, -3,  1, -3,  1, -3,  1,
      -1,  2, -1,  2, -4,  2, -4,  2,
       0,  0,  0,  0,  0, -1,  0, -1
    ]

    from scitbx.array_family import flex
    total_tile_translations = flex.int(corrected_auxiliary_translations)

    TT = list(total_tile_translations)
    working_extract.distl.tile_translations = TT

    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y.  For the
    # XPP CSPAD, this is effectively correcting for the beam center.
    # analysis by hand from run 40+41, xppe0314.
    #working_extract.distl.quad_translations = [3, -1,
    #                                           6, -1,
    #                                           4, -1,
    #                                           4, -5]
    # account for a change in beam center.  offset 7 pixels in the x direction.
    # valid for runs 124, 126
    working_extract.distl.quad_translations = [3,  6,
                                               6,  6,
                                               4,  6,
                                               4,  1]

    return working_extract

  elif cxi_version in ["XPP 11.1"]:
    working_extract = working_phil.command_extractor

    # metrology from trial 14, runs 40, 41, xppe0314
    corrected_auxiliary_translations = [
       1,  2,  1,  2,  1,  2,  1,  2,
       2,  1,  2,  1,  2,  2,  2,  2,
      -1, -1, -1, -1,  0,  0,  0,  0,
      -2,  1, -2,  1, -1,  2, -1,  2,

       0, -1,  0, -1,  0,  0,  0,  0,
      -2, -1, -2, -1, -2, -1, -2, -1,
      -2,  0, -2,  0, -3, -1, -3, -1,
       0, -2,  0, -2,  0, -1,  0, -1,

       1, -1,  1, -1,  1, -1,  1, -1,
       0, -3,  0, -3,  1, -2,  1, -2,
      -2, -3, -2, -3,  0,  0,  0,  0,
      -3,  0, -3,  0,  0, -1,  0, -1,

      -1,  0, -1,  0,  0,  0,  0,  0,
      -3,  1, -3,  1, -3,  1, -3,  1,
      -1,  2, -1,  2, -4,  2, -4,  2,
       0,  0,  0,  0,  0, -1,  0, -1
    ]

    from scitbx.array_family import flex
    total_tile_translations = flex.int(corrected_auxiliary_translations)

    TT = list(total_tile_translations)
    working_extract.distl.tile_translations = TT

    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y.  For the
    # XPP CSPAD, this is effectively correcting for the beam center.
    # analysis by hand from runs 35-38, trial 1, xpph9015.
    working_extract.distl.quad_translations = [1, 10,
                                               2,  9,
                                              -1,  9,
                                              -1,  4]

    return working_extract



  elif cxi_version in ["XPP 7.marccd"]:
    working_extract = working_phil.command_extractor
    working_extract.distl.quad_translations = None
    working_extract.distl.tile_translations = [0, 0]
    return working_extract

  elif cxi_version in ["CXI 8.1inheritedfrom7.1"]:
    working_extract = working_phil.command_extractor
    corrected_auxiliary_translations = [
       2,  1,  1,  1,  1,  2,  0,  2,  3,  0,
       0,  0,  4,  1,  3,  1,  1, -2,  2, -1,
      -1, -1,  1, -1, -1, -2, -1, -1, -1,  0,
      -1,  1, -2,  0,  0,  1,  0,  0,  1,  1,
      -1,  1, -2,  1, -1,  0, -1,  1, -1,  1,
      -1,  1, -2,  1, -1,  1,  1,  1,  0,  2,
      -1,  1,  0,  1,  0,  1,  0,  1,  0,  0,
       0,  1, -1,  0,  0,  1,  0, -1, -1,  1,
      -1,  0, -1,  1, -1,  0, -2,  0,  0,  1,
       0,  1, -1,  0,  1,  0,  2, -1,  2, -1,
       1,  0,  0,  0,  1, -1,  1, -2,  2, -2,
       2, -2,  0, -1,  1, -1,  2,  1,  0,  0,
       0,  0,  0, -2,  1,  1,  1,  0]

    from scitbx.array_family import flex
    total_tile_translations = flex.int(corrected_auxiliary_translations)

    TT = list(total_tile_translations)
    working_extract.distl.tile_translations = TT

    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y
    working_extract.distl.quad_translations = [-22,-4,-14,11,0,8,2,-5]
    return working_extract

  elif cxi_version in ["CXI 8.1"]: # developed against "run 4 optical metrology"
    working_extract = working_phil.command_extractor
    corrected_auxiliary_translations = [
      -2, -3, -2, -3, -3, -3, -3, -3,
      -4, -4, -4, -4, -2, -4, -2, -4,
      -2, -2, -2, -2,  0,  0, -1,  1,
      -4, -2, -4, -2, -3,  0, -3,  0,
      -4,  2, -4,  2, -4,  2, -4,  2,
      -5,  3, -5,  3, -5,  3, -5,  3,
      -2,  3, -2,  3,  3,  1,  3,  1,
      -1,  2, -1,  2, -2,  1, -2,  1,
       3,  3,  3,  3,  2,  3,  2,  3,
       2,  3,  2,  3,  2,  2,  2,  2,
       2,  4,  2,  4,  0, -1, -1,  0,
       3,  3,  3,  3,  2,  2,  2,  2,
       4, -3,  4, -3,  4, -3,  4, -3,
       3, -3,  4, -3,  3, -3,  4, -3,
      -2, -1, -2, -2, -3,  2, -3,  2,
       1, -3,  1, -3,  2, -2,  2, -2,
    ]

    from scitbx.array_family import flex
    total_tile_translations = flex.int(corrected_auxiliary_translations)

    TT = list(total_tile_translations)
    working_extract.distl.tile_translations = TT

    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y
    # Note, these numbers are valid for the CSPAD detector at roughly 100 mm distance
    working_extract.distl.quad_translations = [-14,7,-4,7,-22,0,-12,-4]

    # Calibrated at 201
    #working_extract.distl.quad_translations = [-14,7,-4,6,-23,0,-13,-5]

    # Due to a slight angle in the rail, if the detector is at the back of its
    # stage, around 548 mm, use these numbers instead.  They represent a change
    # of beam center due to a small translation at that distance
    #working_extract.distl.quad_translations = [-10,5,0,5,-18,-2,-8,-6]

    # Recalibration of 548 mm vs run 111 of lysozyme calibration dataset
    #working_extract.distl.quad_translations = [-11,4,-2,6,-18,-3,-9,-4]


    return working_extract

  elif cxi_version in ["CXI 8.2"]:
    working_extract = working_phil.command_extractor

    corrected_auxiliary_translations = [
       1,  7,  1,  7, -1,  0, -1,  0,
       5,  4,  5,  4,  1,  4,  1,  4,
       6, -2,  6, -2,  5,  3,  5,  3,
       5,  0,  5,  0,  0,  3,  0,  3,
       4,  0,  4,  0,  1,  0,  1,  0,
       3, -7,  3, -7,  3, -1,  3, -1,
      -1, -5, -1, -5,  6, -4,  6, -4,
       2, -3,  2, -3,  2, -4,  2, -4,
      -1, -5, -1, -5,  0,  0,  0,  0,
     -10,  0,-10,  0, -1, -1, -1, -1,
      -5,  8, -5,  8, -8,  2, -8,  2,
      -1,  2, -1,  2, -2,  1, -2,  1,
      -7,  0, -7,  0, -1,  0, -1,  0,
      -6,  4, -6,  4, -4, -6, -4, -6,
      -3,  2, -3,  2, -8,  3, -8,  3,
      -4,  6, -4,  6, -4,  2, -4,  2]

    LC06_mark10_001_corrections = [
       0,  0,  0,  0,  0, -1,  0, -1,
      -1,  1, -1,  1, -1,  1, -1,  1,
       0,  0,  0,  0,  0,  0,  0,  0,
       1,  2,  1,  2,  1,  0,  1,  0,

      -1,  0, -1,  0, -1,  0, -1,  0,
       0,  0,  0,  0,  0,  1,  0,  1,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  1,  0,  1,

       1,  2,  1,  2,  0,  1,  0,  1,
       0,  1,  0,  1,  0,  2,  0,  2,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  1,  0,  1,

       1, -1,  1, -1,  1, -1,  1, -1,
       0, -2,  0, -2,  1, -2,  1, -2,
       0,  0,  0,  0,  0,  0,  0,  0,
      -1,  2, -1,  2, -1,  0, -1,  0] # for 79 mm

    LC06_mark10_101_corrections = [
       0,  0,  0,  0,  0, -1,  0, -1,
      -1,  0, -1,  0, -0,  0, -0,  0,
       0,  1,  0,  1,  2,  2,  2,  2,
       0,  1,  0,  1,  0,  0,  0,  0,

      -1,  0, -1,  0, -1,  0, -1,  0,
      -1,  0, -1,  0, -1,  1, -1,  1,
      -1,  0, -1,  0,  0,  1,  0,  1,
      -1,  0, -1,  0,  0,  0,  0,  0,

       0,  2,  0,  2,  0,  1,  0,  1,
       0,  1,  0,  1,  0,  3,  0,  3,
       1,  2,  1,  2,  1,  3,  1,  3,
       0,  1,  0,  1,  0,  1,  0,  1,

       1, -1,  1, -1,  1, -1,  1, -1,
       0, -2,  0, -2,  2, -2,  2, -2,
       0, -1,  0, -1,  0, -2,  0, -2,
      -0,  0, -0,  0, -0,  0, -0,  0] # for 159 mm

    LB67_corrections_post009 = [
       1,  0,  1,  0,  0,  0,  0,  0,
       1,  1,  1,  1,  1,  1,  1,  1,
       2,  3,  2,  3,  2,  4,  2,  4,
       2,  2,  2,  2,  1,  0,  1,  0,

      -1,  0, -1,  0,  0,  0,  0,  0,
       0, -1,  0, -1,  0,  0,  0,  0,
       0, -1,  0, -1,  1,  0,  1,  0,
       0,  0,  0,  0,  0,  0,  0,  0,

       0,  1,  0,  1, -1,  0, -1,  0,
      -1,  1, -1,  1, -1,  2, -1,  2,
       0,  1,  0,  1,  0,  2,  0,  2,
       0,  0,  0,  0,  0,  0,  0,  0,

       1, -1,  1, -1,  0, -1,  0, -1,
       0, -1,  0, -1,  1, -1,  1, -1,
      -2,  0, -2,  0, -2, -1, -2, -1,
      -1,  1, -1,  1, -1, -1, -1, -1,] # for 105 mm

    from scitbx.array_family import flex
    total_tile_translations = flex.int(corrected_auxiliary_translations)  - \
                              flex.int(LB67_corrections_post009)

#                              flex.int(LC06_mark10_001_corrections)
#                              flex.int(LC06_mark10_101_corrections)

    TT = list(total_tile_translations)
    working_extract.distl.tile_translations = TT

    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y. Optimized for 267 mm.
    #working_extract.distl.quad_translations = [-3,  2,
    #                                           -8,  3,
    #                                           -7,  8,
    #                                           -10,  8]

    # determined for LC06_runs15-17. Optimized for 159 mm
    #working_extract.distl.quad_translations = [-2,  1,
    #                                           -9,  3,
    #                                           -5,  8,
    #                                           -9,  10]

    # determined for LC67_run26. Optimized for 101 mm; detz_offset must be changed to move 100 mm distance to 101
    working_extract.distl.quad_translations = [ 0,  3,
                                               -7,  3,
                                               -5,  9,
                                               -9,  9]

    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y. Optimized for 79 mm.
    #working_extract.distl.quad_translations = [-2,  3,
    #                                           -9,  5,
     #                                          -5,  10,
     #                                          -9,  12]
    return working_extract

  elif cxi_version in ["CXI 9.1"]:
    working_extract = working_phil.command_extractor

    from scitbx.array_family import flex
    # Based on the LD91 experiment, runs 95-114, redetermined 3/17/15 NKS
    corrected_auxiliary_translations = flex.int([
       0,  8,  0,  8, -1, -1, -1, -1,
       8,  2,  8,  2,  3,  2,  3,  2,
       5, -6,  5, -6,  3, -1,  3, -1,
       2, -1,  2, -1, -2,  1, -2,  1,

       4,  1,  4,  1,  1,  1,  1,  1,
       0, -7,  0, -7,  2, -2,  2, -2,
      -1, -4, -1, -4,  8, -4,  8, -4,
       7, -7,  7, -7,  7, -1,  7, -1,

       1, -5,  1, -5,  2, -1,  2, -1,
      -8, -3, -8, -3, -1, -5, -1, -5,
      -4,  6, -4,  6, -7, -3, -7, -3,
      -7,  1, -7,  1,  1, -1,  1, -1,

      -8,  3, -8,  3,  0,  1,  0,  1,
      -3,  8, -3,  8, -4,  1, -4,  1,
       0,  6,  0,  6, -7,  7, -7,  7,
      -3,  2, -3,  2, -3,  3, -3,  3])

    working_extract.distl.tile_translations = list(corrected_auxiliary_translations)

    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y.
    # Determined for LD91 run 33. Optimized for 125 mm
    working_extract.distl.quad_translations = [ 8, -3,
                                               -9,  5,
                                                9,  7,
                                               -5, 15]

    return working_extract

  elif cxi_version in ["CXI 8.d"]:
    working_extract = working_phil.command_extractor

    from scitbx.array_family import flex
    total_tile_translations = flex.int(128)

    TT = list(total_tile_translations)
    working_extract.distl.tile_translations = TT



    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y.
    working_extract.distl.quad_translations = [11,  3,
                                               -3, -2,
                                                0,  8,
                                                0,  4]
    print(len(working_extract.distl.tile_translations))

    return working_extract

  elif cxi_version in ["CXI 10.1"]:
    working_extract = working_phil.command_extractor

    from scitbx.array_family import flex

    # Determined from LG36, trial 305
    corrected_auxiliary_translations =flex.int([
       1,  6,  1,  6,  0,  0,  0,  0,
       6,  1,  6,  1,  2,  2,  2,  2,
       5, -6,  5, -6,  3, -3,  3, -3,
       3, -3,  3, -3, -1,  2, -1,  2,

       4,  0,  4,  0,  0,  0,  0,  0,
       2, -7,  2, -7,  2, -1,  2, -1,
      -3, -7, -3, -7,  3, -6,  3, -6,
       1, -9,  1, -9,  0, -2,  0, -2,

       0, -6,  0, -6,  0,  0,  0,  0,
     -10, -1,-10, -1, -1, -3, -1, -3,
      -7,  5, -7,  5,  0,  0,  0,  0,
      -7, -2, -7, -2, -1, -2, -1, -2,

      -7,  0, -7,  0,  0,  0,  0,  0,
      -5,  6, -5,  6, -5, -4, -5, -4,
      -1,  3, -1,  3,  1, -1,  1, -1,
      -3,  5, -3,  5, -2,  1, -2,  1])

    working_extract.distl.tile_translations = list(corrected_auxiliary_translations)

    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y.
    # Determined for LG36 run 82. Optimized for 118 mm
    working_extract.distl.quad_translations = [8,  3,
                                               5,  3,
                                               7, 11,
                                               3, 12]


    return working_extract

  elif cxi_version in ["CXI 10.2"]:
    working_extract = working_phil.command_extractor

    from scitbx.array_family import flex

    corrected_auxiliary_translations =flex.int([
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,

       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,

       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,

       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0])


    working_extract.distl.tile_translations = list(corrected_auxiliary_translations)

    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y.
    # Determined for LG36 run 82. Optimized for 118 mm
    working_extract.distl.quad_translations = [ 0, 0,
                                                0, 0,
                                                0, 0,
                                                0, 0]

    return working_extract

  elif cxi_version in ["CXI 11.1"]:
    working_extract = working_phil.command_extractor

    from scitbx.array_family import flex

    corrected_auxiliary_translations =flex.int([
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,

       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,

       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,

       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0])


    working_extract.distl.tile_translations = list(corrected_auxiliary_translations)

    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y.
    working_extract.distl.quad_translations = [ 0, 0,
                                                0, 0,
                                                0, 0,
                                                0, 0]


    return working_extract

  elif cxi_version in ["CXI 11.2"]:
    working_extract = working_phil.command_extractor

    from scitbx.array_family import flex

    # Determined from LH80, runs 18-23, trial 014_000
    """
    corrected_auxiliary_translations =flex.int([
       2,  1,  2,  1,  1,  2,  1,  2,
       3,  0,  3,  0,  2,  0,  2,  0,
       1,  0,  1,  0,  0,  0,  0,  0,
       1,  1,  1,  1,  0,  1,  0,  1,

      -1,  0, -1,  0,  0,  0,  0,  0,
       0,  0,  0,  0, -1, -1, -1, -1,
       2, -2,  2, -2,  0,  0,  0,  0,
       1, -3,  1, -3,  0,  0,  0,  0,

       0, -3,  0, -3,  0, -2,  0, -2,
      -2, -4, -2, -4, -1, -3, -1, -3,
      -2, -3, -2, -3,  0,  0,  0,  0,
      -3, -1, -3, -1, -1, -2, -1, -2,

      -1,  2, -1,  2, -2,  1, -2,  1,
      -2,  2, -2,  2, -3,  2, -3,  2,
      -2,  2, -2,  2,  0,  0,  0,  0,
      -1, -2, -1, -2, -1,  0, -1,  0])
    """
    # Determined from LH80, runs 46-51, trial 019_000
    corrected_auxiliary_translations =flex.int([
       2,  1,  2,  1,  2,  2,  2,  2,
       2, -2,  2, -2,  2, -1,  2, -1,
      -1, -2, -1, -2, -3, -6, -3, -6,
      -1,  1, -1,  1,  0,  1,  0,  1,

      -2,  1, -2,  1, -1,  0, -1,  0,
     # row from 014_001
     # 0,  0,  0,  0, -1, -1, -1, -1,
     # original row from 019_000
     #-2,  1, -2,  1, -3,  0, -3,  0,
     # the minus 2s and 3s above put this sensor off of the edge.
       0,  1,  0,  1, -1,  0, -1,  0,
       1, -2,  1, -2,  1, -1,  1, -1,
       1, -4,  1, -4,  0,  0,  0,  0,

       1, -3,  1, -3,  0, -2,  0, -2,
      -1, -5, -1, -5,  0, -4,  0, -4,
      -3, -4, -3, -4, -3, -6, -3, -6,
      -1, -1, -1, -1,  0, -2,  0, -2,

      -1,  1, -1,  1, -2,  1, -2,  1,
      -2,  1, -2,  1, -3,  1, -3,  1,
      -1,  0, -1,  0,  0,  0,  0,  0,
       0, -3,  0, -3, -1,  0, -1,  0])

    working_extract.distl.tile_translations = list(corrected_auxiliary_translations)

    # Order: UL x, UL y, UR x, UR y, LL x, LL y, LR x, LR y.
    # Manually determined from LH80 run 12. Optimized for 85 mm
    #working_extract.distl.quad_translations = [-11, 2,
    #                                            -8, 2,
    #                                            -9, 2,
    #                                           -11, 4]
    # Manually determined from LH80 runs 31-37. Optimized for 85 mm. Good for runs 31+
    working_extract.distl.quad_translations = [-16, 5,
                                               -11, 5,
                                               -14, 5,
                                               -16, 6]
    return working_extract


  elif cxi_version in ["Sacla.MPCCD"]:
    working_extract = working_phil.command_extractor
    working_extract.distl.quad_translations = None
    working_extract.distl.tile_translations = None
    return working_extract

  elif cxi_version in ["Sacla.MPCCD.8tile"]:
    working_extract = working_phil.command_extractor
    working_extract.distl.quad_translations = None

    from scitbx.array_family import flex
    corrected_auxiliary_translations =flex.int([
       0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,
    ])

    working_extract.distl.tile_translations = list(corrected_auxiliary_translations)

    #working_extract.distl.tile_translations = None
    return working_extract

  else:
    return working_phil.command_extractor
