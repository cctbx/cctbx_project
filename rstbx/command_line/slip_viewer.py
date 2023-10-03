from __future__ import absolute_import, division, print_function
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cctbx.image_viewer
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
#
# $Id$

import os
import sys
import wx

from rstbx.slip_viewer.frame import XrayFrame
from iotbx import phil
from spotfinder import phil_str
from spotfinder.command_line.signal_strength import additional_spotfinder_phil_defs
from rstbx.phil.phil_preferences import iotbx_defs_viewer_detail

master_str="""
anti_aliasing = False
  .type = bool
  .help = "Enable anti-aliasing"
effective_metrology = None
  .type = path
  .help = "Read effective metrology from PATH"
beam_center = True
  .type = bool
  .help = "Mark beam center"
show_integration_results = False
  .type = bool
  .help = "Show integration results"
show_spotfinder_results = False
  .type = bool
  .help = "Show spotfinder results"
show_effective_tiling = False
  .type = bool
  .help = "Show the effective tiling of the detector"
show_untrusted = False
  .type = bool
  .help = "Always highlyight untrusted pixels, regardless of detector type"
""" + iotbx_defs_viewer_detail

def run(argv=None):
  if (argv is None):
    argv = sys.argv

  # XXX Could/should handle effective metrology the same way, except
  # it does not have a single scope.
  work_phil = phil.process_command_line(
    args=argv[1:],
    master_string=master_str + phil_str + additional_spotfinder_phil_defs)
  work_params = work_phil.work.extract()

  app = wx.App(0)
  wx.SystemOptions.SetOption("osx.openfiledialog.always-show-types", "1")
  frame = XrayFrame(None, -1, "X-ray image display", size=(800,720))
  frame.Show()

  # show settings panel
  frame.OnShowSettings(None)
  frame.settings_frame.panel.center_ctrl.SetValue(
    work_params.beam_center)
  frame.settings_frame.panel.integ_ctrl.SetValue(
    work_params.show_integration_results)
  frame.settings_frame.panel.spots_ctrl.SetValue(
    work_params.show_spotfinder_results)
  frame.settings.show_effective_tiling = work_params.show_effective_tiling
  frame.settings_frame.panel.collect_values()

  if (work_params.effective_metrology is not None):
    from serialtbx.detector.legacy_metrology.metrology import \
      master_phil, metrology_as_transformation_matrices

    stream = open(work_params.effective_metrology)
    metrology_phil = master_phil.fetch(sources=[phil.parse(stream.read())])
    stream.close()
    frame.metrology_matrices = metrology_as_transformation_matrices(
      metrology_phil.extract())

  # Update initial settings with values from the command line.  Needs
  # to be done before image is loaded (but after the frame is
  # instantiated).
  frame.params = work_params
  frame.init_pyslip()
  frame.pyslip.tiles.user_requests_antialiasing = work_params.anti_aliasing
  frame.pyslip.tiles.show_untrusted = frame.params.show_untrusted

  paths = work_phil.remaining_args
  if (len(paths) == 1 and os.path.basename(paths[0]) == "DISTL_pickle"):
    assert os.path.isfile(paths[0])
    frame.load_distl_output(paths[0])
  elif (len(paths) > 0):
    frame.CHOOSER_SIZE = 1500

    # from dxtbx.imageset import ImageSetFactory
    # sets = ImageSetFactory.new(paths)
    from dxtbx.model.experiment_list import ExperimentListFactory
    from rstbx.slip_viewer.frame import chooser_wrapper
    experiments = ExperimentListFactory.from_filenames(paths)
    sets = experiments.imagesets()

    for imgset in sets:
      for idx in imgset.indices():
        frame.add_file_name_or_data(chooser_wrapper(imgset, idx))
    idx = sets[0].indices()[0]
    frame.load_image(chooser_wrapper(sets[0],idx))

  app.MainLoop()

  return 0


if (__name__ == "__main__"):
  sys.exit(run())
