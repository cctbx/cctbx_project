from __future__ import absolute_import, division, print_function
from six.moves import range
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.frame_extractor
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.frame_extractor
#
# $Id: frame_extractor.py idyoung $

from dials.util.options import Importer, flatten_reflections, flatten_experiments, ArgumentParser
import iotbx.phil
import cctbx, os, glob
from libtbx import easy_pickle
from six.moves import zip
from serialtbx.util.construct_frame import ConstructFrame

phil_scope = iotbx.phil.parse("""
  input {
    experiments = None
      .type = path
      .help = path to an experiments.expt file
    reflections = None
      .type = path
      .help = path to a reflection table (integrated.refl) file
  }
  output {
    filename = None
      .type = str
      .help = if set, name of final pickle file
    dirname = None
      .type = path
      .help = if set, path to directory to save the new pickle file
    }
    """)

class ConstructFrameFromFiles(ConstructFrame):
  def __init__(self, refl_name, json_name, outname=None):
    # load the integration.refl file (reflection table) into memory and
    # load the experiments.expt file (json) into memory, piecewise.
    # check_format=False because we don't want to load any imagesets in the
    # experiment list
    importer = Importer([refl_name, json_name], read_experiments=True, read_reflections=True, check_format=False)
    if importer.unhandled:
      print("unable to process:", importer.unhandled)
    reflections_l = flatten_reflections(importer.reflections)
    experiments_l = flatten_experiments(importer.experiments)
    assert len(experiments_l) == 1, "cannot construct a single frame from multiple experiments"
    frame = ConstructFrame.__init__(self, reflections_l[0], experiments_l[0])
    if frame is not None:
      self.frame.make_frame()

def construct_frames_from_files(refl_name, json_name, outname=None, outdir=None):
  importer = Importer([refl_name, json_name], read_experiments=True, read_reflections=True, check_format=False)
  if importer.unhandled:
    print("unable to process:", importer.unhandled)
  reflections_l = flatten_reflections(importer.reflections)[0]
  experiments_l = flatten_experiments(importer.experiments)
  frames = []
  if outdir is None:
    outdir = '.'
  if outname is None:
    outname = 'int-%d' + refl_name.split('.pickle')[0] + '_extracted.pickle'
  elif '%' not in outname:
    outname = outname.split(".pickle")[0] + ("_%d.pickle")
  for i in range(len(experiments_l)):
    refl = reflections_l.select(reflections_l['id'] == i)
    if len(refl) == 0: continue
    expt = experiments_l[i]
    frame = ConstructFrame(refl, expt).make_frame()
    name = outname % i
    easy_pickle.dump(os.path.join(outdir, name), frame)

if __name__ == "__main__":
  parser = ArgumentParser(phil=phil_scope)
  params, options = parser.parse_args(show_diff_phil=True)
  if params.output.dirname is not None:
    assert os.path.isdir(params.output.dirname)
  for refl_name, json_name in zip(sorted(glob.glob(params.input.reflections)), sorted(glob.glob(params.input.experiments))):
    if params.output.filename is None:
      basename = os.path.basename(refl_name)
      name = os.path.splitext(basename)[0] + "_extracted.pickle"
    else:
      name = params.output.filename
    construct_frames_from_files(refl_name, json_name, outname=name, outdir=params.output.dirname)
