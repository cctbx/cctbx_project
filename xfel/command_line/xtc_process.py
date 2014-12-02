from __future__ import division
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.xtc_process
#
from psana import *
import numpy as np
from xfel.cftbx.detector import cspad_cbf_tbx
from xfel.cxi.cspad_ana import cspad_tbx
import pycbf, os, sys
import libtbx.load_env
from libtbx.utils import Sorry, Usage
from dials.util.options import OptionParser
from libtbx.phil import parse
from dxtbx.imageset import MemImageSet
from dxtbx.datablock import DataBlockFactory

phil_scope = parse('''
  dispatch {
    max_events = None
      .type = int
      .help = "If not specified, process all events. Otherwise, only process this many"
    hit_finder = True
      .type = bool
      .help = "If the number of strong reflections is less than"
      .help = "refinement.reflections.minimum_number_of_reflections, hit_filter=True"
      .help = "will discard the image. hit_finder=False: process all images"
    index = True
      .type = bool
      .help = "Attempt to index images"
    integrate = True
      .type = bool
      .help = "Integrated indexed images. Ignored if index=False"
    dump_strong = False
      .type = bool
      .help = "Save strongly diffracting images to cbf format"
    dump_indexed = True
      .type = bool
      .help = "Save indexed images to cbf format"
    dump_all = False
      .type = bool
      .help = "All frames will be saved to cbf format if set to True"
  }
  input {
    cfg = None
      .type = str
      .help = "Path to psana config file"
    experiment = None
      .type = str
      .help = "Experiment identifier, e.g. cxi84914"
    run_num = None
      .type = int
      .help = "Run number or run range to process"
    address = None
      .type = str
      .help = "Detector address, e.g. CxiDs2.0:Cspad.0"
    detz_offset = None
      .type = int
      .help = "Distance from back of detector rail to sample interaction region (CXI)"
      .help = "or actual detector distance (XPP)"
  }
  output {
    output_dir = .
      .type = str
      .help = "Directory output files will be placed"
    datablock_filename = %s_datablock.json
      .type = str
      .help = "The filename for output datablock"
    strong_filename = %s_strong.pickle
      .type = str
      .help = "The filename for strong reflections from spot finder output."
    indexed_filename = %s_indexed.pickle
      .type = str
      .help = "The filename for indexed reflections."
    refined_experiments_filename = %s_refined_experiments.json
      .type = str
      .help = "The filename for saving refined experimental models"
    integrated_filename = %s_integrated.pickle
      .type = str
      .help = "The filename for final integrated reflections."
    profile_filename = None
      .type = str
      .help = "The filename for output reflection profile parameters"
  }

  border_mask {
    include scope dials.command_line.generate_mask.phil_scope
  }
  include scope dials.algorithms.peak_finding.spotfinder_factory.phil_scope
  include scope dials.algorithms.indexing.indexer.master_phil_scope
  include scope dials.algorithms.refinement.refiner.phil_scope
  include scope dials.algorithms.profile_model.factory.phil_scope
  include scope dials.algorithms.integration.integrator.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope
''', process_includes=True)

# work around for copying dxtbx FormatCBFCspad objects
from xfel.cftbx.detector.cspad_cbf_tbx import cbf_wrapper
def __stupid_but_swig_safe__deepcopy__(self, memo):
  pass
cbf_wrapper.__deepcopy__ = __stupid_but_swig_safe__deepcopy__

from dials.command_line.process import Script as DialsProcessScript
class InMemScript(DialsProcessScript):
  """ Script to process XFEL data at LCLS """
  def __init__(self):
    """ Set up the option parser. Arguments come from the command line or a phil file """
    self.usage = \
    """ %s input.cfg=filename.cfg input.experiment=experimentname
    input.run_num=N input.address=address input.detz_offset=N
    """%libtbx.env.dispatcher_name

    self.parser = OptionParser(
      usage = self.usage,
      phil = phil_scope)

  def run(self):
    """ Process all images assigned to this thread """

    params, options = self.parser.parse_args(
      show_diff_phil=True)

    if params.input.cfg is None or \
        params.input.experiment is None or \
        params.input.run_num is None or \
        params.input.address is None or \
        params.input.detz_offset is None:
      raise Usage(self.usage)

    if not os.path.exists(params.output.output_dir):
      raise Sorry("Output path not found:" + params.output.output_dir)

    # The convention is to put %s in the phil parameter to add a time stamp to
    # each output datafile. Save the initial templates here.
    strong_filename_template              = params.output.strong_filename
    indexed_filename_template             = params.output.indexed_filename
    refined_experiments_filename_template = params.output.refined_experiments_filename
    integrated_filename_template          = params.output.integrated_filename

    # Don't allow the strong reflections to be written unless there are enough to
    # process
    params.output.strong_filename = None

    # Save the paramters
    self.params = params
    self.options = options

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank() # each process in MPI has a unique id, 0-indexed
    size = comm.Get_size() # size: number of processes running in this job

    # set up psana
    setConfigFile(params.input.cfg)
    dataset_name = "exp=%s:run=%s:idx"%(params.input.experiment,params.input.run_num)
    ds = DataSource(dataset_name)
    src = Source('DetInfo(%s)'%params.input.address)

    env = ds.env()
    calib_dir = env.calibDir()

    # set this to sys.maxint to analyze all events
    if params.dispatch.max_events is None:
      max_events = sys.maxint
    else:
      max_events = params.dispatch.max_events

    for run in ds.runs():
      # load a header only cspad cbf from the slac metrology
      base_dxtbx = cspad_cbf_tbx.env_dxtbx_from_slac_metrology(run.env(), src)
      if base_dxtbx is None:
        raise Sorry("Couldn't load calibration file for run %d"%run.run())

      # list of all events
      times = run.times()
      nevents = min(len(times),max_events)
      mylength = nevents//size # easy but sloppy. lose few events at end of run.
      # chop the list into pieces, depending on rank
      mytimes= times[rank*mylength:(rank+1)*mylength]
      for i in range(mylength):
        evt = run.event(mytimes[i])
        id = evt.get(EventId)
        print "Event #",i," has id:",id
        # the data needs to have already been processed and put into the event by psana
        data = evt.get(ndarray_float64_3, src, 'image0')
        if data is None:
          print "No data"
          continue
        data = data.astype(np.int32)

        distance = cspad_tbx.env_distance(params.input.address, run.env(), params.input.detz_offset)
        if distance is None:
          print "No distance, skipping shot"
          continue

        wavelength = cspad_tbx.evt_wavelength(evt)
        if wavelength is None:
          print "No wavelength, skipping shot"
          continue

        timestamp = cspad_tbx.evt_timestamp(cspad_tbx.evt_time(evt)) # human readable format
        if timestamp is None:
          print "No timestamp, skipping shot"
          continue

        t = timestamp
        s = t[0:4] + t[5:7] + t[8:10] + t[11:13] + t[14:16] + t[17:19] + t[20:23]
        print "Processing shot", s

        # stitch together the header, data and metadata into the final dxtbx format object
        cspad_img = cspad_cbf_tbx.format_object_from_data(base_dxtbx, data, distance, wavelength, timestamp, params.input.address)

        if params.dispatch.dump_all:
          # save cbf file
          dest_path = os.path.join(params.output.output_dir, "shot-" + s + ".cbf")
          cspad_img._cbf_handle.write_widefile(dest_path, pycbf.CBF,\
            pycbf.MIME_HEADERS|pycbf.MSG_DIGEST|pycbf.PAD_4K, 0)

        imgset = MemImageSet([cspad_img])
        datablock = DataBlockFactory.from_imageset(imgset)[0]

        # before calling DIALS for processing, set output paths according to the templates
        if "%s" in indexed_filename_template:
          self.params.output.indexed_filename = os.path.join(params.output.output_dir, indexed_filename_template%("idx-" + s))
        if "%s" in refined_experiments_filename_template:
          self.params.output.refined_experiments_filename = os.path.join(params.output.output_dir, refined_experiments_filename_template%("idx-" + s))
        if "%s" in integrated_filename_template:
          self.params.output.integrated_filename = os.path.join(params.output.output_dir, integrated_filename_template%("idx-" + s))

        # if border is requested, generate a border only mask
        if params.border_mask.border > 0:
          from dials.command_line.generate_mask import MaskGenerator
          generator = MaskGenerator(params.border_mask)
          mask = generator.generate(imgset)

          params.spotfinder.lookup.mask = mask

        try:
          observed = self.find_spots(datablock)
        except Exception, e:
          print str(e)
          continue

        if params.dispatch.hit_finder and len(observed) < params.refinement.reflections.minimum_number_of_reflections:
          print "Not enough spots to index"
          continue

        # save cbf file
        if params.dispatch.dump_strong:
          dest_path = os.path.join(params.output.output_dir, "hit-" + s + ".cbf")
          try:
            cspad_img._cbf_handle.write_widefile(dest_path, pycbf.CBF,\
              pycbf.MIME_HEADERS|pycbf.MSG_DIGEST|pycbf.PAD_4K, 0)
          except Exception:
            print "Warning, couldn't save cbf:", dest_path

          # save strong reflections.  self.find_spots() would have done this, but we only
          # want to save data if it is enough to try and index it
          if strong_filename_template:
            if "%s" in strong_filename_template:
              strong_filename = strong_filename_template%("hit-" + s)
            else:
              strong_filename = strong_filename_template
            strong_filename = os.path.join(params.output.output_dir, strong_filename)

            from dials.util.command_line import Command
            Command.start('Saving {0} reflections to {1}'.format(
                len(observed), os.path.basename(strong_filename)))
            observed.as_pickle(strong_filename)
            Command.end('Saved {0} observed to {1}'.format(
                len(observed), os.path.basename(strong_filename)))

        if not params.dispatch.index:
          continue

        # index and refine
        try:
          experiments, indexed = self.index(datablock, observed)
        except Exception, e:
          print str(e)
          continue

        if params.dispatch.dump_indexed:
          dest_path = os.path.join(params.output.output_dir, "idx-" + s + ".cbf")
          try:
            cspad_img._cbf_handle.write_widefile(dest_path, pycbf.CBF,\
              pycbf.MIME_HEADERS|pycbf.MSG_DIGEST|pycbf.PAD_4K, 0)
          except Exception:
            print "Warning, couldn't save cbf:", dest_path

        try:
          experiments = self.refine(experiments, indexed)
        except Exception, e:
          print str(e)
          continue

        if not params.dispatch.integrate:
          continue

        # integrate
        try:
          integrated = self.integrate(experiments, indexed)
        except Exception, e:
          print str(e)

if __name__ == "__main__":
  from dials.util import halraiser
  try:
    script = InMemScript()
    script.run()
  except Exception as e:
    halraiser(e)
