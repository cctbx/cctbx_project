from __future__ import division
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.xtc_process
#
try:
  import psana
except ImportError:
  pass # for running at home without psdm build
from xfel.cftbx.detector import cspad_cbf_tbx
from xfel.cxi.cspad_ana import cspad_tbx
import pycbf, os, sys, copy, socket
import libtbx.load_env
from libtbx.utils import Sorry, Usage
from dials.util.options import OptionParser
from libtbx.phil import parse
from dxtbx.imageset import MemImageSet
from dxtbx.datablock import DataBlockFactory

xtc_phil_str = '''
  dispatch {
    max_events = None
      .type = int
      .help = If not specified, process all events. Otherwise, only process this many
    hit_finder{
      enable = True
        .type = bool
        .help = Whether to do hitfinding. hit_finder=False: process all images
      minimum_number_of_reflections = 16
        .type = int
        .help = If the number of strong reflections on an image is less than this, and \
                 the hitfinder is enabled, discard this image.
    }
    index = True
      .type = bool
      .help = Attempt to index images
    integrate = True
      .type = bool
      .help = Integrated indexed images. Ignored if index=False
    dump_strong = False
      .type = bool
      .help = Save strongly diffracting images to cbf format
    dump_indexed = True
      .type = bool
      .help = Save indexed images to cbf format
    dump_all = False
      .type = bool
      .help = All frames will be saved to cbf format if set to True
  }
  debug
    .help = Use these flags to track down problematic events that cause unhandled exceptions. \
            Here, a bad event means it caused an unhandled exception, not that the image \
            failed to index. \
            Examples: \
            Process only unprocessed events (excluding bad events): \
              skip_processed_events=True, skip_unprocessed_events=False skip_bad_events=True \
            Process only bad events (for debugging): \
              skip_processed_events=True, skip_unprocessed_events=True skip_bad_events=False \
            Note, due to how MPI works, if an unhandled exception occurrs, some bad events \
            will be marked as bad that were simply in process when the program terminated \
            due to a bad event. Try processing only bad events single process to find the \
            culprit and alert the program authors.
  {
    skip_processed_events = False
      .type = bool
      .help = If True, will look for diagnostic files in the output directory and use \
              them to skip events that had already been processed (succesfully or not)
    skip_unprocessed_events = False
      .type = bool
      .help = If True, will look for diagnostic files in the output directory and use \
              them to skip events that had haven't been processed
    skip_bad_events = False
      .type = bool
      .help = If True, will look for diagnostic files in the output directory and use \
              them to skip events that had caused unhandled exceptions previously
    event_timestamp = None
      .type = str
      .multiple = True
      .help = List of timestamps. If set, will only process the events that match them
  }
  input {
    cfg = None
      .type = str
      .help = Path to psana config file
    experiment = None
      .type = str
      .help = Experiment identifier, e.g. cxi84914
    run_num = None
      .type = int
      .help = Run number or run range to process
    address = None
      .type = str
      .help = Detector address, e.g. CxiDs2.0:Cspad.0, or detector alias, e.g. Ds1CsPad
    override_trusted_max = None
      .type = int
      .help = During spot finding, override the saturation value for this data. \
              Overloads will not be integrated, but they can assist with indexing.
    override_trusted_min = None
      .type = int
      .help = During spot finding and indexing, override the minimum pixel value \
              for this data. This does not affect integration.
    use_ffb = False
      .type = bool
      .help = Run on the ffb if possible. Only for active users!
    xtc_dir = None
      .type = str
      .help = Optional path to data directory if it's non-standard. Only needed if xtc \
              streams are not in the standard location for your PSDM installation.
  }
  format {
    file_format = *cbf pickle
      .type = choice
      .help = Output file format, 64 tile segmented CBF or image pickle
    pickle {
      cfg = None
        .type = str
        .help = Path to psana config file with a mod_image_dict module
      out_key = cctbx.xfel.image_dict
        .type = str
        .help = Key name that mod_image_dict uses to put image data in each psana event
    }
    cbf {
      detz_offset = None
        .type = int
        .help = Distance from back of detector rail to sample interaction region (CXI) \
                or actual detector distance (XPP)
      override_energy = None
        .type = float
        .help = If not None, use the input energy for every event instead of the energy \
                from the XTC stream
      gain_mask_value = None
        .type = float
        .help = If not None, use the gain mask for the run to multiply the low-gain pixels by this number
    }
  }
  output {
    output_dir = .
      .type = str
      .help = Directory output files will be placed
    logging_dir = None
      .type = str
      .help = Directory output log files will be placed
    datablock_filename = %s_datablock.json
      .type = str
      .help = The filename for output datablock
    strong_filename = %s_strong.pickle
      .type = str
      .help = The filename for strong reflections from spot finder output.
    indexed_filename = %s_indexed.pickle
      .type = str
      .help = The filename for indexed reflections.
    refined_experiments_filename = %s_refined_experiments.json
      .type = str
      .help = The filename for saving refined experimental models
    integrated_filename = %s_integrated.pickle
      .type = str
      .help = The filename for final integrated reflections.
    profile_filename = None
      .type = str
      .help = The filename for output reflection profile parameters
  }
  mp {
    method = *mpi sge
      .type = choice
      .help = Muliprocessing method
  }
'''

dials_phil_str = '''
  verbosity = 1
   .type = int(value_min=0)
   .help = The verbosity level
  border_mask {
    include scope dials.command_line.generate_mask.phil_scope
  }
  include scope dials.algorithms.peak_finding.spotfinder_factory.phil_scope
  include scope dials.algorithms.indexing.indexer.master_phil_scope
  include scope dials.algorithms.refinement.refiner.phil_scope
  include scope dials.algorithms.profile_model.factory.phil_scope
  include scope dials.algorithms.integration.integrator.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope
'''

phil_scope = parse(xtc_phil_str + dials_phil_str, process_includes=True)

# work around for copying dxtbx FormatCBFCspad objects
from xfel.cftbx.detector.cspad_cbf_tbx import cbf_wrapper
def __stupid_but_swig_safe__deepcopy__(self, memo):
  pass
cbf_wrapper.__deepcopy__ = __stupid_but_swig_safe__deepcopy__

from xfel.command_line.xfel_process import Script as DialsProcessScript
class InMemScript(DialsProcessScript):
  """ Script to process XFEL data at LCLS """
  def __init__(self):
    """ Set up the option parser. Arguments come from the command line or a phil file """
    self.usage = """
%s input.cfg=filename.cfg input.experiment=experimentname input.run_num=N input.address=address
 format.file_format=cbf format.cbf.detz_offset=N
%s input.cfg=filename.cfg input.experiment=experimentname input.run_num=N input.address=address
format.file_format=pickle format.pickle.cfg=path
    """%(libtbx.env.dispatcher_name, libtbx.env.dispatcher_name)
    self.parser = OptionParser(
      usage = self.usage,
      phil = phil_scope)

  def run(self):
    """ Process all images assigned to this thread """

    params, options = self.parser.parse_args(
      show_diff_phil=True)

    # Configure the logging
    from dials.util import log
    log.config(params.verbosity)

    # Check inputs
    if params.input.experiment is None or \
       params.input.run_num is None or \
       params.input.address is None:
      raise Usage(self.usage)

    if params.format.file_format == "cbf":
      if params.format.cbf.detz_offset is None:
        raise Usage(self.usage)
    elif params.format.file_format == "pickle":
      if params.format.pickle.cfg is None:
        raise Usage(self.usage)
    else:
      raise Usage(self.usage)

    if not os.path.exists(params.output.output_dir):
      raise Sorry("Output path not found:" + params.output.output_dir)

    # The convention is to put %s in the phil parameter to add a time stamp to
    # each output datafile. Save the initial templates here.
    self.strong_filename_template              = params.output.strong_filename
    self.indexed_filename_template             = params.output.indexed_filename
    self.refined_experiments_filename_template = params.output.refined_experiments_filename
    self.integrated_filename_template          = params.output.integrated_filename

    # Don't allow the strong reflections to be written unless there are enough to
    # process
    params.output.strong_filename = None

    # Save the paramters
    self.params_cache = copy.deepcopy(params)
    self.options = options

    if params.mp.method == "mpi":
      from mpi4py import MPI
      comm = MPI.COMM_WORLD
      rank = comm.Get_rank() # each process in MPI has a unique id, 0-indexed
      size = comm.Get_size() # size: number of processes running in this job
    elif params.mp.method == "sge" and \
        'SGE_TASK_ID'    in os.environ and \
        'SGE_TASK_FIRST' in os.environ and \
        'SGE_TASK_LAST'  in os.environ:
      if 'SGE_STEP_SIZE' in os.environ:
        assert int(os.environ['SGE_STEP_SIZE']) == 1
      if os.environ['SGE_TASK_ID'] == 'undefined' or os.environ['SGE_TASK_ID'] == 'undefined' or os.environ['SGE_TASK_ID'] == 'undefined':
        rank = 0
        size = 1
      else:
        rank = int(os.environ['SGE_TASK_ID']) - int(os.environ['SGE_TASK_FIRST'])
        size = int(os.environ['SGE_TASK_LAST']) - int(os.environ['SGE_TASK_FIRST']) + 1
    else:
      rank = 0
      size = 1

    if params.output.logging_dir is not None:
      log_path = os.path.join(params.output.logging_dir, "log_rank%04d.out"%rank)
      error_path = os.path.join(params.output.logging_dir, "error_rank%04d.out"%rank)
      print "Redirecting stdout to %s"%log_path
      print "Redirecting stderr to %s"%error_path
      assert os.path.exists(log_path)
      sys.stdout = open(log_path,'a', buffering=0)
      sys.stderr = open(error_path,'a',buffering=0)
      print "Should be redirected now"

    debug_dir = os.path.join(params.output.output_dir, "debug")
    if not os.path.exists(debug_dir):
      os.makedirs(debug_dir)

    if params.debug.skip_processed_events or params.debug.skip_processed_events or params.debug.skip_bad_events:
      print "Reading debug files..."
      self.known_events = {}
      for filename in os.listdir(debug_dir):
        # format: hostname,timestamp,status
        for line in open(os.path.join(debug_dir, filename)):
          vals = line.strip().split(',')
          if len(vals) == 2:
            self.known_events[vals[1]] = "unknown"
          elif len(vals) == 3:
            self.known_events[vals[1]] = vals[2]

    debug_file_path = os.path.join(debug_dir, "debug_%d.txt"%rank)
    write_newline = os.path.exists(debug_file_path)
    self.debug_file_handle = open(debug_file_path, 'a', 0) # 0 for unbuffered
    if write_newline: # needed if the there was a crash
      self.debug_file_handle.write("\n")

    # set up psana
    if params.format.file_format=="pickle":
      psana.setConfigFile(params.format.pickle.cfg)
    dataset_name = "exp=%s:run=%s:idx"%(params.input.experiment,params.input.run_num)
    if params.input.xtc_dir is not None:
      if params.input.use_ffb:
        raise Sorry("Cannot specify the xtc_dir and use SLAC's ffb system")
      dataset_name += ":dir=%s"%params.input.xtc_dir
    elif params.input.use_ffb:
      # as ffb is only at SLAC, ok to hardcode /reg/d here
      dataset_name += ":dir=/reg/d/ffb/%s/%s/xtc"%(params.input.experiment[0:3],params.input.experiment)
    ds = psana.DataSource(dataset_name)

    if params.format.file_format == "cbf":
      self.src = psana.Source('DetInfo(%s)'%params.input.address)
      self.psana_det = psana.Detector(params.input.address, ds.env())

    # set this to sys.maxint to analyze all events
    if params.dispatch.max_events is None:
      max_events = sys.maxint
    else:
      max_events = params.dispatch.max_events

    for run in ds.runs():
      if params.format.file_format == "cbf":
        # load a header only cspad cbf from the slac metrology
        self.base_dxtbx = cspad_cbf_tbx.env_dxtbx_from_slac_metrology(run, params.input.address)
        if self.base_dxtbx is None:
          raise Sorry("Couldn't load calibration file for run %d"%run.run())

        if params.format.cbf.gain_mask_value is not None:
          self.gain_mask = self.psana_det.gain_mask(gain=params.format.cbf.gain_mask_value)

      # list of all events
      times = run.times()
      nevents = min(len(times),max_events)
      if params.mp.method == "mpi" and size > 2:
        # use a client/server approach to be sure every process is busy as much as possible
        # only do this if there are more than 2 processes, as one process will be a server
        if rank == 0:
          # server process
          for t in times[:nevents]:
            # a client process will indicate it's ready by sending its rank
            rankreq = comm.recv(source=MPI.ANY_SOURCE)
            comm.send(t,dest=rankreq)
          # send a stop command to each process
          for rankreq in range(size-1):
            rankreq = comm.recv(source=MPI.ANY_SOURCE)
            comm.send('endrun',dest=rankreq)
        else:
          # client process
          while True:
            # inform the server this process is ready for an event
            comm.send(rank,dest=0)
            evttime = comm.recv(source=0)
            if evttime == 'endrun': break
            self.process_event(run, evttime)
      else:
        # chop the list into pieces, depending on rank.  This assigns each process
        # events such that the get every Nth event where N is the number of processes
        mytimes = [times[i] for i in xrange(nevents) if (i+rank)%size == 0]

        for i in xrange(len(mytimes)):
          self.process_event(run, mytimes[i])

      run.end()
    ds.end()

  def process_event(self, run, timestamp):
    """
    Process a single event from a run
    @param run psana run object
    @param timestamp psana timestamp object
    """
    ts = cspad_tbx.evt_timestamp((timestamp.seconds(),timestamp.nanoseconds()/1e6))
    if ts is None:
      print "No timestamp, skipping shot"
      return

    if len(self.params_cache.debug.event_timestamp) > 0 and ts not in self.params_cache.debug.event_timestamp:
      return

    if self.params_cache.debug.skip_processed_events or self.params_cache.debug.skip_unprocessed_events or self.params_cache.debug.skip_bad_events:
      if ts in self.known_events:
        if self.known_events[ts] == "unknown":
          if self.params_cache.debug.skip_bad_events and self.known_events[ts] == "unknown":
            print "Skipping event %s: possibly caused an unknown exception previously"%ts
            return
        elif self.params_cache.debug.skip_processed_events:
          print "Skipping event %s: processed successfully previously"%ts
          return
      else:
        if self.params_cache.debug.skip_unprocessed_events:
          print "Skipping event %s: not processed previously"%ts
          return

    print "Accepted", ts

    self.debug_file_handle.write("%s,%s"%(socket.gethostname(), ts))

    self.params = copy.deepcopy(self.params_cache)

    evt = run.event(timestamp)
    id = evt.get(psana.EventId)
    if evt.get("skip_event"):
      print "Skipping event",id
      self.debug_file_handle.write(",psana_skip\n")
      return

    # the data needs to have already been processed and put into the event by psana
    if self.params.format.file_format == 'cbf':
      # get numpy array, 32x185x388
      data = self.psana_det.calib(evt) # applies psana's complex run-dependent calibrations
      if data is None:
        print "No data"
        self.debug_file_handle.write(",no_data\n")
        return

      if self.params.format.cbf.gain_mask_value is not None:
        # apply gain mask
        data *= self.gain_mask

      distance = cspad_tbx.env_distance(self.params.input.address, run.env(), self.params.format.cbf.detz_offset)
      if distance is None:
        print "No distance, skipping shot"
        self.debug_file_handle.write(",no_distance\n")
        return

      if self.params.format.cbf.override_energy is None:
        wavelength = cspad_tbx.evt_wavelength(evt)
        if wavelength is None:
          print "No wavelength, skipping shot"
          self.debug_file_handle.write(",no_wavelength\n")
          return
      else:
        wavelength = 12398.4187/self.params.format.cbf.override_energy

    if self.params.format.file_format == 'pickle':
      image_dict = evt.get(self.params.format.pickle.out_key)
      data = image_dict['DATA']

    timestamp = t = ts
    s = t[0:4] + t[5:7] + t[8:10] + t[11:13] + t[14:16] + t[17:19] + t[20:23]
    print "Processing shot", s

    if self.params.format.file_format == 'cbf':
      # stitch together the header, data and metadata into the final dxtbx format object
      cspad_img = cspad_cbf_tbx.format_object_from_data(self.base_dxtbx, data, distance, wavelength, timestamp, self.params.input.address)
    elif self.params.format.file_format == 'pickle':
      from dxtbx.format.FormatPYunspecifiedStill import FormatPYunspecifiedStillInMemory
      cspad_img = FormatPYunspecifiedStillInMemory(image_dict)

    cspad_img.timestamp = s

    if self.params.dispatch.dump_all:
      self.save_image(cspad_img, self.params, os.path.join(self.params.output.output_dir, "shot-" + s))

    self.cache_ranges(cspad_img, self.params)

    imgset = MemImageSet([cspad_img])
    datablock = DataBlockFactory.from_imageset(imgset)[0]

    # before calling DIALS for processing, set output paths according to the templates
    if self.indexed_filename_template is not None and "%s" in self.indexed_filename_template:
      self.params.output.indexed_filename = os.path.join(self.params.output.output_dir, self.indexed_filename_template%("idx-" + s))
    if "%s" in self.refined_experiments_filename_template:
      self.params.output.refined_experiments_filename = os.path.join(self.params.output.output_dir, self.refined_experiments_filename_template%("idx-" + s))
    if "%s" in self.integrated_filename_template:
      self.params.output.integrated_filename = os.path.join(self.params.output.output_dir, self.integrated_filename_template%("idx-" + s))

    # if border is requested, generate a border only mask
    if self.params.border_mask.border > 0:
      from dials.command_line.generate_mask import MaskGenerator
      generator = MaskGenerator(self.params.border_mask)
      mask = generator.generate(imgset)

      self.params.spotfinder.lookup.mask = mask

    try:
      observed = self.find_spots(datablock)
    except Exception, e:
      import traceback; traceback.print_exc()
      print str(e), "event", timestamp
      self.debug_file_handle.write(",spotfinding_exception\n")
      return

    print "Found %d bright spots"%len(observed)

    if self.params.dispatch.hit_finder.enable and len(observed) < self.params.dispatch.hit_finder.minimum_number_of_reflections:
      print "Not enough spots to index"
      self.debug_file_handle.write(",not_enough_spots_%d\n"%len(observed))
      return

    self.restore_ranges(cspad_img, self.params)

    # save cbf file
    if self.params.dispatch.dump_strong:
      self.save_image(cspad_img, self.params, os.path.join(self.params.output.output_dir, "hit-" + s))

      # save strong reflections.  self.find_spots() would have done this, but we only
      # want to save data if it is enough to try and index it
      if self.strong_filename_template:
        if "%s" in self.strong_filename_template:
          strong_filename = self.strong_filename_template%("hit-" + s)
        else:
          strong_filename = self.strong_filename_template
        strong_filename = os.path.join(self.params.output.output_dir, strong_filename)

        from dials.util.command_line import Command
        Command.start('Saving {0} reflections to {1}'.format(
            len(observed), os.path.basename(strong_filename)))
        observed.as_pickle(strong_filename)
        Command.end('Saved {0} observed to {1}'.format(
            len(observed), os.path.basename(strong_filename)))

    if not self.params.dispatch.index:
      self.debug_file_handle.write(",strong_shot_%d\n"%len(observed))
      return

    # index and refine
    try:
      experiments, indexed = self.index(datablock, observed)
    except Exception, e:
      import traceback; traceback.print_exc()
      print str(e), "event", timestamp
      self.debug_file_handle.write(",indexing_failed_%d\n"%len(observed))
      return

    if self.params.dispatch.dump_indexed:
      self.save_image(cspad_img, self.params, os.path.join(self.params.output.output_dir, "idx-" + s))

    try:
      experiments = self.refine(experiments, indexed)
    except Exception, e:
      import traceback; traceback.print_exc()
      print str(e), "event", timestamp
      self.debug_file_handle.write(",refine_failed_%d\n"%len(indexed))
      return

    if not self.params.dispatch.integrate:
      self.debug_file_handle.write(",index_ok_%d\n"%len(indexed))
      return

    # integrate
    try:
      integrated = self.integrate(experiments, indexed)
    except Exception, e:
      import traceback; traceback.print_exc()
      print str(e), "event", timestamp
      self.debug_file_handle.write(",integrate_failed_%d\n"%len(indexed))

    self.debug_file_handle.write(",integrate_ok_%d\n"%len(integrated))

  def save_image(self, image, params, root_path):
    """ Save an image, in either cbf or pickle format.
    @param image dxtbx format object
    @param params phil scope object
    @param root_path output file path without extension
    """

    if params.format.file_format == 'cbf':
      dest_path = root_path + ".cbf"
    elif params.format.file_format == 'pickle':
      dest_path = root_path + ".pickle"

    try:
      if params.format.file_format == 'cbf':
        image._cbf_handle.write_widefile(dest_path, pycbf.CBF,\
          pycbf.MIME_HEADERS|pycbf.MSG_DIGEST|pycbf.PAD_4K, 0)
      elif params.format.file_format == 'pickle':
        from libtbx import easy_pickle
        easy_pickle.dump(dest_path, image._image_file)
    except Exception:
      print "Warning, couldn't save image:", dest_path

  def cache_ranges(self, cspad_img, params):
    """ Save the current trusted ranges, and replace them with the given overrides, if present.
    @param cspad_image dxtbx format object
    @param params phil scope
    """
    if params.input.override_trusted_max is None and params.input.override_trusted_min is None:
      return

    detector = cspad_img.get_detector()
    self.cached_ranges = []
    for panel in detector:
      new_range = cached_range = panel.get_trusted_range()
      self.cached_ranges.append(cached_range)
      if params.input.override_trusted_max is not None:
        new_range = new_range[0], params.input.override_trusted_max
      if params.input.override_trusted_min is not None:
        new_range = params.input.override_trusted_min, new_range[1]

      panel.set_trusted_range(new_range)

  def restore_ranges(self, cspad_img, params):
    """ Restore the previously cached trusted ranges, if present.
    @param cspad_image dxtbx format object
    @param params phil scope
    """
    if params.input.override_trusted_max is None and params.input.override_trusted_min is None:
      return

    detector = cspad_img.get_detector()
    for cached_range, panel in zip(self.cached_ranges, detector):
      panel.set_trusted_range(cached_range)

if __name__ == "__main__":
  from dials.util import halraiser
  try:
    script = InMemScript()
    script.run()
  except Exception as e:
    halraiser(e)
