from __future__ import division
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.xtc_dump
#
import psana
from xfel.cftbx.detector import cspad_cbf_tbx
from xfel.cxi.cspad_ana import cspad_tbx
import os, sys
import libtbx.load_env
from libtbx.utils import Sorry, Usage
from dials.util.options import OptionParser
from libtbx.phil import parse
from libtbx import easy_pickle

phil_scope = parse('''
  dispatch {
    max_events = None
      .type = int
      .help = If not specified, process all events. Otherwise, only process this many
    hit_finder = False
      .type = bool
      .help = If the number of strong reflections is less than \
              refinement.reflections.minimum_number_of_reflections, hit_filter=True \
              will discard the image. hit_finder=False: process all images
  }
  input {
    experiment = None
      .type = str
      .help = Experiment identifier, e.g. cxi84914
    run_num = None
      .type = int
      .help = Run number or run range to process
    address = None
      .type = str
      .help = Detector address, e.g. CxiDs2.0:Cspad.0 or detector alias, e.g. Ds1CsPad
    override_energy = None
      .type = float
      .help = If not None, use the input energy for every event instead of the energy \
              from the XTC stream
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
    }
  }
  output {
    output_dir = .
      .type = str
      .help = Directory output files will be placed
  }
''', process_includes=True)

class Script(object):
  """ Script to process dump XFEL data at LCLS """
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

    # Save the paramters
    self.params = params
    self.options = options

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank() # each process in MPI has a unique id, 0-indexed
    size = comm.Get_size() # size: number of processes running in this job

    # set up psana
    if params.format.file_format == "pickle":
      psana.setConfigFile(params.format.pickle.cfg)

    dataset_name = "exp=%s:run=%s:idx"%(params.input.experiment,params.input.run_num)
    ds = psana.DataSource(dataset_name)

    if params.format.file_format == "cbf":
      src = psana.Source('DetInfo(%s)'%params.input.address)
      psana_det = psana.Detector(src, ds.env())

    # set this to sys.maxint to analyze all events
    if params.dispatch.max_events is None:
      max_events = sys.maxint
    else:
      max_events = params.dispatch.max_events

    for run in ds.runs():
      # load a header only cspad cbf from the slac metrology
      base_dxtbx = cspad_cbf_tbx.env_dxtbx_from_slac_metrology(run, src)
      if base_dxtbx is None:
        raise Sorry("Couldn't load calibration file for run %d"%run.run())

      # list of all events
      times = run.times()
      nevents = min(len(times),max_events)
      # chop the list into pieces, depending on rank.  This assigns each process
      # events such that the get every Nth event where N is the number of processes
      mytimes = [times[i] for i in xrange(nevents) if (i+rank)%size == 0]

      for i in xrange(len(mytimes)):
        evt = run.event(mytimes[i])
        id = evt.get(psana.EventId)
        print "Event #",i," has id:",id

        timestamp = cspad_tbx.evt_timestamp(cspad_tbx.evt_time(evt)) # human readable format
        if timestamp is None:
          print "No timestamp, skipping shot"
          continue
        t = timestamp
        s = t[0:4] + t[5:7] + t[8:10] + t[11:13] + t[14:16] + t[17:19] + t[20:23]
        print "Processing shot", s

        if params.format.file_format == "pickle":
          if evt.get("skip_event"):
            print "Skipping event",id
            continue
          # the data needs to have already been processed and put into the event by psana
          data = evt.get(params.format.pickle.out_key)
          if data is None:
            print "No data"
            continue

          # set output paths according to the templates
          path = os.path.join(params.output.output_dir, "shot-" + s + ".pickle")

          print "Saving", path
          easy_pickle.dump(path, data)

        elif params.format.file_format == "cbf":
          # get numpy array, 32x185x388
          data = psana_det.calib(evt) # applies psana's complex run-dependent calibrations

          distance = cspad_tbx.env_distance(params.input.address, run.env(), params.format.cbf.detz_offset)
          if distance is None:
            print "No distance, skipping shot"
            continue

          if self.params.input.override_energy is None:
            wavelength = cspad_tbx.evt_wavelength(evt)
            if wavelength is None:
              print "No wavelength, skipping shot"
              continue
          else:
            wavelength = 12398.4187/self.params.input.override_energy

          # stitch together the header, data and metadata into the final dxtbx format object
          cspad_img = cspad_cbf_tbx.format_object_from_data(base_dxtbx, data, distance, wavelength, timestamp, params.input.address)
          path = os.path.join(params.output.output_dir, "shot-" + s + ".cbf")
          print "Saving", path

          # write the file
          import pycbf
          cspad_img._cbf_handle.write_widefile(path, pycbf.CBF,\
            pycbf.MIME_HEADERS|pycbf.MSG_DIGEST|pycbf.PAD_4K, 0)

      run.end()
    ds.end()

if __name__ == "__main__":
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
