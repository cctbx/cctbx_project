from __future__ import division
from __future__ import print_function
from six.moves import range
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.xtc_dump
#
import psana
from xfel.cftbx.detector import cspad_cbf_tbx
from xfel.cxi.cspad_ana import cspad_tbx, rayonix_tbx
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
    selected_events = False
      .type = bool
      .help = If True, only dump events specified in input.event scopes
  }
  input {
    cfg = None
      .type = str
      .help = Path to psana config file. Genearlly not needed for CBFs. For image pickles, \
              the psana config file should have a mod_image_dict module.
    experiment = None
      .type = str
      .help = Experiment identifier, e.g. cxi84914
    run_num = None
      .type = int
      .help = Run number or run range to process
    address = None
      .type = str
      .help = Detector address, e.g. CxiDs2.0:Cspad.0 or detector alias, e.g. Ds1CsPad
    calib_dir = None
      .type = str
      .help = Non-standard calib directory location
    xtc_dir = None
      .type = str
      .help = Non-standard xtc directory location
    timestamp = None
      .type = str
      .multiple = True
      .help = Event timestamp(s) of event(s) in human-readable format of images to
      .help = dump (must also specify dispatch.selected_events=True.)
  }
  format {
    file_format = *cbf pickle
      .type = choice
      .help = Output file format, 64 tile segmented CBF or image pickle
    pickle {
      out_key = cctbx.xfel.image_dict
        .type = str
        .help = Key name that mod_image_dict uses to put image data in each psana event
    }
    cbf {
      detz_offset = None
        .type = float
        .help = Distance from back of detector rail to sample interaction region (CXI) \
                or actual detector distance (XPP/MFX)
      override_energy = None
        .type = float
        .help = If not None, use the input energy for every event instead of the energy \
                from the XTC stream
      mode = *cspad rayonix
        .type = choice
      cspad {
        gain_mask_value = None
          .type = float
          .help = If not None, use the gain mask for the run to multiply the low-gain pixels by this number
      }
      rayonix {
        bin_size = 2
          .type = int
          .help = Detector binning mode
        override_beam_x = None
          .type = float
          .help = If set, override the beam X position
        override_beam_y = None
          .type = float
          .help = If set, override the beam Y position
      }
    }
  }
  output {
    output_dir = .
      .type = str
      .help = Directory output files will be placed
    tmp_output_dir = None
      .type = str
      .help = Directory for CBFlib tmp output files
  }
''', process_includes=True)

class Script(object):
  """ Script to process dump XFEL data at LCLS """
  def __init__(self):
    """ Set up the option parser. Arguments come from the command line or a phil file """
    self.usage = """
%s input.experiment=experimentname input.run_num=N input.address=address
 format.file_format=cbf format.cbf.detz_offset=N
%s input.experiment=experimentname input.run_num=N input.address=address
 format.file_format=pickle format.pickle.cfg=path
    """%(libtbx.env.dispatcher_name, libtbx.env.dispatcher_name)

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
      if params.input.cfg is None:
        raise Usage(self.usage)
    else:
      raise Usage(self.usage)

    if not os.path.exists(params.output.output_dir):
      raise Sorry("Output path not found:" + params.output.output_dir)

    #Environment variable redirect for CBFLib temporary CBF_TMP_XYZ file output
    if params.format.file_format == "cbf":
      if params.output.tmp_output_dir is None:
        tmp_dir = os.path.join(params.output.output_dir, '.tmp')
      else:
        tmp_dir = os.path.join(params.output.tmp_output_dir, '.tmp')
      if not os.path.exists(tmp_dir):
        try:
          os.makedirs(tmp_dir)
        except Exception as e:
          if not os.path.exists(tmp_dir):
            halraiser(e)
      os.environ['CBF_TMP_DIR'] = tmp_dir

    # Save the paramters
    self.params = params
    self.options = options

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank() # each process in MPI has a unique id, 0-indexed
    size = comm.Get_size() # size: number of processes running in this job

    # set up psana
    if params.input.cfg is not None:
      psana.setConfigFile(params.input.cfg)

    if params.input.calib_dir is not None:
      psana.setOption('psana.calib-dir',params.input.calib_dir)

    dataset_name = "exp=%s:run=%s:idx"%(params.input.experiment,params.input.run_num)
    if params.input.xtc_dir is not None:
      dataset_name = "exp=%s:run=%s:idx:dir=%s"%(params.input.experiment,params.input.run_num,params.input.xtc_dir)

    ds = psana.DataSource(dataset_name)

    if params.format.file_format == "cbf":
      src = psana.Source('DetInfo(%s)'%params.input.address)
      psana_det = psana.Detector(params.input.address, ds.env())

    # set this to sys.maxint to analyze all events
    if params.dispatch.max_events is None:
      max_events = sys.maxint
    else:
      max_events = params.dispatch.max_events

    for run in ds.runs():
      if params.format.file_format == "cbf":
        if params.format.cbf.mode == "cspad":
          # load a header only cspad cbf from the slac metrology
          base_dxtbx = cspad_cbf_tbx.env_dxtbx_from_slac_metrology(run, params.input.address)
          if base_dxtbx is None:
            raise Sorry("Couldn't load calibration file for run %d"%run.run())
        elif params.format.cbf.mode == "rayonix":
          # load a header only rayonix cbf from the input parameters
          detector_size = rayonix_tbx.get_rayonix_detector_dimensions(ds.env())
          base_dxtbx = rayonix_tbx.get_dxtbx_from_params(params.format.cbf.rayonix, detector_size)

      # list of all events
      times = run.times()
      if params.dispatch.selected_events:
        times = [t for t in times if cspad_tbx.evt_timestamp((t.seconds(),t.nanoseconds()/1e6)) in params.input.timestamp]
      nevents = min(len(times),max_events)
      # chop the list into pieces, depending on rank.  This assigns each process
      # events such that the get every Nth event where N is the number of processes
      mytimes = [times[i] for i in range(nevents) if (i+rank)%size == 0]

      for i in range(len(mytimes)):
        evt = run.event(mytimes[i])
        id = evt.get(psana.EventId)
        print("Event #",i," has id:",id)

        timestamp = cspad_tbx.evt_timestamp(cspad_tbx.evt_time(evt)) # human readable format
        if timestamp is None:
          print("No timestamp, skipping shot")
          continue

        if evt.get("skip_event") or "skip_event" in [key.key() for key in evt.keys()]:
          print("Skipping event",timestamp)
          continue

        t = timestamp
        s = t[0:4] + t[5:7] + t[8:10] + t[11:13] + t[14:16] + t[17:19] + t[20:23]
        print("Processing shot", s)

        if params.format.file_format == "pickle":
          if evt.get("skip_event"):
            print("Skipping event",id)
            continue
          # the data needs to have already been processed and put into the event by psana
          data = evt.get(params.format.pickle.out_key)
          if data is None:
            print("No data")
            continue

          # set output paths according to the templates
          path = os.path.join(params.output.output_dir, "shot-" + s + ".pickle")

          print("Saving", path)
          easy_pickle.dump(path, data)

        elif params.format.file_format == "cbf":
          if params.format.cbf.mode == "cspad":
            # get numpy array, 32x185x388
            data = cspad_cbf_tbx.get_psana_corrected_data(psana_det, evt, use_default=False, dark=True,
                                                          common_mode=None,
                                                          apply_gain_mask=params.format.cbf.cspad.gain_mask_value is not None,
                                                          gain_mask_value=params.format.cbf.cspad.gain_mask_value,
                                                          per_pixel_gain=False)

            distance = cspad_tbx.env_distance(params.input.address, run.env(), params.format.cbf.detz_offset)
          elif params.format.cbf.mode == "rayonix":
            data = rayonix_tbx.get_data_from_psana_event(evt, params.input.address)
            distance = params.format.cbf.detz_offset

          if distance is None:
            print("No distance, skipping shot")
            continue

          if self.params.format.cbf.override_energy is None:
            wavelength = cspad_tbx.evt_wavelength(evt)
            if wavelength is None:
              print("No wavelength, skipping shot")
              continue
          else:
            wavelength = 12398.4187/self.params.format.cbf.override_energy

          # stitch together the header, data and metadata into the final dxtbx format object
          if params.format.cbf.mode == "cspad":
            image = cspad_cbf_tbx.format_object_from_data(base_dxtbx, data, distance, wavelength, timestamp, params.input.address)
          elif params.format.cbf.mode == "rayonix":
            image = rayonix_tbx.format_object_from_data(base_dxtbx, data, distance, wavelength, timestamp, params.input.address)
          path = os.path.join(params.output.output_dir, "shot-" + s + ".cbf")
          print("Saving", path)

          # write the file
          import pycbf
          image._cbf_handle.write_widefile(path, pycbf.CBF,\
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
