from __future__ import absolute_import, division, print_function
from six.moves import range
#-*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.mpi_average
#

import psana
import numpy as np
from xfel.cxi.cspad_ana import cspad_tbx, parse_calib
import libtbx
from libtbx import easy_pickle
import libtbx.option_parser
from scitbx.array_family import flex
import sys, os
from libtbx.utils import Sorry
from six.moves import zip

def average(argv=None):
  if argv == None:
    argv = sys.argv[1:]

  try:
    from mpi4py import MPI
  except ImportError:
    raise Sorry("MPI not found")

  command_line = (libtbx.option_parser.option_parser(
    usage="""
%s [-p] -c config -x experiment -a address -r run -d detz_offset [-o outputdir] [-A averagepath] [-S stddevpath] [-M maxpath] [-n numevents] [-s skipnevents] [-v] [-m] [-b bin_size] [-X override_beam_x] [-Y override_beam_y] [-D xtc_dir] [-f] [-g gain_mask_value] [--min] [--minpath minpath]

To write image pickles use -p, otherwise the program writes CSPAD CBFs.
Writing CBFs requires the geometry to be already deployed.

Examples:
cxi.mpi_average -c cxi49812/average.cfg -x cxi49812 -a CxiDs1.0:Cspad.0 -r 25 -d 571

Use one process on the current node to process all the events from run 25 of
experiment cxi49812, using a detz_offset of 571.

mpirun -n 16 cxi.mpi_average -c cxi49812/average.cfg -x cxi49812 -a CxiDs1.0:Cspad.0 -r 25 -d 571

As above, using 16 cores on the current node.

bsub -a mympi -n 100 -o average.out -q psanaq cxi.mpi_average -c cxi49812/average.cfg -x cxi49812 -a CxiDs1.0:Cspad.0 -r 25 -d 571 -o cxi49812

As above, using the psanaq and 100 cores, putting the log in average.out and
the output images in the folder cxi49812.
""" % libtbx.env.dispatcher_name)
                .option(None, "--as_pickle", "-p",
                        action="store_true",
                        default=False,
                        dest="as_pickle",
                        help="Write results as image pickle files instead of cbf files")
                .option(None, "--raw_data", "-R",
                        action="store_true",
                        default=False,
                        dest="raw_data",
                        help="Disable psana corrections such as dark pedestal subtraction or common mode (cbf only)")
                .option(None, "--background_pickle", "-B",
                        default=None,
                        dest="background_pickle",
                        help="")
                .option(None, "--config", "-c",
                        type="string",
                        default=None,
                        dest="config",
                        metavar="PATH",
                        help="psana config file")
                .option(None, "--experiment", "-x",
                        type="string",
                        default=None,
                        dest="experiment",
                        help="experiment name (eg cxi84914)")
                .option(None, "--run", "-r",
                        type="int",
                        default=None,
                        dest="run",
                        help="run number")
                .option(None, "--address", "-a",
                        type="string",
                        default="CxiDs2.0:Cspad.0",
                        dest="address",
                        help="detector address name (eg CxiDs2.0:Cspad.0)")
                .option(None, "--detz_offset", "-d",
                        type="float",
                        default=None,
                        dest="detz_offset",
                        help="offset (in mm) from sample interaction region to back of CSPAD detector rail (CXI), or detector distance (XPP)")
                .option(None, "--outputdir", "-o",
                        type="string",
                        default=".",
                        dest="outputdir",
                        metavar="PATH",
                        help="Optional path to output directory for output files")
                .option(None, "--averagebase", "-A",
                        type="string",
                        default="{experiment!l}_avg-r{run:04d}",
                        dest="averagepath",
                        metavar="PATH",
                        help="Path to output average image without extension. String substitution allowed")
                .option(None, "--stddevbase", "-S",
                        type="string",
                        default="{experiment!l}_stddev-r{run:04d}",
                        dest="stddevpath",
                        metavar="PATH",
                        help="Path to output standard deviation image without extension. String substitution allowed")
                .option(None, "--maxbase", "-M",
                        type="string",
                        default="{experiment!l}_max-r{run:04d}",
                        dest="maxpath",
                        metavar="PATH",
                        help="Path to output maximum projection image without extension. String substitution allowed")
                .option(None, "--numevents", "-n",
                        type="int",
                        default=None,
                        dest="numevents",
                        help="Maximum number of events to process. Default: all")
                .option(None, "--skipevents", "-s",
                        type="int",
                        default=0,
                        dest="skipevents",
                        help="Number of events in the beginning of the run to skip. Default: 0")
                .option(None, "--verbose", "-v",
                        action="store_true",
                        default=False,
                        dest="verbose",
                        help="Print more information about progress")
                .option(None, "--pickle-optical-metrology", "-m",
                        action="store_true",
                        default=False,
                        dest="pickle_optical_metrology",
                        help="If writing pickle files, use the optical metrology in the experiment's calib directory")
                .option(None, "--bin_size", "-b",
                        type="int",
                        default=None,
                        dest="bin_size",
                        help="Rayonix detector bin size")
                .option(None, "--override_beam_x", "-X",
                        type="float",
                        default=None,
                        dest="override_beam_x",
                        help="Rayonix detector beam center x coordinate")
                .option(None, "--override_beam_y", "-Y",
                        type="float",
                        default=None,
                        dest="override_beam_y",
                        help="Rayonix detector beam center y coordinate")
                .option(None, "--calib_dir", "-C",
                        type="string",
                        default=None,
                        dest="calib_dir",
                        metavar="PATH",
                        help="calibration directory")
                .option(None, "--pickle_calib_dir", "-P",
                        type="string",
                        default=None,
                        dest="pickle_calib_dir",
                        metavar="PATH",
                        help="pickle calibration directory specification. Replaces --calib_dir functionality.")
                .option(None, "--xtc_dir", "-D",
                        type="string",
                        default=None,
                        dest="xtc_dir",
                        metavar="PATH",
                        help="xtc stream directory")
                .option(None, "--use_ffb", "-f",
                        action="store_true",
                        default=False,
                        dest="use_ffb",
                        help="Use the fast feedback filesystem at LCLS. Only for the active experiment!")
                .option(None, "--gain_mask_value", "-g",
                        type="float",
                        default=None,
                        dest="gain_mask_value",
                        help="Ratio between low and high gain pixels, if CSPAD in mixed-gain mode. Only used in CBF averaging mode.")
                .option(None, "--min", None,
                        action="store_true",
                        default=False,
                        dest="do_minimum_projection",
                        help="Output a minimum projection")
                .option(None, "--minpath", None,
                        type="string",
                        default="{experiment!l}_min-r{run:04d}",
                        dest="minpath",
                        metavar="PATH",
                        help="Path to output minimum image without extension. String substitution allowed")
                ).process(args=argv)


  if len(command_line.args) > 0 or \
      command_line.options.as_pickle is None or \
      command_line.options.experiment is None or \
      command_line.options.run is None or \
      command_line.options.address is None or \
      command_line.options.detz_offset is None or \
      command_line.options.averagepath is None or \
      command_line.options.stddevpath is None or \
      command_line.options.maxpath is None or \
      command_line.options.pickle_optical_metrology is None:
    command_line.parser.show_help()
    return

  # set this to sys.maxint to analyze all events
  if command_line.options.numevents is None:
    maxevents = sys.maxsize
  else:
    maxevents = command_line.options.numevents

  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  size = comm.Get_size()

  if command_line.options.config is not None:
    psana.setConfigFile(command_line.options.config)
  dataset_name = "exp=%s:run=%d:idx"%(command_line.options.experiment, command_line.options.run)
  if command_line.options.xtc_dir is not None:
    if command_line.options.use_ffb:
      raise Sorry("Cannot specify the xtc_dir and use SLAC's ffb system")
    dataset_name += ":dir=%s"%command_line.options.xtc_dir
  elif command_line.options.use_ffb:
    # as ffb is only at SLAC, ok to hardcode /reg/d here
    dataset_name += ":dir=/reg/d/ffb/%s/%s/xtc"%(command_line.options.experiment[0:3],command_line.options.experiment)
  if command_line.options.calib_dir is not None:
    psana.setOption('psana.calib-dir',command_line.options.calib_dir)
  ds = psana.DataSource(dataset_name)
  address = command_line.options.address
  src = psana.Source('DetInfo(%s)'%address)
  nevent = np.array([0.])

  if command_line.options.background_pickle is not None:
    background = easy_pickle.load(command_line.options.background_pickle)['DATA'].as_numpy_array()

  for run in ds.runs():
    runnumber = run.run()

    if not command_line.options.as_pickle:
      psana_det = psana.Detector(address, ds.env())

    # list of all events
    if command_line.options.skipevents > 0:
      print("Skipping first %d events"%command_line.options.skipevents)
    elif "Rayonix" in command_line.options.address:
      print("Skipping first image in the Rayonix detector") # Shuttering issue
      command_line.options.skipevents = 1

    times = run.times()[command_line.options.skipevents:]
    nevents = min(len(times),maxevents)
    # chop the list into pieces, depending on rank.  This assigns each process
    # events such that the get every Nth event where N is the number of processes
    mytimes = [times[i] for i in range(nevents) if (i+rank)%size == 0]
    for i in range(len(mytimes)):
      if i%10==0: print('Rank',rank,'processing event',rank*len(mytimes)+i,', ',i,'of',len(mytimes))
      evt = run.event(mytimes[i])
      #print "Event #",rank*mylength+i," has id:",evt.get(EventId)
      if 'Rayonix' in command_line.options.address or 'FeeHxSpectrometer' in command_line.options.address or 'XrayTransportDiagnostic' in command_line.options.address:
        data = evt.get(psana.Camera.FrameV1,src)
        if data is None:
          print("No data")
          continue
        data=data.data16().astype(np.float64)
      elif command_line.options.as_pickle:
        data = evt.get(psana.ndarray_float64_3, src, 'image0')
      else:
        # get numpy array, 32x185x388
        from xfel.cftbx.detector.cspad_cbf_tbx import get_psana_corrected_data
        if command_line.options.raw_data:
          data = get_psana_corrected_data(psana_det, evt, use_default=False, dark=False, common_mode=None,
                                          apply_gain_mask=False, per_pixel_gain=False)
        else:
          if command_line.options.gain_mask_value is None:
            data = get_psana_corrected_data(psana_det, evt, use_default=True)
          else:
            data = get_psana_corrected_data(psana_det, evt, use_default=False, dark=True, common_mode=None,
                                            apply_gain_mask=True, gain_mask_value=command_line.options.gain_mask_value,
                                            per_pixel_gain=False)

      if data is None:
        print("No data")
        continue

      if command_line.options.background_pickle is not None:
        data -= background

      if 'FeeHxSpectrometer' in command_line.options.address or 'XrayTransportDiagnostic' in command_line.options.address:
        distance = np.array([0.0])
        wavelength = np.array([1.0])
      else:
        d = cspad_tbx.env_distance(address, run.env(), command_line.options.detz_offset)
        if d is None:
          print("No distance, using distance", command_line.options.detz_offset)
          assert command_line.options.detz_offset is not None
          if 'distance' not in locals():
            distance = np.array([command_line.options.detz_offset])
          else:
            distance += command_line.options.detz_offset
        else:
          if 'distance' in locals():
            distance += d
          else:
            distance = np.array([float(d)])

        w = cspad_tbx.evt_wavelength(evt)
        if w is None:
          print("No wavelength")
          if 'wavelength' not in locals():
            wavelength = np.array([1.0])
        else:
          if 'wavelength' in locals():
            wavelength += w
          else:
            wavelength = np.array([w])

      t = cspad_tbx.evt_time(evt)
      if t is None:
        print("No timestamp, skipping shot")
        continue
      if 'timestamp' in locals():
        timestamp += t[0] + (t[1]/1000)
      else:
        timestamp = np.array([t[0] + (t[1]/1000)])

      if 'sum' in locals():
        sum+=data
      else:
        sum=np.array(data, copy=True)
      if 'sumsq' in locals():
        sumsq+=data*data
      else:
        sumsq=data*data
      if 'maximum' in locals():
        maximum=np.maximum(maximum,data)
      else:
        maximum=np.array(data, copy=True)

      if command_line.options.do_minimum_projection:
        if 'minimum' in locals():
          minimum=np.minimum(minimum,data)
        else:
          minimum=np.array(data, copy=True)

      nevent += 1

  #sum the images across mpi cores
  if size > 1:
    print("Synchronizing rank", rank)
  totevent = np.zeros(nevent.shape)
  comm.Reduce(nevent,totevent)

  if rank == 0 and totevent[0] == 0:
    raise Sorry("No events found in the run")

  sumall = np.zeros(sum.shape).astype(sum.dtype)
  comm.Reduce(sum,sumall)

  sumsqall = np.zeros(sumsq.shape).astype(sumsq.dtype)
  comm.Reduce(sumsq,sumsqall)

  maxall = np.zeros(maximum.shape).astype(maximum.dtype)
  comm.Reduce(maximum,maxall, op=MPI.MAX)

  if command_line.options.do_minimum_projection:
    minall = np.zeros(maximum.shape).astype(minimum.dtype)
    comm.Reduce(minimum,minall, op=MPI.MIN)

  waveall = np.zeros(wavelength.shape).astype(wavelength.dtype)
  comm.Reduce(wavelength,waveall)

  distall = np.zeros(distance.shape).astype(distance.dtype)
  comm.Reduce(distance,distall)

  timeall = np.zeros(timestamp.shape).astype(timestamp.dtype)
  comm.Reduce(timestamp,timeall)

  if rank==0:
    if size > 1:
      print("Synchronized")

    # Accumulating floating-point numbers introduces errors,
    # which may cause negative variances.  Since a two-pass
    # approach is unacceptable, the standard deviation is
    # clamped at zero.
    mean = sumall / float(totevent[0])
    variance = (sumsqall / float(totevent[0])) - (mean**2)
    variance[variance < 0] = 0
    stddev = np.sqrt(variance)

    wavelength = waveall[0] / totevent[0]
    distance = distall[0] / totevent[0]
    pixel_size = cspad_tbx.pixel_size
    saturated_value = cspad_tbx.cspad_saturated_value
    timestamp = timeall[0] / totevent[0]
    timestamp = (int(timestamp), timestamp % int(timestamp) * 1000)
    timestamp = cspad_tbx.evt_timestamp(timestamp)


    if command_line.options.as_pickle:
      extension = ".pickle"
    else:
      extension = ".cbf"

    dest_paths = [cspad_tbx.pathsubst(command_line.options.averagepath + extension, evt, ds.env()),
                  cspad_tbx.pathsubst(command_line.options.stddevpath  + extension, evt, ds.env()),
                  cspad_tbx.pathsubst(command_line.options.maxpath     + extension, evt, ds.env())]
    if command_line.options.do_minimum_projection:
      dest_paths.append(cspad_tbx.pathsubst(command_line.options.minpath + extension, evt, ds.env()))

    dest_paths = [os.path.join(command_line.options.outputdir, path) for path in dest_paths]
    if 'Rayonix' in command_line.options.address:
      all_data = [mean, stddev, maxall]
      if command_line.options.do_minimum_projection:
        all_data.append(minall)
      from xfel.cxi.cspad_ana import rayonix_tbx
      pixel_size = rayonix_tbx.get_rayonix_pixel_size(command_line.options.bin_size)
      beam_center = [command_line.options.override_beam_x,command_line.options.override_beam_y]
      active_areas = flex.int([0,0,mean.shape[1],mean.shape[0]])
      split_address = cspad_tbx.address_split(address)
      old_style_address = split_address[0] + "-" + split_address[1] + "|" + split_address[2] + "-" + split_address[3]
      for data, path in zip(all_data, dest_paths):
        print("Saving", path)
        d = cspad_tbx.dpack(
            active_areas=active_areas,
            address=old_style_address,
            beam_center_x=pixel_size * beam_center[0],
            beam_center_y=pixel_size * beam_center[1],
            data=flex.double(data),
            distance=distance,
            pixel_size=pixel_size,
            saturated_value=rayonix_tbx.rayonix_saturated_value,
            timestamp=timestamp,
            wavelength=wavelength)
        easy_pickle.dump(path, d)
    elif 'FeeHxSpectrometer' in command_line.options.address or 'XrayTransportDiagnostic' in command_line.options.address:
      all_data = [mean, stddev, maxall]
      split_address = cspad_tbx.address_split(address)
      old_style_address = split_address[0] + "-" + split_address[1] + "|" + split_address[2] + "-" + split_address[3]
      if command_line.options.do_minimum_projection:
        all_data.append(minall)
      for data, path in zip(all_data, dest_paths):
        d = cspad_tbx.dpack(
          address = old_style_address,
          data = flex.double(data),
          distance = distance,
          pixel_size = 0.1,
          timestamp=timestamp,
          wavelength=wavelength
        )
        print("Saving", path)
        easy_pickle.dump(path, d)
    elif command_line.options.as_pickle:
      split_address = cspad_tbx.address_split(address)
      old_style_address = split_address[0] + "-" + split_address[1] + "|" + split_address[2] + "-" + split_address[3]

      xpp = 'xpp' in address.lower()
      if xpp:
        evt_time = cspad_tbx.evt_time(evt) # tuple of seconds, milliseconds
        timestamp = cspad_tbx.evt_timestamp(evt_time) # human readable format
        from iotbx.detectors.cspad_detector_formats import detector_format_version, reverse_timestamp
        from xfel.cxi.cspad_ana.cspad_tbx import xpp_active_areas
        version_lookup = detector_format_version(old_style_address, reverse_timestamp(timestamp)[0])
        assert version_lookup is not None
        active_areas = xpp_active_areas[version_lookup]['active_areas']
        beam_center = [1765 // 2, 1765 // 2]
      else:
        if command_line.options.pickle_calib_dir is not None:
          metro_path = command_line.options.pickle_calib_dir
        elif command_line.options.pickle_optical_metrology:
          from xfel.cftbx.detector.cspad_cbf_tbx import get_calib_file_path
          metro_path = get_calib_file_path(run.env(), address, run)
        else:
          metro_path = libtbx.env.find_in_repositories("xfel/metrology/CSPad/run4/CxiDs1.0_Cspad.0")
        sections = parse_calib.calib2sections(metro_path)
        beam_center, active_areas = cspad_tbx.cbcaa(
          cspad_tbx.getConfig(address, ds.env()), sections)

      class fake_quad(object):
        def __init__(self, q, d):
          self.q = q
          self.d = d

        def quad(self):
          return self.q

        def data(self):
          return self.d

      if xpp:
        quads = [fake_quad(i, mean[i*8:(i+1)*8,:,:]) for i in range(4)]
        mean = cspad_tbx.image_xpp(old_style_address, None, ds.env(), active_areas, quads = quads)
        mean = flex.double(mean.astype(np.float64))

        quads = [fake_quad(i, stddev[i*8:(i+1)*8,:,:]) for i in range(4)]
        stddev = cspad_tbx.image_xpp(old_style_address, None, ds.env(), active_areas, quads = quads)
        stddev = flex.double(stddev.astype(np.float64))

        quads = [fake_quad(i, maxall[i*8:(i+1)*8,:,:]) for i in range(4)]
        maxall = cspad_tbx.image_xpp(old_style_address, None, ds.env(), active_areas, quads = quads)
        maxall = flex.double(maxall.astype(np.float64))

        if command_line.options.do_minimum_projection:
          quads = [fake_quad(i, minall[i*8:(i+1)*8,:,:]) for i in range(4)]
          minall = cspad_tbx.image_xpp(old_style_address, None, ds.env(), active_areas, quads = quads)
          minall = flex.double(minall.astype(np.float64))
      else:
        quads = [fake_quad(i, mean[i*8:(i+1)*8,:,:]) for i in range(4)]
        mean = cspad_tbx.CsPadDetector(
          address, evt, ds.env(), sections, quads=quads)
        mean = flex.double(mean.astype(np.float64))

        quads = [fake_quad(i, stddev[i*8:(i+1)*8,:,:]) for i in range(4)]
        stddev = cspad_tbx.CsPadDetector(
          address, evt, ds.env(), sections, quads=quads)
        stddev = flex.double(stddev.astype(np.float64))

        quads = [fake_quad(i, maxall[i*8:(i+1)*8,:,:]) for i in range(4)]
        maxall = cspad_tbx.CsPadDetector(
          address, evt, ds.env(), sections, quads=quads)
        maxall = flex.double(maxall.astype(np.float64))

        if command_line.options.do_minimum_projection:
          quads = [fake_quad(i, minall[i*8:(i+1)*8,:,:]) for i in range(4)]
          minall = cspad_tbx.CsPadDetector(
            address, evt, ds.env(), sections, quads=quads)
          minall = flex.double(minall.astype(np.float64))

      all_data = [mean, stddev, maxall]
      if command_line.options.do_minimum_projection:
        all_data.append(minall)

      for data, path in zip(all_data, dest_paths):
        print("Saving", path)

        d = cspad_tbx.dpack(
          active_areas=active_areas,
          address=old_style_address,
          beam_center_x=pixel_size * beam_center[0],
          beam_center_y=pixel_size * beam_center[1],
          data=data,
          distance=distance,
          pixel_size=pixel_size,
          saturated_value=saturated_value,
          timestamp=timestamp,
          wavelength=wavelength)

        easy_pickle.dump(path, d)
    else:
      # load a header only cspad cbf from the slac metrology
      from xfel.cftbx.detector import cspad_cbf_tbx
      import pycbf
      base_dxtbx = cspad_cbf_tbx.env_dxtbx_from_slac_metrology(run, address)
      if base_dxtbx is None:
        raise Sorry("Couldn't load calibration file for run %d"%run.run())

      all_data = [mean, stddev, maxall]
      if command_line.options.do_minimum_projection:
        all_data.append(minall)

      for data, path in zip(all_data, dest_paths):
        print("Saving", path)
        cspad_img = cspad_cbf_tbx.format_object_from_data(base_dxtbx, data, distance, wavelength, timestamp, address, round_to_int=False)
        cspad_img._cbf_handle.write_widefile(path, pycbf.CBF,\
          pycbf.MIME_HEADERS|pycbf.MSG_DIGEST|pycbf.PAD_4K, 0)


if (__name__ == "__main__"):
  sys.exit(average(sys.argv[1:]))
