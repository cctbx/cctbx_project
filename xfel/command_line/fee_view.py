from __future__ import division
from __future__ import print_function
from six.moves import range
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.fee_view
#
import sys, os
import psana
from psmon import publish
from psmon.plots import Image, XYPlot
from libtbx.phil import parse
from xfel.cxi.spectra_filter import spectra_filter
from xfel.cxi.cspad_ana import cspad_tbx

"""
Example usage:
cd ~/myrelease
sit_setup
psplot SPECTRUM FEE &
cxi.fee_view params.phil selected_filter=best

For full output, use:
psplot SPECTRUM FEE DC_OFFSET ALL_FEE ALL_FEE_RAW &
"""
phil_scope = parse("""
  selected_filter = None
    .type = str
    .help = Which filter to apply
  runs = None
    .type = str
    .help = Set of runs to display, eg 96 or 95-114 or 95,97,99-101, etc.
  experiment = None
    .type = str
    .help = Experiment name, e.g. cxid9114
  skip_events = None
    .type = int
    .help = Optionally skip N events
  include scope xfel.cxi.spectra_filter.phil_scope
""", process_includes = True)

def run(args):
  user_phil = []
  for arg in args:
    if os.path.isfile(arg):
      try:
        user_phil.append(parse(file_name=arg))
      except Exception as e:
        print(str(e))
        raise Sorry("Couldn't parse phil file %s"%arg)
    else:
      try:
        user_phil.append(parse(arg))
      except Exception as e:
        print(str(e))
        raise Sorry("Couldn't parse argument %s"%arg)
  params = phil_scope.fetch(sources=user_phil).extract()

  # cxid9114, source fee: FeeHxSpectrometer.0:Opal1000.1, downstream: CxiDg3.0:Opal1000.0
  # cxig3614, source fee: FeeHxSpectrometer.0:OrcaFl40.0

  src = psana.Source(params.spectra_filter.detector_address)
  dataset_name = "exp=%s:run=%s:idx"%(params.experiment, params.runs)
  print("Dataset string:", dataset_name)
  ds = psana.DataSource(dataset_name)
  spf = spectra_filter(params)

  if params.selected_filter == None:
    filter_name = params.spectra_filter.filter[0].name
  else:
    filter_name = params.selected_filter

  rank = 0
  size = 1
  max_events = sys.maxint

  for run in ds.runs():
    print("starting run", run.run())
    # list of all events
    times = run.times()

    if params.skip_events is not None:
      times = times[params.skip_events:]

    nevents = min(len(times),max_events)

    # chop the list into pieces, depending on rank.  This assigns each process
    # events such that the get every Nth event where N is the number of processes
    mytimes = [times[i] for i in range(nevents) if (i+rank)%size == 0]

    for i, t in enumerate(mytimes):
      evt = run.event(t)
      accepted, data, spectrum, dc_offset, all_data, all_data_raw = spf.filter_event(evt, filter_name)

      if not accepted:
        continue

      print(cspad_tbx.evt_timestamp(cspad_tbx.evt_time(evt)), "Publishing data for event", i)

      #header = "Event %d, m/f: %7.7f, f: %d"%(i, peak_max/flux, flux)
      header = "Event %d"%(i)

      if rank == 0:
        fee = Image(header, "FEE", data) # make a 2D plot
        publish.send("FEE", fee) # send to the display

        spectrumplot = XYPlot(header, 'summed 1D trace', list(range(data.shape[1])), spectrum) # make a 1D plot
        publish.send("SPECTRUM", spectrumplot) # send to the display

        fee = Image(header, "DC_OFFSET", dc_offset) # make a 2D plot
        publish.send("DC_OFFSET", fee) # send to the display
        fee = Image(header, "ALL_FEE", all_data) # make a 2D plot
        publish.send("ALL_FEE", fee) # send to the display
        fee = Image(header, "ALL_FEE_RAW", all_data_raw) # make a 2D plot
        publish.send("ALL_FEE_RAW", fee) # send to the display

if __name__ == "__main__":
  run(sys.argv[1:])
