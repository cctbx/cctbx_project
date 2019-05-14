from __future__ import division
from __future__ import print_function
from six.moves import range
from psana import *
import numpy as np
from libtbx import easy_pickle
import iotbx.phil, sys
import libtbx.load_env
from libtbx.utils import Sorry, Usage

master_phil = """
  dispatch{
    events_begin = None
    .type = int
    .help = If not specified, process all events. Otherwise, process events beginning at this number.
    events_end = None
    .type = int
    .help = If not specified, process all events. Otherwise, process events ending at this number.
    max_events = None
    .type = int
    .help = If not specified, process all events. Otherwise, only process this many
    events_accepted = False
    .type = bool
    .help = Plot average of filtered events
    events_rejected = False
    .type = bool
    .help = Plot average of rejected events
    events_all = False
    .type = bool
    .help = Plot average of all events
  }
  input {
    cfg = None
    .type = str
    .help = Path to psana config file
    experiment = None
    .type = str
    help = Experiment identifier, e.g. cxi84914
    run_num = None
    .type = int
    .help = Run number or run range to process
    address = None
    .type = str
    .help = FEE detector address, e.g. FEE-SPEC0
    dark = None
    .type = str
    .help = Path to FEE dark pickle file
    pixel_to_eV{
      energy_per_px = None
      .type = float
      .help = Energy per pixel conversion if known
      x_coord_one = None
      .type = int
      .help = Pixel valued x coordinate of known energy y.
       y_coord_one = None
      .type = int
      .help = Energy in eV of y coordinate of known x pixel position.
      x_coord_two = None
      .type = int
      .help = Pixel valued x coordinate of known energy y.
       y_coord_two = None
      .type = int
      .help = Energy in eV of y coordinate of known x pixel position.
    }
  }
  output {
    output_dir = .
    .type = str
    .help = Directory output files will be placed
  }
"""

def run(args):
  phil = iotbx.phil.process_command_line(args=args, master_string=master_phil).show()
  usage = \
  """ %s input.experiment=experimentname input.run_num=N input.address=address
  """%libtbx.env.dispatcher_name

  params = phil.work.extract()
  if not os.path.exists(params.output.output_dir):
    raise Sorry("Output path not found:" + params.output.output_dir)

  if params.input.experiment is None or \
  params.input.run_num is None or \
  params.input.address is None:
    raise Usage(usage)
  # set up psana
  if params.dispatch.events_accepted or params.dispatch.events_rejected:
    assert params.input.cfg is not None
    setConfigFile(params.input.cfg)

  dataset_name = "exp=%s:run=%s:idx"%(params.input.experiment,params.input.run_num)
  ds = DataSource(dataset_name)
  src = Source('DetInfo(%s)'%params.input.address)
  # set up multiprocessing with MPI
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank() # each process in MPI has a unique id, 0-indexed
  size = comm.Get_size() # size: number of processes running in this job

  if params.dispatch.max_events is None:
    max_events = sys.maxsize
  else:
    max_events = params.dispatch.max_events
  if params.input.dark is not None:
    dark = easy_pickle.load('%s'%params.input.dark)
  for run in ds.runs():
    times = run.times()
    if (params.dispatch.events_begin is None and params.dispatch.events_end is None):
      times = times[:]
    elif (params.dispatch.events_begin is not None and params.dispatch.events_end is None):
      times = times[params.dispatch.events_begin:]
    elif (params.dispatch.events_begin is None and params.dispatch.events_end is not None):
      times = times[:params.dispatch.events_end]
    elif (params.dispatch.events_begin is not None and params.dispatch.events_end is not None):
      times = times[params.dispatch.events_begin:params.dispatch.events_end]
    nevents = min(len(times),max_events)
  # chop the list into pieces, depending on rank.  This assigns each process
  # events such that the get every Nth event where N is the number of processes
    mytimes = [times[i] for i in range(nevents) if (i+rank)%size == 0]
    print(len(mytimes))
    #mytimes = mytimes[len(mytimes)-1000:len(mytimes)]
    totals = np.array([0.0])
    print("initial totals", totals)

    for i, t in enumerate(mytimes):
      print("Event", i, "of", len(mytimes), end=' ')
      evt = run.event(t)
      if params.dispatch.events_accepted or params.dispatch.events_all:
        if evt.get("skip_event")==True:
          continue
      elif params.dispatch.events_rejected:
        if evt.get("skip_event")==False:
          continue
      try:
        data = evt.get(Camera.FrameV1,src)
      except ValueError as e:
        src = Source('BldInfo(%s)'%params.input.address)
        data = evt.get(Bld.BldDataSpectrometerV1, src)
      if data is None:
        print("No data")
        continue
      #set default to determine FEE data type
      two_D=False
      #check attribute of data for type
      try:
        data = np.array(data.data16().astype(np.int32))
        two_D=True
      except AttributeError as e:
        data = np.array(data.hproj().astype(np.float64))

      if two_D:
        if 'dark' in locals():
          data = data - dark
        one_D_data = np.sum(data,0)/data.shape[0]
        two_D_data = np.double(data)
      else:
      #used to fix underflow problem that was present in earlier release of psana and pressent for LH80
        for i in range(len(data)):
          if data[i]>1000000000:
            data[i]=data[i]-(2**32)
        if 'dark' in locals():
          data = data - dark
        one_D_data = data

      totals[0] += 1
      print("total good:", totals[0])

      if not 'fee_one_D' in locals():
        fee_one_D = one_D_data
      else:
        fee_one_D += one_D_data
      if ('two_D_data' in locals() and not 'fee_two_D' in locals()):
        fee_two_D = two_D_data
      elif 'fee_two_D' in locals():
        fee_two_D += two_D_data

    acceptedtotals = np.zeros(totals.shape)
    acceptedfee1 = np.zeros((fee_one_D.shape))
    if 'fee_two_D' in locals():
      acceptedfee2 = np.zeros((fee_two_D.shape))
    print("Synchronizing rank", rank)
  comm.Reduce(fee_one_D,acceptedfee1)
  comm.Reduce(totals,acceptedtotals)
  if 'acceptedfee2' in locals():
    comm.Reduce(fee_two_D,acceptedfee2)
  print("number averaged", acceptedtotals[0])
  if rank == 0:
    if acceptedtotals[0] > 0:
      acceptedfee1 /= acceptedtotals[0]
      if 'acceptedfee2' in locals():
        acceptedfee2 /= acceptedtotals[0]

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from pylab import savefig,close
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    from matplotlib import cm

    if params.dispatch.events_accepted:
      easy_pickle.dump(os.path.join(params.output.output_dir,"fee_avg_1_D_"+'r%s'%params.input.run_num+"_accepted.pickle"), acceptedfee1)
      pp1 = PdfPages(os.path.join(params.output.output_dir,"fee_avg_1_D_"+'r%s'%params.input.run_num+"_accepted.pdf"))
      if 'acceptedfee2' in locals():
        easy_pickle.dump(os.path.join(params.output.output_dir,"fee_avg_2_D_"+'r%s'%params.input.run_num+"_accepted.pickle"), acceptedfee2)
        pp2 = PdfPages(os.path.join(params.output.output_dir,"fee_avg_2_D_"+'r%s'%params.input.run_num+"_accepted.pdf"))
    if params.dispatch.events_all:
      easy_pickle.dump(os.path.join(params.output.output_dir,"fee_avg_1_D_"+'r%s'%params.input.run_num+"_all.pickle"), acceptedfee1)
      pp1 = PdfPages(os.path.join(params.output.output_dir,"fee_avg_1_D_"+'r%s'%params.input.run_num+"_all.pdf"))
      if 'acceptedfee2' in locals():
        easy_pickle.dump(os.path.join(params.output.output_dir,"fee_avg_2_D_"+'r%s'%params.input.run_num+"_all.pickle"), acceptedfee2)
        pp2 = PdfPages(os.path.join(params.output.output_dir,"fee_avg_2_D_"+'r%s'%params.input.run_num+"_all.pdf"))
    if params.dispatch.events_rejected:
      easy_pickle.dump(os.path.join(params.output.output_dir,"fee_avg_1_D_"+'r%s'%params.input.run_num+"_rejected.pickle"), acceptedfee1)
      pp1 = PdfPages(os.path.join(params.output.output_dir,"fee_avg_1_D_"+'r%s'%params.input.run_num+"_rejected.pdf"))
      if 'acceptedfee2' in locals():
        easy_pickle.dump(os.path.join(params.output.output_dir,"fee_avg_2_D_"+'r%s'%params.input.run_num+"_rejected.pickle"), acceptedfee2)
        pp2 = PdfPages(os.path.join(params.output.output_dir,"fee_avg_2_D_"+'r%s'%params.input.run_num+"_rejected.pdf"))
    print("Done")
   #plotting result
   # matplotlib needs a different backend when run on the cluster nodes at SLAC
    # these two lines not needed when working interactively at SLAC, or on mac or on viper

    if params.input.pixel_to_eV.energy_per_px is not None:
      xvals = (np.array(range(acceptedfee1.shape[0]))-params.input.pixel_to_eV.x_coord_one)*params.input.pixel_to_eV.energy_per_px+params.input.pixel_to_eV.y_coord_one
      xvals = xvals[::-1]

    if params.input.pixel_to_eV.x_coord_two is not None:
      eV_per_px = (params.input.pixel_to_eV.y_coord_two-params.input.pixel_to_eV.y_coord_one)/(params.input.pixel_to_eV.x_coord_two-params.input.pixel_to_eV.x_coord_one)
      xvals = (np.array(range(acceptedfee1.shape[0]))-params.input.pixel_to_eV.x_coord_one)*eV_per_px+params.input.pixel_to_eV.y_coord_one
      xvals = xvals[::-1]

    if params.input.pixel_to_eV.x_coord_two is None and params.input.pixel_to_eV.energy_per_px is None:
      xvals=np.arange(0,len(acceptedfee1),1)

    yvals = acceptedfee1
    def OneD_plot(X,Y):
      plt.figure()
      plt.clf()
      plt.plot(X,Y)
      if params.dispatch.events_accepted:
        plt.title('Accepted Shots FEE Spectrum Run %s'%params.input.run_num)
      elif params.dispatch.events_all:
        plt.title('All Shots FEE Spectrum Run %s'%params.input.run_num)
      elif params.dispatch.events_rejected:
        plt.title('Rejected Shots FEE Spectrum Run %s'%params.input.run_num)
      if params.input.pixel_to_eV.x_coord_one is not None:
        plt.xlabel('eV', fontsize = 13)
      else:
        plt.xlabel('pixels', fontsize = 13)
      plt.ylabel('pixels', fontsize = 13)
      pp1.savefig()

    def TwoD_plot(data):
      plt.figure()
      ax = plt.gca()
     #  use specified range 0, 50 to plot runs 117 - 201
     #min=0, vmax=50
      cax=ax.imshow(data, interpolation='nearest',origin='lower',cmap=cm.coolwarm)
      plt.colorbar(cax, fraction=0.014, pad=0.04)
      if params.dispatch.events_accepted:
        ax.set_title('Accepted 2-D FEE Spectrum Run %s'%params.input.run_num)
      elif params.dispatch.events_all:
        ax.set_title('All 2-D FEE Spectrum Run %s'%params.input.run_num)
      elif params.dispatch.events_rejected:
        ax.set_title('Rejected 2-D FEE Spectrum Run %s'%params.input.run_num)
        pp2.savefig()

    OneD_plot(xvals,yvals)
    pp1.close()
    if 'acceptedfee2' in locals():
      TwoD_plot(acceptedfee2)
      pp2.close()



if __name__ == "__main__":
  run(sys.argv[1:])
