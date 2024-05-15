from __future__ import absolute_import, division, print_function
from six.moves import range

from psana import *
import sys
import numpy as np
from xfel.amo.pnccd_ana                 import pnccd_tbx
from xfel.amo.pnccd_ana                 import fxs
import matplotlib.pyplot as plt
from six.moves import zip

plt.ion()
########################################
# Due to the mask sometimes having zero values
# we're bound to get divisions with zeros at
#times. Here ignoring those errors.
np.seterr(divide='ignore', invalid='ignore')

from libtbx.mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def h5gen(run,timestamps = None, first = None, last = None):

    # Singel CPU
    if size == 1:
       nom   = rank
       denom = size
    # MPI
    else:
       nom   = rank - 1
       denom = size - 1


    times     = timestamps
    nevents   = len(times)
    mytimes,myevents  = zip(*[(times[i],i) for i in range(nevents) if (i+nom)%denom == 0])

    for j in range(len(mytimes)):
         yield myevents[j],mytimes[j]


def idxgen(run,timestamps = None, first = None, last = None):
    #print "idx mode"
    #  Use timestamps from index file
    if timestamps is  None:
       timestamps      = run.times()

    if first is None :
       first   = 0

    if last is None :
       last    = len(timestamps)
    else:
       last    = min(last,len(timestamps))      # Check that last time-stamp exists

    # Singel CPU
    if size == 1:
       nom   = rank
       denom = size
    # MPI
    else:
       nom   = rank - 1
       denom = size - 1


    times     = timestamps[first:last]
    nevents   = len(times)
    mytimes,myevents  = zip(*[(times[i],i) for i in range(nevents) if (i+nom)%denom == 0])

    for j in range(len(mytimes)):
         yield myevents[j],run.event(mytimes[j])



def smdgen(run,timestamps = None, first = None, last = None):
    #print "smd mode"
    if first is None :
       first   = 0

    if last is None :
       last    = 1e20                           # We typically don't know what the last events is. So for now use a large number


    # Singel CPU
    if size == 1:
       nom   = rank
       denom = size
    # MPI
    else:
       nom   = rank - 1
       denom = size - 1


    if timestamps is None :

       for nevent,evt in enumerate(run.events()):
           if   nevent <  first : continue
           elif nevent == last  : return
           elif nevent%denom == nom:
             yield nevent-first,evt

    else :  # Only applicable for xtc format

       ct = 0
       for nevent,evt in enumerate(run.events()):
            t = pnccd_tbx.get_psana_time(evt)
            # Check if event exists in timestamps
            if np.equal(t, timestamps).all(axis=1).any() :
               if   ct <  first : continue
               elif ct == last  : return
               elif ct%denom == nom:
                 yield ct,evt
               ct += 1


def compute_calib(argv=None) :

  """Function to compute the average image without applied geometry from FXS images
       extracted from xtc (smd,idx,xtc format) or h5 files.
       Works for Single CPU, Multi-Processor interactive jobs and MPI batch jobs

       For a definition of input arguments argv and batch processing instructions see  ***  mpi_fxs_launch.py ***

       compute_calib produces the following output files:

       * Index file : Information about the events processed including time-stamps, beam center, total intensity, particle size etc
       * Ave        : Average image, as 1D array, in cartesian coordinates without geometry applied

  """

  if argv == None:
    argv = sys.argv[1:]

  try:
     from libtbx.mpi4py import MPI
  except ImportError:
     raise Sorry("MPI not found")

  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  size = comm.Get_size()


  if argv.hit is None :
     hit        = -1.0e20       # Process everything
  else:
     hit        = argv.hit      # Process everything > hit

  ftype  = argv.ftype

  if argv.param_path is not None :
     if ftype == 'h5' :
        param_file            = np.genfromtxt(argv.param_path,skiprows=1,dtype=None)
        timestamps,filestamps = pnccd_tbx.get_h5_event(param_file)
     elif ftype == 'xtc' :
        param_file            = np.genfromtxt(argv.param_path,skiprows=1,dtype=None)
        timestamps            = pnccd_tbx.get_time(param_file)
     else :
        param_file            = np.genfromtxt(argv.param_path,skiprows=1)
        timestamps            = pnccd_tbx.get_psana_event(param_file)
  else:
     timestamps = None

  # The first and last events to processed
  first = argv.first
  last  = argv.last


  # Check data format

  if ftype == 'h5' :
       import h5py
       run          = int(argv.run)

       # Get time-stamps from all h5-files
       if argv.param_path is None :
          timestamps = []
          filestamps = []
          # Loop over all h5-files and store the time-stamps
          for i in os.listdir(argv.xtc_dir):
              if i.endswith(".h5"):
                 f  = h5py.File(i,'r')
                 filestamps.append(i[-7:-4])
                 timestamps.append(list(f.keys()))
                 continue
              else:
                 continue

       dataset_name = "%s-r%s"%(argv.experiment, str(argv.run).zfill(4)) # Ascert 4 digit run number
       exprun       = os.path.join(argv.xtc_dir,dataset_name)

       if argv.first is None :
          first   = 0

       if argv.last is None :
          last    = len(timestamps)
       else:
          last    = min(last,len(timestamps))      # Check that last time-stamp exists

       timestamps = timestamps[first:last]
       filestamps = filestamps[first:last]


       evtgen       = h5gen

  else :

       exprun = "exp=%s:run=%d"%(argv.experiment, argv.run)
       if (ftype == 'xtc') :
           dataset_name = exprun+':xtc'
       elif  (ftype == 'idx') :
           dataset_name = exprun+':idx'
       elif(ftype == 'idx_ffb') :
           dataset_name = exprun+':idx'
           # as ffb is only at SLAC, ok to hardcode /reg/d here
           dataset_name += ":dir=/reg/d/ffb/%s/%s/xtc"%(argv.experiment[0:3],argv.experiment)
       elif(ftype == 'smd') :
           dataset_name = exprun+':smd'
       elif(ftype == 'smd_ffb') :
           dataset_name = exprun+':smd'
           # as ffb is only at SLAC, ok to hardcode /reg/d here ADD live!
           dataset_name += ":dir=/reg/d/ffb/%s/%s/xtc:live"%(argv.experiment[0:3],argv.experiment)
           exprun = dataset_name

       ds           = DataSource(dataset_name)
       run          = next(ds.runs())

       # Select event generator
       if    (ftype=='smd') or (ftype == 'smd_ffb') or (ftype == 'xtc'):
         evtgen = smdgen
       elif  (ftype=='idx') or (ftype == 'idx_ffb'):
         evtgen = idxgen


  if size == 1:
     plot = argv.plot
  else:
     plot = 0


  FXS  = fxs.fluctuation_scattering(dataset_name                     = exprun,
                                    detector_address                 = argv.address,
                                    data_type                        = argv.ftype,
                                    mask_path                        = argv.mask_path,
                                    mask_angles                      = None, #np.array([88, 270]),    # static masking at 88 and 270 deg
                                    mask_widths                      = None, #np.array([6,  10]),     # +/- degrees
                                    backimg_path                     = argv.bg_img_path,
                                    backmsk_path                     = argv.bg_msk_path,
                                    geom_path                        = argv.geom_path,
                                    det_dist                         = argv.det_distance,
                                    det_pix                          = argv.det_pixel,
                                    beam_l                           = argv.lambda_b,
                                    mask_thr                         = argv.thr,
                                    nQ                               = argv.nQ,
                                    nPhi                             = argv.nPhi,
                                    dQ                               = argv.dQ,
                                    dPhi                             = argv.dP,
                                    cent0                            = [argv.x,argv.y],
                                    r_max                            = argv.r_max,
                                    dr                               = argv.dr,
                                    dx                               = argv.dx,
                                    dy                               = argv.dy,
                                    r_0                              = argv.r0,
                                    q_bound                          = argv.q_bound)


  # Initialize iterator
  FXS.cnt       = np.array([0.])

  # Initialize Index variables
  if argv.param_path is None :
     maxevents = 400000          # We don't always know the total nr of events. Therefore set to large value
  else:
     maxevents = min(len(timestamps),len(timestamps[first:last]))


  FXS.get_index2(maxevents)
  # chop the list into pieces, depending on rank.  This assigns each process
  # events such that the get every Nth event where N is the number of processes

  if size > 1 :
     if rank > 0 :

        # MPI process. Here we set rank 0 to work as a listening server only.
        for j,evt in evtgen(run,timestamps = timestamps, first = first, last = last):
            #print '***',rank,j,evt.get(EventId).fiducials()
            if j%10==0: print('Rank',rank,'processing event',j)

            if ftype == 'h5' :
               FXS.get_h5(filestamps[j],evt)
            else :
               FXS.get_calib(evt)


            # Process hits
            if (FXS.image is not None) and (float(FXS.image.sum()) > hit)  :

               if ftype == 'h5' :
                  FXS.store_index_h5(evt, j)
               else:
                  ######################################
                  # Ugly way to get the time-stamps. Fix!!
                  time = evt.get(EventId).time()
                  fid = evt.get(EventId).fiducials()
                  sec  = time[0]
                  nsec = time[1]
                  et = EventTime(int((sec<<32)|nsec),fid)
                  #######################################
                  FXS.store_index2(et, j)                                         # Store index
               FXS.cnt  += 1


  else :


     # Single CPU
     for j,evt in evtgen(run,timestamps = timestamps, first = first, last = last):
         #print '***',rank,j,evt.get(EventId).fiducials()
         if j%10==0: print('Rank',rank,'processing event',j)


         if ftype == 'h5' :
            FXS.get_h5(filestamps[j],evt)
         else :
            FXS.get_calib(evt)

         # Process hits
         if (FXS.image is not None) and (float(FXS.image.sum()) > hit) :


            if ftype == 'h5' :
               FXS.store_index_h5(evt, j)
            else:
               ######################################
               # Ugly way to get the time-stamps. Fix!!
               time = evt.get(EventId).time()
               fid = evt.get(EventId).fiducials()
               sec  = time[0]
               nsec = time[1]
               et = EventTime(int((sec<<32)|nsec),fid)
               #######################################
               FXS.store_index2(et, j)                                         # Store index
            FXS.cnt  += 1

     print('Rank',rank,'total events:   ', int(FXS.cnt),' * ')



  #sum the images across mpi cores
  if size > 1:
    print("Synchronizing rank", rank)

  Tot         = np.zeros(FXS.cnt.shape)
  comm.Reduce(FXS.cnt,Tot)


  # Collect Indexing variables

  AVE_all           = np.zeros(FXS.ave.shape)
  comm.Reduce(FXS.ave,AVE_all)

  Tot_t       = np.zeros(FXS.tot_t.shape)
  comm.Reduce(FXS.tot_t,Tot_t)

  Tot_s       = np.zeros(FXS.tot_s.shape)
  comm.Reduce(FXS.tot_s,Tot_s)

  Tot_ns      = np.zeros(FXS.tot_ns.shape)
  comm.Reduce(FXS.tot_ns,Tot_ns)

  Tot_fd      = np.zeros(FXS.tot_fd.shape)
  comm.Reduce(FXS.tot_fd,Tot_fd)

  Tot_int     = np.zeros(FXS.tot_int.shape)
  comm.Reduce(FXS.tot_int,Tot_int)

  Tot_cx     = np.zeros(FXS.tot_cx.shape)
  comm.Reduce(FXS.tot_cx,Tot_cx)

  Tot_cy     = np.zeros(FXS.tot_cy.shape)
  comm.Reduce(FXS.tot_cy,Tot_cy)

  Tot_size   = np.zeros(FXS.tot_size.shape)
  comm.Reduce(FXS.tot_size,Tot_size)

  Tot_score  = np.zeros(FXS.tot_score.shape)
  comm.Reduce(FXS.tot_score,Tot_score)


  # Reduce results

  if rank==0:

    if size > 1:
      print("Synchronized")

    # Write out data

    if argv.outputdir is None:
        opath = os.getcwd()
    else:
        opath = argv.outputdir

    f_ave       = os.path.join(opath,'Average_run' + str(argv.run) + '_'+ str(int(Tot)) + '.dat')

    f_index     = os.path.join(opath,'Index_run' + str(argv.run) + '.dat')
    stamps      = ['Time','Seconds','Nanoseconds','Fiducial','Total Intensity','Beam X','Beam Y','Radius [Ang]','Score']
    head        ="                 ".join(stamps)


    Ave          = AVE_all / Tot


    # Get rid of zero lines add the end
    # Last non-zero intensity
    nz   = np.nonzero(Tot_t)
    fend = nz[0][-1]+1

    f              = open(f_index,'w')
    np.savetxt(f,np.c_[Tot_t[:fend],Tot_s[:fend],Tot_ns[:fend],Tot_fd[:fend],Tot_int[:fend],Tot_cx[:fend],Tot_cy[:fend],Tot_size[:fend],Tot_score[:fend]],header = head, comments='' )
    f.close()


    f              = open(f_ave,'w')
    np.savetxt(f_ave,Ave.reshape((-1)))
    f.close()
