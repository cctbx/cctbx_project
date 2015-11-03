from __future__ import division

from psana import *
import sys
import numpy as np
from xfel.amo.pnccd_ana                 import pnccd_tbx
from xfel.amo.pnccd_ana                 import pnccd_hit
from xfel.amo.pnccd_ana                 import fxs
import matplotlib.pyplot as plt

plt.ion()
########################################
# Due to the mask sometimes having zero values
# we're bound to get divisions with zeros at
#times. Here ignoring those errors.
np.seterr(divide='ignore', invalid='ignore')


def compute_index(argv=None) :
  if argv == None:
    argv = sys.argv[1:]

  try:
     from mpi4py import MPI
  except ImportError:
     raise Sorry("MPI not found")

  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  size = comm.Get_size()

  if argv.hit is None :
     hit        = -1.0e6                        # Process everything
  else:
     hit        = argv.hit      # Process everything > hit


  if argv.param_path is not None :
     param_file = np.genfromtxt(argv.param_path,skiprows=1)


  if (argv.ftype == 'xtc') or (argv.ftype == 'ffb') or (argv.ftype == 'smd')  :
       if  (argv.ftype == 'xtc') :
           dataset_name = "exp=%s:run=%d:idx"%(argv.experiment, argv.run)
       elif(argv.ftype == 'ffb') :
           dataset_name = "exp=%s:run=%d:idx"%(argv.experiment, argv.run)
           # as ffb is only at SLAC, ok to hardcode /reg/d here
           dataset_name += ":dir=/reg/d/ffb/%s/%s/xtc"%(argv.experiment[0:3],argv.experiment)
       elif(argv.ftype == 'smd') :
           dataset_name = "exp=%s:run=%d:smd"%(argv.experiment, argv.run)
           # as ffb is only at SLAC, ok to hardcode /reg/d here ADD live!
           dataset_name += ":dir=/reg/d/ffb/%s/%s/xtc:live"%(argv.experiment[0:3],argv.experiment)


       ds           = DataSource(dataset_name)

       # Get run
       for run in ds.runs():
           if rank == 0:
              print "Processing run ", run.run()

       # Get timestamps
       if argv.param_path is not None :
          timestamps      = pnccd_tbx.get_psana_event(param_file)
       else:
          timestamps      = run.times()

  elif argv.ftype == 'h5' :
       import h5py
       dataset_name = "%s_run_%d.h5"%(argv.experiment, argv.run)
       dataset_name = os.path.join(argv.xtc_dir,dataset_name)
       run          = int(argv.run)
       f               = h5py.File(dataset_name,'r')
       timestamps      = f.keys()

       if rank == 0:
          print "Processing run ", argv.run


  if argv.first is None :
     first   = 0
  else:
     first   = argv.first
  if argv.last is None :
     last    = len(timestamps)
  else:
     last    = min(argv.last,len(timestamps))


  times     = timestamps[first:last]
  nevents   = len(times)

  if rank == 0 :
     print "Processing events " +str(first)+ " to " +str(last)


  if size == 1:
     plot = argv.plot
  else:
     plot = 0


  FXS  = fxs.fluctuation_scattering(dataset_name                     = dataset_name,
                                    detector_address                 = argv.address,
                                    data_type                        = argv.ftype,
                                    mask_path                        = argv.mask_path,
                                    mask_angles                      = None,#np.array([90, 270])    # static masking at 90 and 270
                                    mask_widths                      = None,#np.array([10,  10])    # +/- degrees
                                    backimg_path                     = argv.bg_img_path,
                                    backmsk_path                     = argv.bg_msk_path,
                                    param_path                       = argv.param_path,
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
  FXS.get_index(nevents)

  # chop the list into pieces, depending on rank.  This assigns each process
  # events such that the get every Nth event where N is the number of processes

  if size > 1 :
     if rank > 0 :

        hd=pnccd_hit.hit()

        # MPI process. Here we set rank 0 to work as a listening server only.
        mytimes,myevents  = zip(*[(times[i],i) for i in xrange(nevents) if (i+rank)%(size-1) == 0])

        for j in xrange(len(mytimes)):
            if j%10==0: print 'Rank',rank,'processing event',rank*len(mytimes)+j,', ',j,'of',len(mytimes)


            FXS.get_image(run,mytimes[j])

            # Process hits
            if float(FXS.img.sum()) >= hit :

               FXS.get_beam(plot = plot)                                        # Beam center refinement
               FXS.get_polar(plot = plot)                                       # Polar transform
               FXS.get_streak_mask(plot = plot)                                 # Mask out streaks
               FXS.get_pixel_mask(plot = plot)                                  # Mask out pixels
               FXS.get_norm(plot = plot)                                        # Normalize image, get SAXS


               if FXS.r_0 is not None :
                  FXS.get_size()

               FXS.store_index(mytimes[j], myevents[j])


               # Send partial results to master (rank 0)
               if (int(FXS.cnt) > 0) and (int(FXS.cnt) % 100 == 0):             # Send every 100 events

                  tmp_n    = int(FXS.cnt)

                  # Average image
                  tmp_im   = FXS.ave / tmp_n

                  # Total intensity, Size and Score
                  tmp_ind = np.column_stack((FXS.tot_int,FXS.tot_size,FXS.tot_score))

                  hd.send(tmp_n, image = tmp_im, ind=tmp_ind)

            FXS.cnt  += 1


        hd.endrun()

     else:

        FXS.run_nr      = int(run.run())


        hd              = pnccd_hit.hit()
        adim            = FXS.ave.shape
        idim            = (nevents,3)

        hd.total_ave    = [np.zeros(adim)]*(size-1)
        hd.total_ind    = [np.zeros(idim)]*(size-1)
        hd.total_ev_a   = [0.0]*(size-1)
        hd.total_ev_i   = [0.0]*(size-1)

        nClients = size - 1

        while nClients > 0:
            # Remove client if the run ended
            if hd.recv():
               nClients -= 1
            else:
               na = sum(hd.total_ev_a)
               ni = sum(hd.total_ev_i)

               if  (na == ni) and  (na % 100 == 0) :                                            # Publish every 100 events


                  AVE     = np.zeros(adim)
                  IND     = np.zeros(idim)

                  for i in range(size-1) :
                      AVE     = AVE     + (hd.total_ave[i] * (hd.total_ev_a[i] /na))
                      IND     = IND     + hd.total_ind[i]

                  FXS.publish(image = AVE, ind=IND, n_a=na, n_i=ni)


  else :


     # Single CPU
     mytimes,myevents  = zip(*[(times[i],i) for i in xrange(nevents) if (i+rank)%size == 0])

     for j in xrange(len(mytimes)):
         if j%10==0: print 'Rank',rank,'processing event',rank*len(mytimes)+j,', ',j,'of',len(mytimes)


         FXS.get_image(run,mytimes[j])

         # Process hits
         if float(FXS.img.sum()) >= hit :

             FXS.get_beam(plot = plot)                                      # Beam center refinement
             FXS.get_polar(plot = plot)                                     # Polar transform
             FXS.get_streak_mask(plot = plot)                               # Mask out streaks
             FXS.get_pixel_mask(plot = plot)                                # Mask out pixels
             FXS.get_norm(plot = plot)                                      # Normalize image, get SAXS

             if FXS.r_0 is not None :
                FXS.get_size()

             FXS.store_index(mytimes[j], myevents[j])

             FXS.cnt  += 1


  #sum the images across mpi cores
  if size > 1:
    print "Synchronizing rank", rank

  Tot         = np.zeros(FXS.cnt.shape)
  comm.Reduce(FXS.cnt,Tot)



  if rank == 0 and Tot[0] == 0 :
    raise Sorry("No events found in the run")


  # Collect Indexing variables

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
      print "Synchronized"

    # Write out data

    if argv.outputdir is None:
        opath = os.getcwd()
    else:
        opath = argv.outputdir

    f_index     = os.path.join(opath,'Index_run' + str(argv.run) + '.dat')
    stamps      = ['Time','Seconds','Nanoseconds','Fiducial','Total Intensity','Beam X','Beam Y','Radius [Ang]','Score']
    head        ="                 ".join(stamps)


    f              = open(f_index,'w')
    np.savetxt(f,np.c_[Tot_t,Tot_s,Tot_ns,Tot_fd,Tot_int,Tot_cx,Tot_cy,Tot_size,Tot_score],header = head, comments='' )
    f.close()
