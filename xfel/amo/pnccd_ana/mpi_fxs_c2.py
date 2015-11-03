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


def compute_c2(argv=None) :
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
  FXS.cnt_0 = np.array([0.])
  FXS.cnt_1 = np.array([0.])

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
               FXS.get_c2(plot = plot)                                          # Compute C2

               # Split upp into 2 half sets
               if int(FXS.cnt_0 + FXS.cnt_1) % 2 == 0:
                  FXS.sum_c2(flag = 0)                                          # Accumulate C2 for Half I
               else:
                  FXS.sum_c2(flag = 1)                                          # Accumulate C2 for Half II

               if FXS.r_0 is not None :
                  FXS.get_size()

               FXS.store_index(mytimes[j], myevents[j])



               # Send partial results to master (rank 0)
               if int(FXS.cnt_0 + FXS.cnt_1) % 100 == 0:                        # Send every 100 events


                  # C2 and Saxs data
                  tmp_n    = int(FXS.cnt_0   +  FXS.cnt_1)
                  tmp_saxs = (FXS.Isaxs_0 +  FXS.Isaxs_1) / tmp_n
                  tmp_c2   = (FXS.C2_0    +  FXS.C2_1) / (FXS.C2m_0   +  FXS.C2m_1)


                  # Average image
                  tmp_im   = FXS.ave / tmp_n

                  # Ascert no nan values
                  tmp_saxs[np.isnan(tmp_saxs)] = 0
                  tmp_c2[np.isnan(tmp_c2)]     = 0

                  # Total intensity, Size and Score

                  tmp_ind = np.column_stack((FXS.tot_int,FXS.tot_size,FXS.tot_score))

                  hd.send(tmp_n,image = tmp_im, saxs=tmp_saxs,c2=tmp_c2,ind=tmp_ind)

        hd.endrun()

     else:

        FXS.run_nr      = int(run.run())


        hd              = pnccd_hit.hit()

        adim            = FXS.ave.shape
        sdim            = len(FXS.q)
        cdim            = (len(FXS.q),len(FXS.phi))
        idim            = (nevents,3)


        hd.total_ave    = [np.zeros(adim)]*(size-1)
        hd.total_c2     = [np.zeros(cdim)]*(size-1)
        hd.total_ind    = [np.zeros(idim)]*(size-1)
        hd.total_saxs   = [np.zeros(sdim)]*(size-1)
        hd.total_ev_a   = [0.0]*(size-1)
        hd.total_ev_c   = [0.0]*(size-1)
        hd.total_ev_s   = [0.0]*(size-1)
        hd.total_ev_i   = [0.0]*(size-1)

        nClients = size - 1

        while nClients > 0:
            # Remove client if the run ended
            if hd.recv():
               nClients -= 1
            else:
               na = sum(hd.total_ev_a)
               ns = sum(hd.total_ev_s)
               nc = sum(hd.total_ev_c)
               ni = sum(hd.total_ev_i)

               if (na == ns ==  nc == ni) and (ns % 100 == 0) : # Publish every 100 events

                  AVE     = np.zeros(adim)
                  C2ave   = np.zeros(cdim)
                  SAXSave = np.zeros(sdim)
                  IND     = np.zeros(idim)

                  for i in range(size-1) :
                      AVE     = AVE     + (hd.total_ave[i] * (hd.total_ev_a[i] /na))
                      C2ave   = C2ave   + (hd.total_c2[i] * (hd.total_ev_c[i] /nc))
                      SAXSave = SAXSave + (hd.total_saxs[i] * (hd.total_ev_s[i] /ns))
                      IND     = IND     + hd.total_ind[i]

                  FXS.publish(image = AVE, saxs=SAXSave, c2=C2ave, ind=IND, n_a=na, n_saxs=ns, n_c2=nc, n_i=ni)


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
             FXS.get_c2(plot = plot)                                        # Compute C2

             # Split upp into 2 half sets
             if int(FXS.cnt_0 + FXS.cnt_1) % 2 == 0:
                FXS.sum_c2(flag = 0)                                        # Accumulate C2 for Half I
             else:
                FXS.sum_c2(flag = 1)                                        # Accumulate C2 for Half II

             if FXS.r_0 is not None :
                FXS.get_size()

             FXS.store_index(mytimes[j], myevents[j])

  #sum the images across mpi cores
  if size > 1:
    print "Synchronizing rank", rank

  Tot_0          = np.zeros(FXS.cnt_0.shape)
  Tot_1          = np.zeros(FXS.cnt_1.shape)

  comm.Reduce(FXS.cnt_0,Tot_0)
  comm.Reduce(FXS.cnt_1,Tot_1)


  if rank == 0 and Tot_0[0] == 0 and Tot_1[0]:
    raise Sorry("No events found in the run")


  if not hasattr(FXS, 'Isaxs_0'):
     FXS.Isaxs_0    = np.zeros(FXS.q.shape)
  if not hasattr(FXS, 'Vsaxs_0'):
     FXS.Vsaxs_0    = np.zeros(FXS.q.shape)
  if not hasattr(FXS, 'C2_0'):
     FXS.C2_0       = np.zeros((len(FXS.q),len(FXS.phi)))
  if not hasattr(FXS, 'C2m_0'):
     FXS.C2m_0      = np.zeros((len(FXS.q),len(FXS.phi)))

  if not hasattr(FXS, 'Isaxs_1'):
     FXS.Isaxs_1    = np.zeros(FXS.q.shape)
  if not hasattr(FXS, 'Vsaxs_1'):
     FXS.Vsaxs_1    = np.zeros(FXS.q.shape)
  if not hasattr(FXS, 'C2_1'):
     FXS.C2_1       = np.zeros((len(FXS.q),len(FXS.phi)))
  if not hasattr(FXS, 'C2m_1'):
     FXS.C2m_1      = np.zeros((len(FXS.q),len(FXS.phi)))


  # Collect  Variables

  SAXS_0_all     = np.zeros(FXS.Isaxs_0.shape)
  comm.Reduce(FXS.Isaxs_0,SAXS_0_all)

  VAR_0_all      = np.zeros(FXS.Vsaxs_0.shape)
  comm.Reduce(FXS.Vsaxs_0,VAR_0_all)

  C2_0_all       = np.zeros(FXS.C2_0.shape)
  comm.Reduce(FXS.C2_0,C2_0_all)

  C2m_0_all      = np.zeros(FXS.C2m_0.shape)
  comm.Reduce(FXS.C2m_0,C2m_0_all)

  SAXS_1_all     = np.zeros(FXS.Isaxs_1.shape)
  comm.Reduce(FXS.Isaxs_1,SAXS_1_all)

  VAR_1_all     = np.zeros(FXS.Vsaxs_1.shape)
  comm.Reduce(FXS.Vsaxs_1,VAR_1_all)

  C2_1_all       = np.zeros(FXS.C2_1.shape)
  comm.Reduce(FXS.C2_1,C2_1_all)

  C2m_1_all      = np.zeros(FXS.C2m_1.shape)
  comm.Reduce(FXS.C2m_1,C2m_1_all)


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

    f_saxs0     = os.path.join(opath,'Saxs_run' + str(argv.run) + '_0_'+ str(int(Tot_0)) + '.dat')
    f_saxs1     = os.path.join(opath,'Saxs_run' + str(argv.run) + '_1_'+ str(int(Tot_1)) + '.dat')
    stamps_s    = ['q','Mean','Std']
    head_s      ="                 ".join(stamps_s)

    f_c0        = os.path.join(opath,'C2_run' + str(argv.run) + '_0_'+ str(int(Tot_0)) + '.dat')
    f_c1        = os.path.join(opath,'C2_run' + str(argv.run) + '_1_'+ str(int(Tot_1)) + '.dat')
    stamps_c    = ['C2']
    head_c      ="                 ".join(stamps_c)

    f_index     = os.path.join(opath,'Index_run' + str(argv.run) + '.dat')
    stamps      = ['Time','Seconds','Nanoseconds','Fiducial','Total Intensity','Beam X','Beam Y','Radius [Ang]','Score']
    head        ="                 ".join(stamps)


    Tot_0          = int(Tot_0)
    Isaxs_ave_0    = SAXS_0_all / Tot_0
    Isaxs_std_0    = np.sqrt( VAR_0_all / Tot_0 )
    C2_ave_0       = ( C2_0_all / Tot_0 ) / ( C2m_0_all / Tot_0 )

    f              = open(f_saxs0,'w')
    np.savetxt(f,np.c_[FXS.q,Isaxs_ave_0,Isaxs_std_0],header = head_s, comments='')
    f.close()

    f              = open(f_c0,'w')
    np.savetxt(f,C2_ave_0,header = head_c, comments='')
    f.close()

    Tot_1          = int(Tot_1)
    Isaxs_ave_1    = SAXS_1_all / Tot_1
    Isaxs_std_1    = np.sqrt( VAR_1_all / Tot_1 )
    C2_ave_1       = ( C2_1_all / Tot_1 ) / ( C2m_1_all / Tot_1 )

    f              = open(f_saxs1,'w')
    np.savetxt(f,np.c_[FXS.q,Isaxs_ave_1,Isaxs_std_1],header = head_s, comments='')
    f.close()

    f              = open(f_c1,'w')
    np.savetxt(f,C2_ave_1,header = head_c, comments='')
    f.close()

    f              = open(f_index,'w')
    np.savetxt(f,np.c_[Tot_t,Tot_s,Tot_ns,Tot_fd,Tot_int,Tot_cx,Tot_cy,Tot_size,Tot_score],header = head, comments='' )
    f.close()
