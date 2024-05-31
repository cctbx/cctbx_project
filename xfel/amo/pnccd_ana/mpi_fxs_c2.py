from __future__ import absolute_import, division, print_function
from six.moves import range

from psana import *
import sys
import numpy as np
from xfel.amo.pnccd_ana             import pnccd_tbx
from xfel.amo.pnccd_ana             import pnccd_hit
from xfel.amo.pnccd_ana             import fxs
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



def compute_c2(argv=None) :

  """Function to compute the average angular intensity function, <C2(q,dPhi)>, from FXS images
       extracted from xtc (smd,idx,xtc format) or h5 files.
       Works for Single CPU, Multi-Processor interactive jobs and MPI batch jobs

       For a definition of input arguments argv and batch processing instructions see  ***  mpi_fxs_launch.py ***

       compute_c2 produces the following output files:

       * Index file : Information about the events processed including time-stamps, beam center, total and peak intensities, streak locations, particle size etc
       * I(q)       : Average Azimuthal intensities (SAXS)
       * C2(q,dPhi) : Average 2-point correlation functions (FXS)
       * Bg_img     : Average image in polar coordinates
       * Bg_norm    : Mean subtracted average image in polar coordinates
       * Bg_msk     : Average mask image in polar coordinates

       The I and C2 output data are split in 2 ~equal bins [0 and 1] for statistical comparison to be made

       example:
                C2_run111_0_29793.dat  is the average C2(q,dPhi) from run111 bin 0 where 29,793 images have been processed
                C2_run111_1_29718.dat  is the average C2(q,dPhi) from run111 bin 1 where 29,718 images have been processed

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
       elif (ftype == 'idx') :
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
                                    mask_angles                      = None,#np.array([88, 270]),    # static masking at 88 and 270 deg
                                    mask_widths                      = None,#np.array([6,  10]),     # +/- degrees
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
                                    q_bound                          = argv.q_bound,
                                    peak                             = [0.037, 0.064],              # Location in q for peaks to be integrated
                                    dpeak                            = [0.002, 0.002])              # Width    in q for peaks to be integrated


  # Initialize iterator
  FXS.cnt_0 = np.array([0.])
  FXS.cnt_1 = np.array([0.])

  # Initialize Index variables
  if argv.param_path is None :
     maxevents = 400000          # We don't always know the total nr of events. Therefore set to large value
  else:
     maxevents = min(len(timestamps),len(timestamps[first:last]))


  FXS.get_index(maxevents)
  # chop the list into pieces, depending on rank.  This assigns each process
  # events such that the get every Nth event where N is the number of processes

  if size > 1 :
     if rank > 0 :

        hd=pnccd_hit.hit()

        # MPI process. Here we set rank 0 to work as a listening server only.
        for j,evt in evtgen(run,timestamps = timestamps, first = first, last = last):
            #print '***',rank,j,evt.get(EventId).fiducials()
            if j%10==0: print('Rank',rank,'processing event',j)

            if ftype == 'h5' :
               FXS.get_h5(filestamps[j],evt)
            else :
               FXS.get_image(evt)

            # Process hits

            if (FXS.img is not None) and (float(FXS.img.sum()) > hit)  :

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

               FXS.sum_bg2()                                                    # Sum Background


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
                  FXS.store_index(et, j)                                       # Store index

               if int(FXS.cnt_0 + FXS.cnt_1)%10==0: print('Rank',rank,'processed events: ', int(FXS.cnt_0 + FXS.cnt_1))


               # Send partial results to master (rank 0)
               if int(FXS.cnt_0 + FXS.cnt_1) % 50 == 0:                        # Send every 50 events


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

        print('Rank',rank,'total events:     ', int(FXS.cnt_0 + FXS.cnt_1),' * ')

     else:

        if ftype == 'h5' :
           FXS.run_nr      = run
        else:
           FXS.run_nr      = int(run.run())


        hd              = pnccd_hit.hit()

        adim            = FXS.ave.shape
        sdim            = len(FXS.q)
        cdim            = (len(FXS.q),len(FXS.phi))
        idim            = (maxevents,3)


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


     for j,evt in evtgen(run,timestamps = timestamps, first = first, last = last):
         #print '***',rank,j,evt.get(EventId).fiducials()
         if j%10==0: print('Rank',rank,'processing event',j)


         if ftype == 'h5' :
            FXS.get_h5(filestamps[j],evt)
         else :
            FXS.get_image(evt)


         # Process hits
         if (FXS.img is not None) and (float(FXS.img.sum()) > hit) :

             FXS.get_beam(plot=plot)                                      # Beam center refinement
             FXS.get_polar()                                              # Polar transform
             FXS.get_streak_mask(plot=plot)                               # Mask out streaks
             FXS.get_pixel_mask(plot=plot)                                # Mask out pixels
             FXS.get_norm(plot=plot)                                      # Normalize image, get SAXS
             FXS.get_c2(plot=plot)                                        # Compute C2


             # Split upp into 2 half sets
             if int(FXS.cnt_0 + FXS.cnt_1) % 2 == 0:
                FXS.sum_c2(flag = 0)                                     # Accumulate C2 for Half I
             else:
                FXS.sum_c2(flag = 1)                                     # Accumulate C2 for Half II

             if FXS.r_0 is not None :
                FXS.get_size()


             FXS.sum_bg2()                                              # Sum Background

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
                FXS.store_index(et, j)                                  # Store index

     print('Rank',rank,'total events:   ', int(FXS.cnt_0 + FXS.cnt_1),' * ')


  #sum the images across mpi cores
  if size > 1:
    print("Synchronizing rank", rank)

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

  if not hasattr(FXS, 'Back_img'):
     FXS.Back_img   = np.zeros((len(FXS.q),len(FXS.phi)))
  if not hasattr(FXS, 'Back_norm'):
     FXS.Back_norm  = np.zeros((len(FXS.q),len(FXS.phi)))
  if not hasattr(FXS, 'Back_msk'):
     FXS.Back_msk   = np.zeros((len(FXS.q),len(FXS.phi)))



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

  Tot_peak1   = np.zeros(FXS.tot_peak1_int.shape)
  comm.Reduce(FXS.tot_peak1_int,Tot_peak1)

  Tot_peak2   = np.zeros(FXS.tot_peak2_int.shape)
  comm.Reduce(FXS.tot_peak2_int,Tot_peak2)

  Tot_s_m    = np.zeros(FXS.tot_streak_m.shape)
  comm.Reduce(FXS.tot_streak_m,Tot_s_m)

  Tot_s_s    = np.zeros(FXS.tot_streak_s.shape)
  comm.Reduce(FXS.tot_streak_s,Tot_s_s)

  Tot_cx     = np.zeros(FXS.tot_cx.shape)
  comm.Reduce(FXS.tot_cx,Tot_cx)

  Tot_cy     = np.zeros(FXS.tot_cy.shape)
  comm.Reduce(FXS.tot_cy,Tot_cy)

  Tot_size   = np.zeros(FXS.tot_size.shape)
  comm.Reduce(FXS.tot_size,Tot_size)

  Tot_score  = np.zeros(FXS.tot_score.shape)
  comm.Reduce(FXS.tot_score,Tot_score)


  # Collect background variables

  BG_img_all        = np.zeros(FXS.Back_img.shape)
  comm.Reduce(FXS.Back_img,BG_img_all)

  BG_norm_all       = np.zeros(FXS.Back_norm.shape)
  comm.Reduce(FXS.Back_norm,BG_norm_all)

  BG_msk_all        = np.zeros(FXS.Back_msk.shape)
  comm.Reduce(FXS.Back_msk,BG_msk_all)


  # Reduce results

  if rank==0:

    if size > 1:
      print("Synchronized")

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
    stamps      = ['Time','Seconds','Nanoseconds','Fiducial','Total Intensity','Peak1, q='+str(FXS.peak[0])+'+/-'+str(FXS.dpeak[0]),'Peak2, q='+str(FXS.peak[1])+'+/-'+str(FXS.dpeak[1]),'Mean streak angle','Std streak angle','Beam X','Beam Y','Radius [Ang]','Score']
    head        ="                 ".join(stamps)



    Tot_0          = int(Tot_0)
    Isaxs_ave_0    = SAXS_0_all / Tot_0
    Isaxs_std_0    = np.sqrt( VAR_0_all / Tot_0 )
    C2_ave_0       = ( C2_0_all / Tot_0 )   / ( C2m_0_all / Tot_0 )

    Tot_1          = int(Tot_1)
    Isaxs_ave_1    = SAXS_1_all / Tot_1
    Isaxs_std_1    = np.sqrt( VAR_1_all / Tot_1 )
    C2_ave_1       = ( C2_1_all / Tot_1 )   / ( C2m_1_all / Tot_1 )

    # Prepare background
    tmp2           = np.copy(BG_msk_all)
    ind            = tmp2 == 0
    tmp2[ind]      = 1.0
    Bg_img         = BG_img_all / tmp2
    Bg_norm        = BG_norm_all / tmp2

    Bg_msk         = np.ones(tmp2.shape)
    Bg_msk[ind]    = 0.0


    f              = open(f_saxs0,'w')
    np.savetxt(f,np.c_[FXS.q,Isaxs_ave_0,Isaxs_std_0],header = head_s, comments='')
    f.close()

    f              = open(f_c0,'w')
    np.savetxt(f,C2_ave_0,header = head_c, comments='')
    f.close()


    f              = open(f_saxs1,'w')
    np.savetxt(f,np.c_[FXS.q,Isaxs_ave_1,Isaxs_std_1],header = head_s, comments='')
    f.close()

    f              = open(f_c1,'w')
    np.savetxt(f,C2_ave_1,header = head_c, comments='')
    f.close()

    f_bg_im     = os.path.join(opath,'Bg_img_' + str(argv.run) + '_'+ str(Tot_0+Tot_1) + '.dat')
    f_bg_norm   = os.path.join(opath,'Bg_norm_' + str(argv.run) + '_'+ str(Tot_0+Tot_1) + '.dat')
    f_bg_ms     = os.path.join(opath,'Bg_msk_' + str(argv.run) + '_'+ str(Tot_0+Tot_1) + '.dat')
    stamps_s    = ['q','Mean','Std']
    head_s      ="                 ".join(stamps_s)


    # Get rid of zero lines at the end
    # Last non-zero intensity
    nz   = np.nonzero(Tot_t)
    fend = nz[0][-1]+1

    f              = open(f_index,'w')
    np.savetxt(f,np.c_[Tot_t[:fend],Tot_s[:fend],Tot_ns[:fend],Tot_fd[:fend],Tot_int[:fend],Tot_peak1[:fend],Tot_peak2[:fend],Tot_s_m[:fend],Tot_s_s[:fend],Tot_cx[:fend],Tot_cy[:fend],Tot_size[:fend],Tot_score[:fend]],header = head, comments='' )
    f.close()


    f              = open(f_bg_im,'w')
    np.savetxt(f,Bg_img)
    f.close()

    f              = open(f_bg_norm,'w')
    np.savetxt(f,Bg_norm)
    f.close()

    f           = open(f_bg_ms,'w')
    np.savetxt(f,Bg_msk)
    f.close()
