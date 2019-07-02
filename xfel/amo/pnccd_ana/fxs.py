
"""Class for processing Fluctuation X-ray Scattering data (FXS)

"""
from __future__ import absolute_import, division, print_function
from six.moves import range
from psana                              import *

import numpy             as np
import matplotlib.pyplot as plt
import time
import h5py

from xfel.cxi.cspad_ana                 import cspad_tbx
from xfel.amo.pnccd_ana                 import pnccd_tbx
from psmon                              import publish
from psmon.plots                        import Image, XYPlot, Hist, MultiPlot


class fluctuation_scattering(object):
  """Class for processing of 2D fluctuation scattering images.

  """


  def __init__(self,
               dataset_name     = None,
               detector_address = None,
               data_type        = 'idx',
               mask_path        = None,
               mask_angles      = None,
               mask_widths      = None,
               backimg_path     = None,
               backmsk_path     = None,
               geom_path        = None,
               det_dist         = None,
               det_pix          = 0.075,
               beam_l           = None,
               mask_thr         = None,
               nQ               = None,
               nPhi             = None,
               dQ               = 1,
               dPhi             = 1,
               cent0            = None,
               r_max            = None,
               dr               = None,
               dx               = None,
               dy               = None,
               r_0              = None,
               q_bound          = None,
               peak             = None,
               dpeak            = None):

    """The fluctuation scattering class stores processing parameters,
       initiates mask and background data and retrieves 2D images
       from events. Processing options of the 2D images include:
       * Transform from cartesian to polar coordinates
       * Beam center refinement
       * Dynamic masking
       * Normalization % SAXS calculation
       * Particle sizing
       * Computation of in-frame 2-point angular auto-correlations using FFTs

    @param dataset_name         Experiment name and run number
    @param detector_address     Adress to back or front detector
    @param data_type            Type of data file format (h5 or xtc in the formats: idx|idx_ffb|smd|smd_ffb|h5)
    @param mask_path            Full path to static image mask
    @param mask_angles          Center of angluar slices (deg) that should be masked out (due to jet streaks etc), [Ang1 Ang2 ...]
    @param mask_widths          Width  of angular slices (deg) that should be masked out (due to jet streaks etc), [delta1 delta2 ...]
    @param backimg_path         Full path to background image
    @param backmsk_path         Full path to background mask
    @param geom_path            Full path to geometry file (for h5 format)
    @param det_dist             Override of detecor distance (in mm)
    @param det_pix              Pixel size (in mm)
    @param beam_l               Override of beam wavelength (in Angstrom)
    @param mask_thr             Threshold for dynamic masking
    @param nQ                   Number of Q-bins to consider   (in pixels)
    @param nPhi                 Number of Phi-bins to consider (in pixels)
    @param dQ                   Stepsize in Q   (in pixels)
    @param dPhi                 Stepsize in Phi (in pixels)
    @param cent0                Initial beam center coordinates [xc,yc]
    @param r_max                Maximum radial value to use for beamcenter refinement (in pixels)
    @param dr                   Stepsize in r (in pixels)
    @param dx                   Gridsize for beam center refinement in x, i.e xc+/-dx (in pixels)
    @param dy                   Gridsize for beam center refinement in y, i.e yc+/-dy (in pixles)
    @param r_0                  Starting value for particle radius refinement [in Ang]
    @param q_bound              Upper and Lower boundaries of q for Particle radius refinement [in Ang^-1]
    @param peak                 Q-values for peak maxima [q_peak1 q_peak2 ...]
    @param dpeak                Delta Q used for peak integration of peak maxima [delta_q1 delta_q2 ...]
    """


    # Initialize parameters and configuration files once

    self.data_type          = data_type
    self.dataset_name       = dataset_name
    self.detector_address   = detector_address


    if (self.data_type == 'idx') or (self.data_type == 'idx_ffb') or (self.data_type == 'smd')  or (self.data_type == 'smd_ffb') or (self.data_type == 'xtc') :


       self.ds       = DataSource(self.dataset_name)
       self.src      = Detector(self.detector_address, self.ds.env())


       if mask_path is None :                   # Create a binary mask of ones, default mask only works for xtc/ffb

          evt               = next(self.ds.events())
          self.mask_address = self.src.mask(evt,calib=True,status=True)
          self.msk          = self.src.image(evt,self.mask_address)
          self.mask         = np.copy(self.msk)

       else:

          self.msk          = np.loadtxt(mask_path)
          self.mask         = np.copy(self.msk)


    if geom_path is not None :
       geom          = np.genfromtxt(geom_path,skiprows=1)
       self.gap      = geom[0]
       self.shift    = geom[1]
       self.orient   = geom[2
]
       self.comm     = geom[3]
       self.param1   = geom[4]
       self.param2   = geom[5]
       self.param3   = geom[6]
       self.param4   = geom[7]

    if (self.data_type == 'h5'):

       if mask_path is None :                      # Create a binary mask of ones

         ## Add default binary mask here ## Default pnCCD dimensions
         dim1      = 1024
         dim2      = 1024
         self.mask = np.ones((dim1,dim2))

       else:

          self.mask = np.loadtxt(mask_path)

       # Apply geometry
       self.msk   = pnccd_tbx.get_geometry(img    = self.mask,
                                          gap    = self.gap,
                                          shift  = self.shift,
                                          orient = self.orient)


    if self.detector_address == 'pnccdFront' :

       evt          = next(self.ds.events())
       gain         = self.src.gain(evt)
       self.gain    = self.src.image(evt,gain)

    self.cart = 0
    self.flat = 0

    if backimg_path is None :
       self.backimg     = None
    else :
       self.backimg     = np.loadtxt(backimg_path).astype(np.float64)

       # Check if background image in cartesian coordinates exists
       if (self.backimg.shape == self.msk.shape):
          self.cart     = 1                             # Remove bg before transform to polar coordinates


    if backmsk_path is None :
       self.backmsk     = None
    else :
       self.backmsk     = np.loadtxt(backmsk_path).astype(np.float64)

    # Check if flat-field image exists
    if (self.backmsk is None) and (self.backimg is not None) and (self.cart==0):
       self.pcflat      = pnccd_tbx.dynamic_flatfield(self.backimg)
       self.flat        = 1


    if det_dist is None :                       # Get detector distance from events
       for run in self.ds.runs():
           self.det_dist = cspad_tbx.env_distance(self.detector_address, run.env(), 577)
    else :
       self.det_dist = det_dist

    self.det_pix = det_pix

    if beam_l is None :                        # Get wavelength from event, note it can change slightly between events. So in the future use average.
       self.beam_l   = cspad_tbx.evt_wavelength(next(self.ds.events()))
    else :
       self.beam_l   = beam_l

    if mask_thr is None :                      # No dynamic masking
       self.thr      = None
    else :
       self.thr      = mask_thr

    if nQ is None :                            # Use image dimensions as a guide, leave room for offset beamC
       if self.msk.shape[0] > self.msk.shape[1] :  # nQ determined by smallest dimension
          self.nQ   = int(self.msk.shape[1]/2)-20
       else :
          self.nQ   = int(self.msk.shape[0]/2)-20
    else :
       self.nQ       = nQ


    if (self.nQ % 10):                        # Ascert even number, speeds things up massively for FFT
        self.nQ  = np.floor(self.nQ/10)*10

    if (self.nQ % dQ):                        # Ascert clean divisor
        self.nQ  = np.floor(self.nQ/dQ)*dQ

    if nPhi is None :                         # Estimate based on 2*pi*nQ
       self.nPhi     = np.ceil(2*np.pi*self.nQ)
    else :
       self.nPhi     = nPhi

    if (self.nPhi % 10):                      # Ascert even number, speeds things up massively for FFT
        self.nPhi  = np.ceil(self.nPhi/10)*10

    if (self.nPhi % dPhi):                    # Ascert clean divisor
        self.nPhi  = np.ceil(self.nPhi/dPhi)*dPhi

    self.dQ          = dQ
    self.dPhi        = dPhi

    self.mask_angles        = mask_angles
    self.mask_widths        = mask_widths

    # Compute slices that should be masked in static mask
    if (self.mask_angles is not None) and (self.mask_widths is not None) :
       self.mask_angles        = (self.mask_angles/360) * self.nPhi
       self.mask_widths        = (self.mask_widths/360) * self.nPhi

    # Check if nQ > smallest dimension/2 then we extend the image with zero values
    if self.nQ > min(self.msk.shape[0]/2,self.msk.shape[1]/2) :
       self.msk   = pnccd_tbx.extend_image(img = self.msk)
       self.mask  = np.copy(self.msk)


    if (cent0 is None) or (sum(cent0) == 0):  # Use center of gravity to estimate starting beamC
       self.cent0    = [int(round(self.msk.shape[1]/2)) , int(round(self.msk.shape[0]/2))]
    else :
       self.cent0    = cent0

    self.cent = self.cent0                    # Default center

    if r_max is None :                        # Default, Use half of nQ
       self.r_max    = int(self.nQ*(3/4))
    else :
       self.r_max    = r_max

    if (self.r_max % dr):                     # Ascert clean divisor
        self.r_max  = np.floor(self.r_max/dr)*dr

    self.dr          = dr
    self.dx          = dx
    self.dy          = dy

    if r_0 is None :
       self.radius   = 0
       self.score    = 0

    self.r_0       = r_0

    if q_bound is None or sum(q_bound)==0 :
       self.q_bound    = [None,None]
    else :
       self.q_bound    = [None,self.q_bound]


    # Compute q-spacing
    self.q           = np.arange(0, self.nQ, self.dQ)
    self.q           = self.q*self.det_pix/self.det_dist*4*np.pi/self.beam_l/2

    # Compute Phi (Not accounting for curvature)
    self.phi         = np.linspace(0, 2*np.pi, self.nPhi/self.dPhi,endpoint=False)


    # Compute indices for Peak maxima
    if (peak is not None) and (dpeak is not None) :
       self.peak       = peak
       self.dpeak      = dpeak
       self.ind1       = (self.q >= (self.peak[0] - self.dpeak[0])) & (self.q <= (self.peak[0] + self.dpeak[0]) )
       self.ind2       = (self.q >= (self.peak[1] - self.dpeak[1])) & (self.q <= (self.peak[1] + self.dpeak[1]) )
    else:
       self.peak       = None
       self.dpeak      = None

    ###################################################################################
    # Define functions



  def publish(self, image = None, saxs = None, c2 = None, ind = None, n_a = None, n_saxs = None, n_c2 = None, n_i = None, n_q = None, n_bin = None) :
      """Publish Intermediate results:
         @image    Average image
         @saxs     Averaged saxs data
         @c2       Averaged c2 data
         @ind      Indexed data
         @n_a      Nr of averaged images
         @n_saxs   Nr of averaged saxs curves
         @n_c2     Nr of averaged c2 data
         @n_i      Nr of indexed images
         @n_q      Nr of q-rings to plot
         @n_bin    Nr of bins for size histogram


         KEYWORDS FOR PLOTS

         AVE        : Average image
         C2_IMAGE   : Heat plot of C2
         C2         : Individual C2 plots
         SAXS       : Saxs curve
         IND        : Index data

         ex: psplot -s psanaXXXX AVE C2 SAXS IND C2_IMAGE

      """

      if n_q is None :
         n_q   = min(15,len(self.q))

      if n_q > len(self.q) :            # Ascert that there is enough q's
         n_q   = len(self.q)

      if n_bin is None :
         n_bin = n_i / 10


      if image is not None :

         # Average Image
         title    = 'AVERAGE  Run ' + str(self.run_nr)
         AVEimg   = Image(n_a,title,image)
         publish.send('AVE',AVEimg)


      if saxs is not None :

         # SAXS plot
         title   = 'SAXS Run ' + str(self.run_nr)
         SAXSimg = XYPlot(n_saxs,title,self.q,saxs,xlabel='q (1/A)', formats='b')
         publish.send('SAXS',SAXSimg)

      if c2 is not None :

         # C2 plots
         title  = 'C2  Run ' + str(self.run_nr)
         # C2 heatmap plot
         C2img   = Image(n_c2,title,c2)
         publish.send('C2_IMAGE',C2img)
         # Multiplot, plot C2 for 10 q-points
         multi   = MultiPlot(n_c2,title,ncols=5)
         step    = round(len(self.q) / (n_q + 1))
         for p in range(n_q):
             R    = XYPlot(n_c2,'q = '+ str(np.around(self.q[(p+1)*step],decimals=3)),self.phi,c2[(p+1)*step],xlabel='dPhi')
             multi.add(R)
         publish.send('C2',multi)

      if ind is not None :
         if n_bin is None :
            n_bin = n_i / 10

         # Last non-zero intensity
         nz   = np.nonzero(ind[:,0])
         last = nz[0][-1]
         ind  = ind[0:last,:]

         # Check if we manged to estimate sizes
         sind    = ind[:,2] > 0.98
         if sind.any() :
            title   = 'INDEX Run ' + str(self.run_nr)
            # INDEX plot
            multi2  = MultiPlot(n_i,title,ncols=1)
            # Intensity plot
            title   = 'Intensity Run ' + str(self.run_nr)
            I       = XYPlot(n_i,title,np.arange(last),ind[:,0],xlabel='N',formats='rs')
            multi2.add(I)
            # Size plot
            title   = 'Size Run ' + str(self.run_nr)
            diam      = ind[sind,1]*(2/10) # Diameter in nm
            hist,bins = np.histogram(diam, n_bin)
            S         = Hist(n_i,title,bins,hist,xlabel='Size [nm]')
            multi2.add(S)
            publish.send('IND',multi2)
         else:
            title   = 'Intensity Run ' + str(self.run_nr)
            I       = XYPlot(n_i,title,np.arange(last),ind[:,0],xlabel='N',formats='rs')
            publish.send('IND',I)


  def store_index_h5(self, time, index,flag = 1) :
      """Store information about:
         * Time-stamp
         * Total intensity
         * Beam center
         * Estimated Particle Size
         * Estimated particle nr

      """

      self.tot_t[index]         = time.split('_')[1]

      self.tot_int[index]       = float(self.img.sum())
      self.tot_peak1_int[index] = self.peak1
      self.tot_peak2_int[index] = self.peak2
      self.tot_streak_m[index]  = self.streak_m
      self.tot_streak_s[index]  = self.streak_s

      self.tot_cx[index]        = self.cent[0]
      self.tot_cy[index]        = self.cent[1]
      self.tot_size[index]      = self.radius
      self.tot_score[index]     = self.score

      if flag :
         self.ave                  += self.img


  def store_index(self, time, index, flag = 1) :
      """Store information about:
         * Time-stamp
         * Total intensity
         * Beam center
         * Estimated Particle Size
         * Estimated particle nr

      """

      if not hasattr(self, 'streak_m'):
         self.streak_m    = 0
      if not hasattr(self, 'streak_s'):
         self.streak_s    = 0


      self.tot_t[index]         = time.time()
      self.tot_s[index]         = time.seconds()
      self.tot_ns[index]        = time.nanoseconds()
      self.tot_fd[index]        = time.fiducial()

      self.tot_int[index]       = float(self.img.sum())

      if (self.peak is not None):
          self.tot_peak1_int[index] = self.peak1
          self.tot_peak2_int[index] = self.peak2

      self.tot_streak_m[index]  = self.streak_m
      self.tot_streak_s[index]  = self.streak_s

      self.tot_cx[index]        = self.cent[0]
      self.tot_cy[index]        = self.cent[1]
      self.tot_size[index]      = self.radius
      self.tot_score[index]     = self.score

      if flag :
         self.ave                  += self.img


  def store_index2(self, time, index, flag = 1) :
      """Store information about:
         * Time-stamp
         * Total intensity
         * Beam center
         * Estimated Particle Size
         * Estimated particle nr

      """

      self.tot_t[index]         = time.time()
      self.tot_s[index]         = time.seconds()
      self.tot_ns[index]        = time.nanoseconds()
      self.tot_fd[index]        = time.fiducial()

      self.tot_int[index]       = float(self.image.sum())

      self.tot_cx[index]        = self.cent[0]
      self.tot_cy[index]        = self.cent[1]
      self.tot_size[index]      = self.radius
      self.tot_score[index]     = self.score

      if flag :
         self.ave                  += self.image



  def store_image(self,index) :
      """Store raw images

      """

      self.images[:,:,index]      = self.image
      self.cnt                   += 1


  def sum_c2(self,  flag = 0) :
      """Sum up SAXS, C2 and other quantaties continusly

         @flag    Flag as solvent [0] or signal [1]               [0/1]

      """
      # Initialize
      if (self.cnt_0 == 0.0)  and (self.cnt_1 == 0.0) :


         self.Isaxs_0     = np.zeros(self.saxs_m.shape)
         self.Vsaxs_0     = np.zeros(self.saxs_m.shape)

         self.C2_0        = np.zeros(self.c2.shape)
         self.C2sqr_0     = np.zeros(self.c2.shape)

         self.C2m_0       = np.zeros(self.c2.shape)
         self.C2msqr_0    = np.zeros(self.c2.shape)


         self.Isaxs_1     = np.zeros(self.saxs_m.shape)
         self.Vsaxs_1     = np.zeros(self.saxs_m.shape)

         self.C2_1        = np.zeros(self.c2.shape)
         self.C2sqr_1     = np.zeros(self.c2.shape)

         self.C2m_1       = np.zeros(self.c2.shape)
         self.C2msqr_1    = np.zeros(self.c2.shape)


      if flag == 0:

         self.Isaxs_0    += self.saxs_m
         self.Vsaxs_0    += (self.saxs_s)**2

         self.C2_0       += self.c2
         self.C2sqr_0    += (self.c2)**2

         self.C2m_0      += self.c2msk
         self.C2msqr_0   += (self.c2msk)**2

         self.cnt_0 += 1

      else :

         self.Isaxs_1    += self.saxs_m
         self.Vsaxs_1    += (self.saxs_s)**2

         self.C2_1       += self.c2
         self.C2sqr_1    += (self.c2)**2

         self.C2m_1      += self.c2msk
         self.C2msqr_1   += (self.c2msk)**2

         self.cnt_1 += 1



  def sum_bg2(self) :
      """Sum up BG

      """
      # Initialize
      if (self.cnt_0 == 0.0)  or (self.cnt_1 == 0.0) :

         self.Back_img    = np.zeros(self.pcimg.shape)
         self.Back_norm   = np.zeros(self.pcnorm.shape)
         self.Back_msk    = np.zeros(self.pcmsk.shape)


      self.Back_norm     +=  self.pcnorm
      self.Back_img      +=  self.pcimg
      self.Back_msk      +=  self.pcmsk



  def sum_bg(self) :
      """Sum up BG

      """
      # Initialize
      if self.cnt == 0.0  :

         self.Isaxs       = np.zeros(self.saxs_m.shape)
         self.Vsaxs       = np.zeros(self.saxs_m.shape)


         self.Back_img    = np.zeros(self.pcimg.shape)
         self.Back_norm   = np.zeros(self.pcnorm.shape)
         self.Back_msk    = np.zeros(self.pcmsk.shape)

      else:

         self.Isaxs       += self.saxs_m
         self.Vsaxs       += (self.saxs_s)**2


         self.Back_img    +=  self.pcimg
         self.Back_norm   +=  self.pcnorm
         self.Back_msk    +=  self.pcmsk



  def get_index(self, nevents, flag = 0) :
      """Generate array for storing indeces

         @nevents    Total nr of events

      """

      # Time stamp identifiers
      self.tot_t          = np.zeros(nevents)   # Time
      self.tot_s          = np.zeros(nevents)   # Second
      self.tot_ns         = np.zeros(nevents)   # Nanosecond
      self.tot_fd         = np.zeros(nevents)   # Fiducial

      # Image specific identifiers
      self.tot_int        = np.zeros(nevents)   # Total image intensity
      self.tot_peak1_int  = np.zeros(nevents)   # Radial peak 1 intensity
      self.tot_peak2_int  = np.zeros(nevents)   # Radial peak 2 intensity
      self.tot_streak_m   = np.zeros(nevents)   # Mean streak angle
      self.tot_streak_s   = np.zeros(nevents)   # Std streak angle
      self.tot_cx         = np.zeros(nevents)   # Beam x
      self.tot_cy         = np.zeros(nevents)   # Beam y
      self.tot_size       = np.zeros(nevents)   # Estimated size
      self.tot_score      = np.zeros(nevents)   # Size score (best is 1.0)

      if flag :
         # Initialize image array
         self.images      = np.zeros((self.mask.shape[0],self.mask.shape[1],nevents))
      else:
         # Initialize average image
         self.ave         = np.zeros(self.msk.shape) # Zero Image


  def get_index2(self, nevents) :
      """Generate array for storing indeces

         @nevents    Total nr of events

      """

      # Time stamp identifiers
      self.tot_t      = np.zeros(nevents)       # Time
      self.tot_s      = np.zeros(nevents)       # Second
      self.tot_ns     = np.zeros(nevents)       # Nanosecond
      self.tot_fd     = np.zeros(nevents)       # Fiducial

      # Image specific identifiers
      self.tot_int    = np.zeros(nevents)       # Total image intensity
      self.tot_cx     = np.zeros(nevents)       # Beam x
      self.tot_cy     = np.zeros(nevents)       # Beam y
      self.tot_size   = np.zeros(nevents)       # Estimated size
      self.tot_score  = np.zeros(nevents)       # Size score (best is 1.0)

      # Initialize average image

      self.ave        = np.zeros((4,512,512))   # Zero Image


  def get_size(self, plot = 0) :
      """Spherical Besselfit of SAXS data starting from radius, r_0
         Returns optimized radius in Angstrom

         @plot    Display fit between data & theory               [0/1]

      """

      self.radius, self.score    = pnccd_tbx.get_size(saxs              = self.saxs_m,
                                                      q                 = self.q,
                                                      r_i               = self.r_0,
                                                      q_i               = self.q_bound[0],
                                                      q_f               = self.q_bound[1],
                                                      plot              = plot)



  def get_c2(self, plot = 0) :
      """Computes 2-point auto-correlation using FFTs
      """

      if (self.backimg is not None) and (self.backmsk is not None) :    # Ascerts that only pixles values are considered that are  defined in both image and background

         self.pcmsk  = self.pcmsk*self.backmsk
         self.pcnorm = (self.pcnorm - self.backimg)*self.pcmsk


      self.F      = np.fft.rfft(self.pcnorm,axis=1)
      self.c2     = np.fft.irfft(self.F * self.F.conjugate()) / self.nPhi

      self.Fm     = np.fft.rfft(self.pcmsk,axis=1)
      self.c2msk  = np.fft.irfft(self.Fm * self.Fm.conjugate()) / self.nPhi


      if plot :

         m1 = np.nanmean(self.pcnorm.reshape((-1)))
         s1 = np.std(self.pcnorm.reshape((-1)))

         c  = self.c2/self.c2msk
         m2 = np.nanmean(c.reshape((-1)))
         s2 = np.std(c.reshape((-1)))

         plt.figure(5000)
         plt.clf()
         ax = plt.subplot(211)
         plt.imshow(self.pcnorm)
         plt.axis('tight')
         plt.clim(m1-3*s1,m1+3*s1)
         plt.colorbar()
         ax.set_title("Polar norm")
         ax = plt.subplot(212)
         plt.imshow(self.c2/self.c2msk)
         plt.axis('tight')
         plt.clim(m2-3*s2,m2+3*s2)
         plt.colorbar()
         ax.set_title("Angular correlation (C2)")
         plt.draw()

         if (self.cnt_0 + self.cnt_1)> 0 :

            k = np.arange(0, self.nQ, 50)

            plt.figure(6000,figsize=(20,20))
            plt.clf()
            sfig = 331
            for i in range(len(k)) :
                ax = plt.subplot(int(sfig + i))
                plt.plot(self.phi,self.c2[k[i],:]/self.c2msk[k[i],:],'b',self.phi,self.C2_0[k[i],:]/self.C2m_0[k[i],:],'r',self.phi,self.C2_1[k[i],:]/self.C2m_1[k[i],:],'g')
                plt.xlim(self.phi[0],self.phi[-1])
                m = np.median(self.C2_0[k[i],:]/self.C2m_0[k[i],:])
                s = np.std(self.C2_0[k[i],:]/self.C2m_0[k[i],:])
                plt.ylim(m-3*s,m+3*s)
                ax.set_title("q = " + str(self.q[k[i]]))

            plt.draw()

            plt.figure(7000,figsize=(10,10))
            plt.clf()
            plt.semilogy(self.q,self.saxs_m,'b',self.q,self.Isaxs_0/self.cnt_0,'r',self.q,self.Isaxs_1/self.cnt_1,'g')



  def get_norm(self,flag = 1, plot = 0) :
      """Normalizes image and mask in polar coordinates by subtracting the
         the azimuthal mean intensity. Stores the mean and std saxs intensity.

         @flag     return image norm          [0/1]

      """


      self.saxs_m           = np.ndarray(self.pcimg.shape[0])
      self.saxs_s           = np.ndarray(self.pcimg.shape[0])

      for q in range(self.pcimg.shape[0]) :

          ind               = np.nonzero(self.pcmsk[q,:])

          if len(ind[0]) == 0 :
             self.saxs_m[q]    = 0
             self.saxs_s[q]    = 0
          else:
             self.saxs_m[q]    = np.nanmean(self.pcimg[q,ind],axis=1)
             self.saxs_s[q]    = np.nanstd(self.pcimg[q,ind],axis=1)


      if flag :
          self.pcnorm         = self.pcimg - self.saxs_m[:,None]
          self.pcnorm         = self.pcnorm*self.pcmsk

      if (self.peak is not None):
          # Get Peak maxima
          self.peak1 = np.mean(self.saxs_m[self.ind1])
          self.peak2 = np.mean(self.saxs_m[self.ind2])

      if plot :

         m1 = np.nanmean(self.pcimg.reshape((-1)))
         s1 = np.std(self.pcimg.reshape((-1)))

         m2 = np.nanmean(self.pcnorm.reshape((-1)))
         s2 = np.std(self.pcnorm.reshape((-1)))

         plt.figure(4000)
         plt.clf()
         ax = plt.subplot(211)
         plt.imshow(self.pcimg)
         plt.axis('tight')
         plt.clim(0,m1+3*s1)
         plt.colorbar()
         ax.set_title("Polar image")
         ax = plt.subplot(212)
         plt.imshow(self.pcnorm)
         plt.axis('tight')
         plt.clim(m2-3*s2,m2+3*s2)
         plt.colorbar()
         ax.set_title("Polar norm")
         plt.draw()


  def get_streak_mask(self,thr = None, flag = 1, plot = 0):
      """Compute dynamic streak mask based on
         systematic intensity deviations along Phi

         @thr     Threshold (thr * std)          [pos. integer]
         @flag    Output streak angle            [0/1]
         @plot    Plot result                    [0/1]

      """

      self.pcimg        =  self.pcimg * self.pcmsk


      if thr is None :
         thr = self.thr

      if thr is not None :

         s_p                   = np.ndarray(self.pcimg.shape[1])

         self.get_norm()


         for p in range(self.pcnorm.shape[1]):

             ind               = np.nonzero(self.pcmsk[:,p])

             if len(ind[0]) == 0 :
                s_p[p]         = 0
             else:
                s_p[p]         = np.nanmean(self.pcnorm[ind,p],axis=1)

         u                     = np.nanmean(s_p)
         s                     = np.nanstd(s_p)

         A                     = ( s_p - u ) <= thr*s
         T                     = A.astype(int)

         self.pcmsk            = self.pcmsk*T
         self.pcimg            = self.pcimg * self.pcmsk


         if flag :
            # Get masked angles
            ind      = np.where(T == 0)[0]*(360/self.nPhi)
            # Select the angles for the first streak (Friedel symmetry at +180)
            if ind.size == 0 :
               self.streak_m = 0
               self.streak_s = 0
            else:
               streaks  = ind[ind < (ind[0]+90)]
               # Compute mean and std angle of the streak
               self.streak_m = np.mean(streaks)
               self.streak_s = np.std(streaks)


         if plot :

            m1 = np.nanmean(self.pcimg.reshape((-1)))
            s1 = np.std(self.pcimg.reshape((-1)))

            plt.figure(3000)
            plt.clf()
            ax = plt.subplot(211)
            plt.imshow(self.pcimg)
            plt.axis('tight')
            plt.clim(0,m1+3*s1)
            plt.colorbar()
            ax.set_title("Polar image")
            ax = plt.subplot(212)
            plt.imshow(self.pcmsk)
            plt.axis('tight')
            plt.clim(0,1)
            plt.colorbar()
            ax.set_title("Polar mask")
            plt.draw()


  def get_pixel_mask(self,thr = None, split = None, plot = 0):
      """Compute dynamic pixel mask based on
         systematic intensity deviations along Q

         @thr     Threshold (thr * std)                         [pos. integer]
         @split   Angle for splitting detector in 2 halves      [pos. integer]
         @plot    Plot output                                   [0/1]

      """

      self.pcimg        =  self.pcimg * self.pcmsk

      if thr is None :
         thr = self.thr

      if thr is not None :
         if split is not None:

             # Get index for each half panel

             i0     = int((split/360)*self.nPhi)
             i1     = i0 + int((180/360)*self.nPhi)

             pp     = np.arange(self.nPhi)
             a      = pp[i0:i1]
             b      = np.setdiff1d(pp,a)

             msk_a  = self.pcmsk[:,a]
             img_a  = self.pcimg[:,a]
             msk_b  = self.pcmsk[:,b]
             img_b  = self.pcimg[:,b]


             for q in range(self.pcimg.shape[0]) :

                 ind_a             = np.nonzero(msk_a[q,:])
                 if len(ind_a[0]) == 0 :
                    u_a            = 0
                    s_a            = 0
                 else:
                    u_a            = np.nanmean(img_a[q,ind_a],axis=1)
                    s_a            = np.nanstd(img_a[q,ind_a],axis=1)

                 ind_b             = np.nonzero(msk_b[q,:])
                 if len(ind_b[0]) == 0 :
                     u_b           = 0
                     s_b           = 0
                 else:
                     u_b           = np.nanmean(img_b[q,ind_b],axis=1)
                     s_b           = np.nanstd(img_b[q,ind_b],axis=1)

                 # Find out which panel has the smallest std, i.e which is the reference panel.

                 if s_a > s_b :
                    u   = u_b   # Use mean from panel B
                    s   = s_b   # Use std from panel B
                 else:
                    u   = u_a   # Use mean from panel A
                    s   = s_a   # Use std from panel A


                 A                 = ( img_a[q,:] - u ) <= thr*s
                 T                 = A.astype(int)
                 self.pcmsk[q,a]   = self.pcmsk[q,a]*T

                 A                 = ( img_b[q,:] - u ) <= thr*s
                 T                 = A.astype(int)
                 self.pcmsk[q,b]   = self.pcmsk[q,b]*T

             self.pcimg        = self.pcimg * self.pcmsk

         else:

             for q in range(self.pcimg.shape[0]) :

                 ind               = np.nonzero(self.pcmsk[q,:])

                 if len(ind[0])   == 0 :
                    u              = 0
                    s              = 0
                 else:
                    u              = np.nanmean(self.pcimg[q,ind],axis=1)
                    s              = np.nanstd(self.pcimg[q,ind],axis=1)

                 A                 = ( self.pcimg[q,:] - u ) <= thr*s
                 T                 = A.astype(int)
                 self.pcmsk[q,:]   = self.pcmsk[q,:]*T

             self.pcimg        = self.pcimg * self.pcmsk



         if plot :

            m1 = np.nanmean(self.pcimg.reshape((-1)))
            s1 = np.std(self.pcimg.reshape((-1)))

            plt.figure(3000)
            plt.clf()
            ax = plt.subplot(211)
            plt.imshow(self.pcimg)
            plt.axis('tight')
            plt.clim(0,m1+3*s1)
            plt.colorbar()
            ax.set_title("Polar image")
            ax = plt.subplot(212)
            plt.imshow(self.pcmsk)
            plt.axis('tight')
            plt.clim(0,1)
            plt.colorbar()
            ax.set_title("Polar mask")
            plt.draw()



  def get_polar(self, plot = 0) :
      """Returns cartesian images img & msk in polar coordinates pcimg[r,Phi] & pcmsk[r,Phi]

         @plot    Display result of grid search                   [0/1]

      """
      self.pcimg,self.pcmsk = pnccd_tbx.get_polar(img                = self.img,
                                                  msk                = self.msk,
                                                  cent               = self.cent,
                                                  r_max              = self.nQ,
                                                  r_min              = 0,
                                                  dr                 = self.dQ,
                                                  nPhi               = self.nPhi,
                                                  dPhi               = self.dPhi,
                                                  msk_a              = self.mask_angles,
                                                  msk_w              = self.mask_widths,
                                                  plot               = plot)
      # Get norm
#      scale      = sum(sum(self.pcimg))/sum(sum(self.pcmsk))
#      self.pcimg = self.pcimg/scale

      # Apply Flat-field if available
      if self.flat :
         self.pcimg = self.pcimg * self.pcflat


  def get_beam(self, angle = 45, dangle = 10, plot = 0 ) :
      """Returns estimated beam center coordintes (cent) refined using
         a grid search assuming Friedel symmetry in the image

         @ang     Center of angular slice in degrees              [pos. integer]
         @dang    Size of angular slice, ang+/-dang               [pos. integer]
         @plot    Display result of grid search                   [0/1]

      """

      self.cent           = pnccd_tbx.get_beam(img              = self.img,
                                               msk              = self.msk,
                                               r_max            = self.r_max,
                                               dr               = self.dr,
                                               cent0            = self.cent0,
                                               dx               = self.dx,
                                               dy               = self.dy,
                                               ang              = angle,
                                               dang             = dangle,
                                               plot             = plot)


  def get_image(self,evt) :
      """Retrives 2D image from event.

         @run     run number from Psana
         @time    time-stamp from Psana
      """
      self.run_nr         = int(evt.run())

      self.evt            = evt
      self.img            = self.src.image(self.evt)

      # Check if nQ > smallest dimension/2 then we extend the image with zero values
      if self.nQ > min(self.img.shape[0]/2,self.img.shape[1]/2) :
         self.img   = pnccd_tbx.extend_image(img = self.img)

      # Check if carttesian coord background image exists
      if self.cart :
         self.img = self.img - self.backimg

      # CM correction, Make sure that no other CM is already implemented

      self.img            = pnccd_tbx.common_mode(self.img, 150 , 100 , plot = 0)
      self.img            = self.img*self.msk

#      if self.detector_address == 'pnccdFront' :
#           self.img            = self.img * self.gain

  def get_calib(self,evt) :
      """Retrives 2D image from event.

         @run     run number from Psana
         @time    time-stamp from Psana
      """
      self.run_nr         = int(evt.run())

      self.evt            = evt

      self.image          = self.src.calib(self.evt)


  def get_h5(self,dfile,time) :
      """Retrives 2D image from event.

         @run     run number
         @time    time-stamp
      """

      if dfile == 0 :
         filename = self.dataset_name + '.h5'
      else :
         filename = self.dataset_name + '__' + str(dfile).zfill(4) + '.h5'

      f         = h5py.File(filename,'r')

      # Ascert that time-stamp exist in file
      if time in f:
         self.image  = f[time][self.detector_address]['HistData'].value
      else:
        self.image   = None
        self.img     = None
        return


      # Apply common mode correction

      if self.comm == 1:

         self.image  = pnccd_tbx.common_mode(img     = self.image,
                                             edge    = self.param1,
                                             side    = self.param2,
                                             plot    = self.param4)

      elif self.comm == 2:

         self.image = pnccd_tbx.common_mode_hart(img     = self.image,
                                                 msk     = self.mask,
                                                 max_int = self.param1,
                                                 max_com = self.param2,
                                                 length  = self.param3,
                                                 orient  = self.orient,
                                                 plot    = self.param4)

      # Apply geometry
      self.img  = pnccd_tbx.get_geometry(img    = self.image,
                                         gap    = self.gap,
                                         shift  = self.shift,
                                         orient = self.orient)


      # Check if nQ > smallest dimension/2 then we extend the image with zero values
      if self.nQ > min(self.img.shape[0]/2,self.img.shape[1]/2) :
         self.img   = pnccd_tbx.extend_image(img = self.img)

      # Check if carttesian coord background image exists
      if self.cart :
         self.img = self.img - self.backimg

    ###################################################################################
