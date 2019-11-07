"""Main idea:  having already done
      1) indexing & integration --> allresults
      2) metrology assessment --> mysql store tag_spotfinder
      3) cxi.merge --> mysql store tag_frame
   ...now do a detailed metrology refinement to simultaneously optimize
   metrology and crystal orientation.
"""
from __future__ import absolute_import, division, print_function
from six.moves import range
from cctbx.array_family import flex
import iotbx.phil
import math
from scitbx import matrix
from xfel import get_radial_tangential_vectors
from xfel.merging.database.merging_database import manager
from xfel import correction_vector_store
from libtbx.development.timers import Timer

from xfel.merging.database.merging_database import mysql_master_phil
from six.moves import zip
master_phil="""
bravais_setting_id = None
  .type = int
  .help = ID number for the Bravais setting of interest (Labelit format).  eg, 1=triclinic, 12=hexagonal
show_plots = False
  .type = bool
  .help = Show graphical plots using matplotlib
show_consistency = False
  .type = bool
  .help = Run the consistency controls
effective_tile_boundaries = None
  .type = ints
  .help = effective integer tile boundaries applied to convert xtc stream to pickled image files. Must be 64 * 4 integers.
max_frames = 0
  .type = int
  .help = for the SQL query, maximum number of frames to return (generally for development testing only). 0->return all frames
""" + mysql_master_phil

def consistency_controls(DATA,params,annotate=False):#DATA is an instance of correction_vectors()
  PIXEL_SZ = 0.11 # mm/pixel
  CART = manager(params)
  db = CART.connection()
  cursor = db.cursor()

  for iframe in range(len(DATA.FRAMES["frame_id"])):
    frame = DATA.FRAMES["frame_id"][iframe]
    selection = (DATA.frame_id == frame)
    match_count = selection.count(True)
    if match_count>0:
      print(frame, DATA.frame_id.select(selection)[0], end=' ') # frame number
      frame_beam_x = DATA.FRAMES["beam_x"][iframe]
      obs_beam_x = DATA.refined_cntr_x.select(selection)[0] * PIXEL_SZ
      print("%7.3f"%(frame_beam_x - obs_beam_x), end=' ') # agreement of beam_x in mm
      frame_beam_y = DATA.FRAMES["beam_y"][iframe]
      obs_beam_y = DATA.refined_cntr_y.select(selection)[0] * PIXEL_SZ
      print("%7.3f"%(frame_beam_y - obs_beam_y), end=' ') # agreement of beam_y in mm
      #...The labelit-refined direct beam position agrees with CV_listing logfile output

      file_name = DATA.FRAMES["unique_file_name"][iframe]

      cursor.execute("SELECT COUNT(*) FROM %s_observation WHERE frame_id_0_base=%d-1;"%(params.mysql.runtag,frame))
      integrated_observations = cursor.fetchall()[0][0]
      print("%4d <? %4d"%(match_count,integrated_observations),file_name, end=' ')

      cursor.execute(
        """SELECT t1.detector_x, t1.detector_y, t1.original_h, t1.original_k, t1.original_l
           FROM %s_observation AS t1
           WHERE t1.frame_id_0_base=%d-1
        """%(
           params.mysql.runtag,frame))
      fetched = cursor.fetchall()
      detector_x = [a[0] for a in fetched]
      detector_y = [a[1] for a in fetched]
      spotfx = DATA.spotfx.select(selection)
      spotfy = DATA.spotfy.select(selection)
      spotcx = DATA.spotcx.select(selection)
      spotcy = DATA.spotcy.select(selection)
      hkl = DATA.HKL.select(selection)
      integrated_hkl = [(int(a[2]),int(a[3]),int(a[4])) for a in fetched]


      # Now compute the displacement between integrated and calculated spot position.
      # presumably tells us about the empirical nudge factor.
      sq_displace = flex.double()
      sq_cv = flex.double()
      for icalc,calc_hkl in enumerate(hkl):
        try:
          jinteg = integrated_hkl.index(calc_hkl)
          sq_displace.append(  (spotcx[icalc]-detector_x[jinteg])**2 + (spotcy[icalc]-detector_y[jinteg])**2 )
        except ValueError: pass
        sq_cv.append( (spotcx[icalc]-spotfx[icalc])**2 + (spotcy[icalc]-spotfy[icalc])**2 )
      if len(sq_displace) > 2:
        print("rmsd=%7.3f"%math.sqrt(flex.mean(sq_displace)), end=' ')
      else:
        print("rmsd    None", end=' ')
      rmsd_cv = math.sqrt(flex.mean(sq_cv))
      print("cv%7.3f"%rmsd_cv)

      if params.show_plots is True:
        import os
        os.environ["BOOST_ADAPTBX_FPE_DEFAULT"]="1"
        from matplotlib import pyplot as plt
        plt.figure(figsize=(9,9))
        plt.plot(spotcx,spotcy,
          markerfacecolor="g",marker=".",markeredgewidth=0,linestyle="None")
        plt.plot(spotfx, spotfy,
          markerfacecolor="r",marker=".",markeredgewidth=0,linestyle="None")
        plt.plot(detector_x,detector_y,
          markerfacecolor="b",marker=".",markeredgewidth=0,linestyle="None")
        if annotate:
          for idx in range(len(spotfx)):
            plt.annotate("%s"%str(hkl[idx]), xy=(spotfx[idx],spotfy[idx]),
                         xytext=None, xycoords="data", textcoords="data", arrowprops=None,
                         color="red",size=8)
            plt.annotate("%s"%str(hkl[idx]), xy=(spotcx[idx],spotcy[idx]),
                         xytext=None, xycoords="data", textcoords="data", arrowprops=None,
                         color="green",size=8)
          for idx in range(len(fetched)):
            plt.annotate("%s"%str(integrated_hkl[idx]), xy=(detector_x[idx],detector_y[idx]),
                         xytext=None, xycoords="data", textcoords="data", arrowprops=None,
                         color="blue",size=8)
        plt.axes().set_aspect("equal")
        plt.show()
        #...confirms that integrated spot position (observation table) aligns with
        #   spotfinder spot position (spotfinder table) in both approximate position
        #   (match is not exact due to empirical positional nudge factor that is applied
        #   after spotfinder table and before observation table)
        #   and Miller index tag.

class correction_vectors(correction_vector_store):

 def get_obs_from_mysql(self,params):
   T = Timer("database")
   CART = manager(params)
   db = CART.connection()
   cursor = db.cursor()
   cursor.execute("SELECT DISTINCT frame_id FROM %s_spotfinder;"%params.mysql.runtag)
   AAA = cursor.fetchall()
   print("From the CV log file text output there are %d distinct frames with spotfinder spots"%len(AAA))

   if params.max_frames==0:
     cursor.execute("SELECT * FROM %s_spotfinder;"%params.mysql.runtag)
   else:
     cursor.execute("SELECT * FROM %s_spotfinder WHERE frame_id<%d;"%(
       params.mysql.runtag, params.max_frames))
   return cursor.fetchall()

 def get_frames_from_mysql(self,params):
   T = Timer("frames")
   CART = manager(params)
   db = CART.connection()
   cursor = db.cursor()
   cursor.execute("SELECT * FROM %s_frame;"%params.mysql.runtag)
   ALL = cursor.fetchall()
   from cctbx.crystal_orientation import crystal_orientation
   orientations = [crystal_orientation(
     (a[8],a[9],a[10],a[11],a[12],a[13],a[14],a[15],a[16]),False) for a in ALL]
   return dict( frame_id = flex.int( [a[0] for a in ALL] ),
               wavelength = flex.double( [a[1] for a in ALL] ),
                   beam_x = flex.double( [a[2] for a in ALL] ),
                   beam_y = flex.double( [a[3] for a in ALL] ),
                 distance = flex.double( [a[4] for a in ALL] ),
              orientation = orientations,
          rotation100_rad = flex.double([a[17] for a in ALL] ),
          rotation010_rad = flex.double([a[18] for a in ALL] ),
          rotation001_rad = flex.double([a[19] for a in ALL] ),
       half_mosaicity_deg = flex.double([a[20] for a in ALL] ),
              wave_HE_ang = flex.double([a[21] for a in ALL] ),
              wave_LE_ang = flex.double([a[22] for a in ALL] ),
          domain_size_ang = flex.double([a[23] for a in ALL] ),
         unique_file_name = [a[24] for a in ALL],
  )

 def delrsq_functional(self,calcx,calcy):
   delx = calcx - self.spotfx
   dely = calcy - self.spotfy
   delrsq = delx*delx + dely*dely
   self.last_functional = delrsq
   return delrsq

 INCIDENT_BEAM = (0.,0.,-1.)
 DETECTOR_NORMAL = (0.,0.,-1.)

 @staticmethod
 def standalone_check(self,setting_id,entry,d,cutoff):

    wavelength = (d['wavelength'])
    beam_x = (d['xbeam'])
    beam_y = (d['ybeam'])
    distance = (d['distance'])
    orientation = (d['current_orientation'][0])

    print("testing frame....................",entry)

    for cv in d['correction_vectors'][0]:

      from rstbx.bandpass import use_case_bp3, parameters_bp3
      from scitbx.matrix import col
      from math import hypot, pi
      indices = flex.miller_index()
      indices.append(cv['hkl'])
      parameters = parameters_bp3(
        indices=indices,
        orientation=orientation,
        incident_beam=col(self.INCIDENT_BEAM),
        packed_tophat=col((1.,1.,0.)),
        detector_normal=col(self.DETECTOR_NORMAL),
        detector_fast=col((0.,1.,0.)),detector_slow=col((1.,0.,0.)),
        pixel_size=col((0.11,0.11,0)), # XXX hardcoded, twice!
        pixel_offset=col((0.,0.,0.0)),
        distance=distance,
        detector_origin=col((-beam_x,-beam_y,0))
      )
      ucbp3 = use_case_bp3(parameters=parameters)
      ucbp3.set_active_areas(self.tiles)
      integration_signal_penetration=0.5
      ucbp3.set_sensor_model(thickness_mm=0.5,
                             mu_rho=8.36644, # CS_PAD detector at 1.3 Angstrom
                             signal_penetration=integration_signal_penetration)

      ucbp3.set_mosaicity(0.)
      ucbp3.set_bandpass(wavelength,
                         wavelength)
      ucbp3.set_orientation(orientation)
      ucbp3.set_domain_size(5000.)

      ucbp3.picture_fast_slow_force()

      ucbp3_prediction = 0.5 * (ucbp3.hi_E_limit + ucbp3.lo_E_limit)
      diff = hypot(ucbp3_prediction[0][0] - cv['predspot'][1],
                   ucbp3_prediction[0][1] - cv['predspot'][0])

      if diff > cutoff:
        print("Correction vector too long: %6.2f pixels; ignore image or increase diff_cutoff (current value=%5.1f)"%(diff,cutoff))
        return False

      # For some reason, the setting_id is recorded for each
      # correction vector as well--assert that it is consistent.
      #if cv['setting_id'] != setting_id:
      #  print "HATTNE BIG SLIPUP 2"
      if not cv['setting_id'] == setting_id: return False

      # For each observed spot, figure out what tile it is on, and
      # store in itile.  XXX This is probably not necessary here, as
      # correction_vector_store::register_line() does the same thing.
      obstile = None
      for i in range(0, len(self.tiles), 4):
        if     cv['obsspot'][0] >= self.tiles[i + 0] \
           and cv['obsspot'][0] <= self.tiles[i + 2] \
           and cv['obsspot'][1] >= self.tiles[i + 1] \
           and cv['obsspot'][1] <= self.tiles[i + 3]:
          obstile = i
          break
      if obstile is None: return False

      spotfx = (cv['obsspot'][0])
      spotfy = (cv['obsspot'][1])
      spotcx = (cv['predspot'][0])
      spotcy = (cv['predspot'][1])
      correction_vector_x = spotcx - spotfx
      correction_vector_y = spotcy - spotfy
      length = hypot(correction_vector_x, correction_vector_y)
      if length > 8:
        print("LENGTH SLIPUP",length)
        return False

    return True

 def read_data(self,params):
  from os import listdir, path
  from libtbx import easy_pickle
  from cctbx.crystal_orientation import crystal_orientation # XXX Necessary later?

  #directory = "/net/viper/raid1/hattne/L220/merging/05fs"
  #directory = "/reg/d/psdm/cxi/cxib8113/scratch/sauter/metrology/008"
  #directory = "/reg/d/psdm/xpp/xpp74813/scratch/sauter/metrology/204"
  #directory = "/net/viper/raid1/hattne/L220/merging/test"
  #directory = "/reg/d/psdm/xpp/xppb4313/scratch/brewster/results/r0243/003/integration"
  #directory = "/reg/d/psdm/cxi/cxic0614/scratch/sauter/metrology/004/integration"
  #directory = "/reg/d/psdm/cxi/cxic0614/scratch/sauter/metrology/150/integration"
  #directory = "/reg/d/psdm/cxi/cxic0614/scratch/sauter/metrology/152/integration"
  directory = "/reg/d/psdm/cxi/cxib6714/scratch/sauter/metrology/009/integration"
  dir_glob = "/reg/d/psdm/CXI/cxib6714/scratch/sauter/results/r*/009/integration"
  dir_glob = "/reg/d/psdm/CXI/cxib6714/scratch/sauter/results/r*/801/integration"
  dir_glob = "/reg/d/psdm/xpp/xpp74813/scratch/sauter/r*/216/integration"
  dir_glob = "/reg/d/psdm/xpp/xpp74813/ftc/sauter/result/r*/104/integration"
  dir_glob = "/reg/d/psdm/cxi/cxid9114/scratch/sauter/metrology/001/integration"
  dir_glob = "/reg/d/psdm/CXI/cxid9114/ftc/brewster/results/r00[3-4]*/003/integration"
  dir_glob = "/reg/d/psdm/CXI/cxid9114/ftc/sauter/results/r00[3-4]*/004/integration"
  dir_glob = "/reg/d/psdm/CXI/cxid9114/ftc/sauter/results/r00[3-4]*/006/integration"
  dir_list = ["/reg/d/psdm/CXI/cxid9114/ftc/brewster/results/r%04d/006/integration"%seq for seq in range(95,115)]
  dir_list = ["/reg/d/psdm/CXI/cxid9114/ftc/sauter/results/r%04d/018/integration"%seq for seq in range(102,115)]
  dir_list = params.data

  T = Timer("populate C++ store with register line")

  itile = flex.int()
  self.spotfx = flex.double()
  self.spotfy = flex.double()
  self.spotcx = flex.double()
  self.spotcy = flex.double()
  self.observed_cntr_x = flex.double()
  self.observed_cntr_y = flex.double()
  self.refined_cntr_x = flex.double()
  self.refined_cntr_y = flex.double()
  self.HKL = flex.miller_index()
  self.radial = flex.double()
  self.azimut = flex.double()

  self.FRAMES = dict(
    frame_id=flex.int(),
    wavelength=flex.double(),
    beam_x=flex.double(),
    beam_y=flex.double(),
    distance=flex.double(),
    orientation=[],
    rotation100_rad=flex.double(),
    rotation010_rad=flex.double(),
    rotation001_rad=flex.double(),
    half_mosaicity_deg=flex.double(),
    wave_HE_ang=flex.double(),
    wave_LE_ang=flex.double(),
    domain_size_ang=flex.double(),
    unique_file_name=[]
  )

  self.frame_id = flex.int()
  import glob
  #for directory in glob.glob(dir_glob):
  for directory in dir_list:
   if self.params.max_frames is not None and len(self.FRAMES['frame_id']) >= self.params.max_frames:
      break
   for entry in listdir(directory):
    tttd = d = easy_pickle.load(path.join(directory, entry))

    # XXX Hardcoded, should honour the phil!  And should be verified
    # to be consistent for each correction vector later on!
    #import pdb; pdb.set_trace()
    setting_id = d['correction_vectors'][0][0]['setting_id']

    #if setting_id != 5:
    #if setting_id != 12:
    if setting_id != self.params.bravais_setting_id:
    #if setting_id != 22:
      #print "HATTNE BIG SLIPUP 1"
      continue

    # Assert that effective_tiling is consistent, and a non-zero
    # multiple of eight (only whole sensors considered for now--see
    # mark10.fit_translation4.print_table()).  self.tiles is
    # initialised to zero-length in the C++ code.  XXX Should now be
    # able to retire the "effective_tile_boundaries" parameter.
    #
    # XXX Other checks from correction_vector plot, such as consistent
    # setting?
    if hasattr(self, 'tiles') and len(self.tiles) > 0:
      assert (self.tiles == d['effective_tiling']).count(False) == 0
    else:
      assert len(d['effective_tiling']) > 0 \
        and  len(d['effective_tiling']) % 8 == 0
      self.tiles = d['effective_tiling']

    if not self.standalone_check(self,setting_id,entry,d,params.diff_cutoff): continue

    # Reading the frame data.  The frame ID is just the index of the
    # image.
    self.FRAMES['frame_id'].append(len(self.FRAMES['frame_id']) + 1) # XXX try zero-based here
    self.FRAMES['wavelength'].append(d['wavelength'])
    self.FRAMES['beam_x'].append(d['xbeam'])
    self.FRAMES['beam_y'].append(d['ybeam'])
    self.FRAMES['distance'].append(d['distance'])
    self.FRAMES['orientation'].append(d['current_orientation'][0])
    self.FRAMES['rotation100_rad'].append(0) # XXX FICTION
    self.FRAMES['rotation010_rad'].append(0) # XXX FICTION
    self.FRAMES['rotation001_rad'].append(0) # XXX FICTION
    self.FRAMES['half_mosaicity_deg'].append(0) # XXX FICTION
#    self.FRAMES['wave_HE_ang'].append(0.995 * d['wavelength']) # XXX FICTION -- what does Nick use?
#    self.FRAMES['wave_LE_ang'].append(1.005 * d['wavelength']) # XXX FICTION
    self.FRAMES['wave_HE_ang'].append(d['wavelength'])
    self.FRAMES['wave_LE_ang'].append(d['wavelength'])
    self.FRAMES['domain_size_ang'].append(5000) # XXX FICTION
    self.FRAMES['unique_file_name'].append(path.join(directory, entry))

    print("added frame", self.FRAMES['frame_id'][-1],entry)


    for cv in d['correction_vectors'][0]:

      # Try to reproduce every predicition using the model from the
      # frame -- skip CV if fail.  Could be because of wrong HKL:s?
      #
      # Copy these two images to test directory to reproduce:
      #  int-s01-2011-02-20T21:27Z37.392_00000.pickle
      #  int-s01-2011-02-20T21:27Z37.725_00000.pickle
      from rstbx.bandpass import use_case_bp3, parameters_bp3
      from scitbx.matrix import col
      from math import hypot, pi
      indices = flex.miller_index()
      indices.append(cv['hkl'])
      parameters = parameters_bp3(
        indices=indices,
        orientation=self.FRAMES['orientation'][-1],
        incident_beam=col(self.INCIDENT_BEAM),
        packed_tophat=col((1.,1.,0.)),
        detector_normal=col(self.DETECTOR_NORMAL),
        detector_fast=col((0.,1.,0.)),detector_slow=col((1.,0.,0.)),
        pixel_size=col((0.11,0.11,0)), # XXX hardcoded, twice!
        pixel_offset=col((0.,0.,0.0)),
        distance=self.FRAMES['distance'][-1],
        detector_origin=col((-self.FRAMES['beam_x'][-1],
                             -self.FRAMES['beam_y'][-1],
                             0))
      )
      ucbp3 = use_case_bp3(parameters=parameters)
      ucbp3.set_active_areas(self.tiles)
      integration_signal_penetration=0.5
      ucbp3.set_sensor_model(thickness_mm=0.5,
                             mu_rho=8.36644, # CS_PAD detector at 1.3 Angstrom
                             signal_penetration=integration_signal_penetration)
      half_mosaicity_rad = self.FRAMES['half_mosaicity_deg'][-1] * pi/180.
      ucbp3.set_mosaicity(half_mosaicity_rad)
      ucbp3.set_bandpass(self.FRAMES['wave_HE_ang'][-1],
                         self.FRAMES['wave_LE_ang'][-1])
      ucbp3.set_orientation(self.FRAMES['orientation'][-1])
      ucbp3.set_domain_size(self.FRAMES['domain_size_ang'][-1])

      ucbp3.picture_fast_slow_force()

      ucbp3_prediction = 0.5 * (ucbp3.hi_E_limit + ucbp3.lo_E_limit)
      diff = hypot(ucbp3_prediction[0][0] - cv['predspot'][1],
                   ucbp3_prediction[0][1] - cv['predspot'][0])

      if diff > self.params.diff_cutoff:
        print("HATTNE INDEXING SLIPUP")
        continue

      # For some reason, the setting_id is recorded for each
      # correction vector as well--assert that it is consistent.
      #if cv['setting_id'] != setting_id:
      #  print "HATTNE BIG SLIPUP 2"
      assert cv['setting_id'] == setting_id

      # For each observed spot, figure out what tile it is on, and
      # store in itile.  XXX This is probably not necessary here, as
      # correction_vector_store::register_line() does the same thing.
      obstile = None
      for i in range(0, len(self.tiles), 4):
        if     cv['obsspot'][0] >= self.tiles[i + 0] \
           and cv['obsspot'][0] <= self.tiles[i + 2] \
           and cv['obsspot'][1] >= self.tiles[i + 1] \
           and cv['obsspot'][1] <= self.tiles[i + 3]:
          obstile = i
          break
      assert obstile is not None
      itile.append(obstile) # XXX unused variable?

      # ID of current frame.
      self.frame_id.append(self.FRAMES['frame_id'][-1])

      self.spotfx.append(cv['obsspot'][0])
      self.spotfy.append(cv['obsspot'][1])
      self.spotcx.append(cv['predspot'][0])
      self.spotcy.append(cv['predspot'][1])

      self.observed_cntr_x.append(cv['obscenter'][0])
      self.observed_cntr_y.append(cv['obscenter'][1])
      self.refined_cntr_x.append(cv['refinedcenter'][0])
      self.refined_cntr_y.append(cv['refinedcenter'][1])

      self.HKL.append(cv['hkl'])

      self.azimut.append(cv['azimuthal'])
      self.radial.append(cv['radial'])
    #print self.FRAMES['frame_id'][-1]
    # Should honour the max_frames phil parameter
    #if len(self.FRAMES['frame_id']) >= 1000:
    if self.params.max_frames is not None and \
      len(self.FRAMES['frame_id']) >= self.params.max_frames:
      break

    """
For 5000 first images:
STATS FOR TILE 14
  sel_delx           -6.59755265524 -4.41676757746e-10 5.7773557278
  sel_dely           -6.30796620634 -8.3053734774e-10 6.3362200841
  symmetric_offset_x -6.5975526548 -2.73229417105e-15 5.77735572824
  symmetric_offset_y -6.30796620551 1.16406818748e-15 6.33622008493
  symmetric rsq      0.000255199593417 2.95803352999 56.1918083904
  rmsd               1.71989346472

For 10000 first images:
STATS FOR TILE 14
  sel_delx           -6.92345292727 6.9094552919e-10 611.497770006
  sel_dely           -6.39690476093 1.1869355797e-09 894.691806871
  symmetric_offset_x -6.92345292796 1.28753258216e-14 611.497770005
  symmetric_offset_y -6.39690476212 -2.10251420168e-15 894.69180687
  symmetric rsq      1.58067791823e-05 30.3331143761 1174402.952
  rmsd               5.50755066941
    """


  # This is mark3.fit_translation2.nominal_tile_centers()
  self.To_x = flex.double(len(self.tiles) // 4)
  self.To_y = flex.double(len(self.tiles) // 4)
  for x in range(len(self.tiles) // 4):
    self.To_x[x] = (self.tiles[4 * x + 0] + self.tiles[4 * x + 2]) / 2
    self.To_y[x] = (self.tiles[4 * x + 1] + self.tiles[4 * x + 3]) / 2


  delx = self.spotcx - self.spotfx
  dely = self.spotcy - self.spotfy
  self.delrsq = self.delrsq_functional(calcx = self.spotcx, calcy = self.spotcy)

  self.initialize_per_tile_sums()
  self.tile_rmsd = [0.]*(len(self.tiles) // 4)
  self.asymmetric_tile_rmsd = [0.]*(len(self.tiles) // 4)


  # XXX Is (beam1x, beam1y) really observed center and (beamrx,
  # beamry) refined center?  Nick thinks YES!
  #
  #itile2 = flex.int([self.register_line(a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9]) for a in ALL])
  itile2 = flex.int(
    [self.register_line(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7])
     for a in zip(self.observed_cntr_x, self.observed_cntr_y,
                  self.refined_cntr_x, self.refined_cntr_y,
                  self.spotfx, self.spotfy,
                  self.spotcx, self.spotcy)])
  if params.show_consistency: consistency_controls(self,params)

  T = Timer("calcs based on C++ store")
  self.selections = []
  self.selection_counts = []
  for x in range(len(self.tiles) // 4):
      if self.tilecounts[x]==0:
        self.radii[x] = 0
        self.mean_cv[x] = matrix.col((0, 0))
      else:
        self.radii[x]/=self.tilecounts[x]
        self.mean_cv[x] = matrix.col(self.mean_cv[x]) / self.tilecounts[x]

      selection = (self.master_tiles == x)
      self.selections.append(selection)
      selected_cv = self.master_cv.select(selection)
      self.selection_counts.append(selected_cv.size()) # for curvatures

      if len(selected_cv)>0:
        self.asymmetric_tile_rmsd[x] = math.sqrt(flex.mean (self.delrsq.select(selection)))
        sel_delx = delx.select(selection)
        sel_dely = dely.select(selection)
        symmetric_offset_x = sel_delx - self.mean_cv[x][0]
        symmetric_offset_y = sel_dely - self.mean_cv[x][1]
        symmetricrsq = symmetric_offset_x*symmetric_offset_x + symmetric_offset_y*symmetric_offset_y

        self.tile_rmsd[x] =math.sqrt(flex.mean(symmetricrsq))
      else:
        self.asymmetric_tile_rmsd[x]=0.
        self.tile_rmsd[x]=0.

  self.overall_N = flex.sum(flex.int( [int(t) for t in self.tilecounts] ))
  self.overall_cv = matrix.col(self.overall_cv)/self.overall_N
  self.overall_rmsd = math.sqrt( self.sum_sq_cv / self.overall_N )

  # master weights for mark3 calculation takes 0.3 seconds
  self.master_weights = flex.double(len(self.master_tiles))
  self.largest_sample = max(self.tilecounts)
  for x in range(len(self.tiles) // 4):
    self.master_weights.set_selected( self.selections[x], self.tile_weight(x))

  print("AFTER read     cx,     cy", flex.mean(self.spotcx), flex.mean(self.spotcy))
  print("AFTER read     fx,     fy", flex.mean(self.spotfx), flex.mean(self.spotfy))
  print("AFTER read rmsd_x, rmsd_y", math.sqrt(flex.mean(flex.pow(self.spotcx - self.spotfx, 2))), \
                                     math.sqrt(flex.mean(flex.pow(self.spotcy - self.spotfy, 2))))

  return


 def tile_weight(self,idx):
   # at most, a tile is allow to be upweighted 10x due to low sample count, but no more.
   # A tile with no observations carries zero weight.
   if self.tilecounts[idx] == 0:
     return 0
   return min(10.,self.largest_sample / self.tilecounts[idx])

 def print_table(self):
  from libtbx import table_utils
  from libtbx.str_utils import format_value

  table_header = ["Tile","Dist","Nobs","aRmsd","Rmsd","delx","dely","disp","rotdeg","Rsigma","Tsigma"]
  table_data = []
  table_data.append(table_header)
  sort_radii = flex.sort_permutation(flex.double(self.radii))
  tile_rmsds = flex.double()
  radial_sigmas = flex.double(len(self.tiles) // 4)
  tangen_sigmas = flex.double(len(self.tiles) // 4)
  for idx in range(len(self.tiles) // 4):
    x = sort_radii[idx]
    if self.tilecounts[x] < 3:
      wtaveg = 0.0
      radial = (0,0)
      tangential = (0,0)
      rmean,tmean,rsigma,tsigma=(0,0,1,1)
    else:
      wtaveg = self.weighted_average_angle_deg_from_tile(x)
      radial,tangential,rmean,tmean,rsigma,tsigma = get_radial_tangential_vectors(self,x)

    radial_sigmas[x]=rsigma
    tangen_sigmas[x]=tsigma
    table_data.append(  [
      format_value("%3d",   x),
      format_value("%7.2f", self.radii[x]),
      format_value("%6d",  self.tilecounts[x]),
      format_value("%5.2f", self.asymmetric_tile_rmsd[x]),
      format_value("%5.2f", self.tile_rmsd[x]),
      format_value("%5.2f", self.mean_cv[x][0]),
      format_value("%5.2f", self.mean_cv[x][1]),
      format_value("%5.2f", matrix.col(self.mean_cv[x]).length()),
      format_value("%6.2f", wtaveg),
      format_value("%6.2f", rsigma),
      format_value("%6.2f", tsigma),
    ])
  table_data.append([""]*len(table_header))
  rstats = flex.mean_and_variance(radial_sigmas,self.tilecounts.as_double())
  tstats = flex.mean_and_variance(tangen_sigmas,self.tilecounts.as_double())
  table_data.append(  [
      format_value("%3s",   "ALL"),
      format_value("%s", ""),
      format_value("%6d",  self.overall_N),
      format_value("%5.2f", math.sqrt(flex.mean(self.delrsq))),
      format_value("%5.2f", self.overall_rmsd),
      format_value("%5.2f", self.overall_cv[0]),
      format_value("%5.2f", self.overall_cv[1]),
      format_value("%5.2f", flex.mean(flex.double([matrix.col(cv).length() for cv in self.mean_cv]))),
      format_value("%s", ""),
      format_value("%6.2f", rstats.mean()),
      format_value("%6.2f", tstats.mean()),
    ])

  print()
  print(table_utils.format(table_data,has_header=1,justify='center',delim=" "))

#-----------------------------------------------------------------------
def get_phil(args):
  phil = iotbx.phil.process_command_line(args=args, master_string=master_phil).show()
  work_params = phil.work.extract()
  if ("--help" in args) :
    libtbx.phil.parse(master_phil.show())
    return

  if work_params.show_plots is True:
    from matplotlib import pyplot as plt # special import
  return work_params

def run(args):

  work_params = get_phil(args)
  C = correction_vectors()
  C.read_data(work_params)
  C.print_table()

  return None

if (__name__ == "__main__"):

  result = run(args=["mysql.runtag=for_may060corner","mysql.passwd=sql789",
                     "mysql.user=nick","mysql.database=xfelnks",
#                     "effective_tile_boundaries=700, 447, 894, 632, 505, 447, 699, 632, 700, 658, 894, 843, 505, 658, 699, 843, 496, 26, 681, 220, 495, 224, 680, 418, 706, 27, 891, 221, 706, 224, 891, 418, 75, 241, 269, 426, 272, 241, 466, 426, 74, 27, 268, 212, 270, 28, 464, 213, 95,  454, 280, 648, 96, 650, 281, 844, 309, 455, 494, 649, 308, 650, 493, 844, 433, 860, 618, 1054, 433, 1056, 618, 1250, 643, 860, 828, 1054, 643, 1055, 828, 1249, 15, 1077, 209, 1262, 212, 1076, 406, 1261, 14, 865, 208, 1050, 212, 865, 406, 1050, 228, 1485, 413, 1679, 228, 1288, 413, 1482, 15, 1494, 200, 1688, 15, 1290, 200, 1484, 439, 1473, 633, 1658, 635, 1473, 829, 1658, 440, 1261, 634, 1446, 636, 1261, 830, 1446, 842, 1133, 1036, 1318, 1037, 1133, 1231, 1318, 841, 922, 1035, 1107, 1036, 922, 1230, 1107, 1055, 1541, 1240, 1735, 1055, 1344, 1240, 1538, 843, 1542, 1028, 1736, 843, 1345, 1028, 1539, 1464, 1338, 1658, 1523, 1267, 1337, 1461, 1522, 1467, 1550, 1661, 1735, 1269, 1550, 1463, 1735, 1454, 1117, 1639, 1311, 1453, 921, 1638, 1115, 1241, 1117, 1426, 1311, 1241, 922, 1426, 1116, 1120, 716, 1305, 910, 1121, 521, 1306, 715, 910, 717, 1095, 911, 910, 522, 1095, 716, 1533, 515, 1727, 700, 1336, 514, 1530, 699, 1532, 726, 1726, 911, 1335, 725, 1529, 910, 1329, 93, 1514, 287, 1327, 291, 1512, 485, 1543, 91, 1728, 285, 1540, 290, 1725, 484, 1106, 113, 1300, 298, 911, 113, 1105, 298, 1105, 325, 1299, 510, 910, 325, 1104, 510".replace(" ",""),
#                                               713, 437, 907, 622, 516, 438, 710, 623, 713, 650, 907, 835, 516, 650, 710, 835, 509, 18, 694, 212, 507, 215, 692, 409, 721, 19, 906, 213, 720, 215, 905, 409, 86, 231, 280, 416, 283, 231, 477, 416, 85, 19, 279, 204, 283, 19, 477, 204, 106, 444, 291, 638, 106, 640, 291, 834, 318, 443, 503, 637, 318, 640, 503, 834, 434, 849, 619, 1043, 436, 1046, 621, 1240, 647, 848, 832, 1042, 648, 1045, 833, 1239, 18, 1066, 212, 1251, 214, 1065, 408, 1250, 17, 853, 211, 1038, 213, 853, 407, 1038, 229, 1474, 414, 1668, 229, 1277, 414, 1471, 15, 1474, 200, 1668, 16, 1278, 201, 1472, 442, 1464, 636, 1649, 638, 1464, 832, 1649, 441, 1252, 635, 1437, 638, 1252, 832, 1437, 846, 1134, 1040, 1319, 1042, 1133, 1236, 1318, 845, 922, 1039, 1107, 1042, 922, 1236, 1107, 1060, 1542, 1245, 1736, 1060, 1346, 1245, 1540, 848, 1543, 1033, 1737, 847, 1348, 1032, 1542, 1469, 1336, 1663, 1521, 1272, 1337, 1466, 1522, 1472, 1550, 1666, 1735, 1274, 1549, 1468, 1734, 1460, 1117, 1645, 1311, 1460, 921, 1645, 1115, 1247, 1117, 1432, 1311, 1248, 921, 1433, 1115, 1130, 718, 1315, 912, 1130, 522, 1315, 716, 918, 719, 1103, 913, 917, 523, 1102, 717, 1543, 514, 1737, 699, 1346, 513, 1540, 698, 1543, 725, 1737, 910, 1346, 725, 1540, 910, 1338, 94, 1523, 288, 1339, 290, 1524, 484, 1552, 93, 1737, 287, 1551, 289, 1736, 483, 1115, 114, 1309, 299, 918, 113, 1112, 298, 1115, 326, 1309, 511, 918, 326, 1112, 511".replace(" ",""),
  ])
