from __future__ import absolute_import, division, print_function
from six.moves import range
import math,os
from six.moves import cStringIO as StringIO
from six.moves import cPickle as pickle
from labelit.dptbx.status import cellstr
from libtbx.utils import Sorry
from libtbx.test_utils import approx_equal
from rstbx.dials_core.integration_core import show_observations

class IntegrateCharacters:
  def __init__(self,Characters,process_dictionary,horizons_phil,files,spotfinder_results):
    self.M = Characters
    self.process_dictionary = process_dictionary
    self.horizons_phil = horizons_phil
    self.files = files
    self.spotfinder_results = spotfinder_results
    self.triclinic = self.M.best()[-1]
    fres = ResLimitControl(self.process_dictionary,self.horizons_phil)

    self.triclinic['integration'] = self.integrate_one_character(
      setting=self.triclinic,
      integration_limit=fres.current_limit)
    #return # Enforces legacy behavior--no recycling to expand the integration limit
            # Comment this "return" in for testing without the macrocycle
    #With appropriate safeguards, macrocycle gives better resolution estimate:
    A = ResolutionAnalysisMetaClass(self.triclinic['integration'], self.horizons_phil)
    print(A)
    safety_counter = 2
    while A.retest_required() and safety_counter > 0 and \
          A.target_resol() < fres.current_limit:
      safety_counter -= 1
      tolerance = 2.E-2
      if (not approx_equal(A.target_resol(),fres.outer_limit,
             eps=tolerance,out=None)) and A.target_resol() < fres.outer_limit:
        break
      fres.current_limit = A.target_resol()

      trial = self.integrate_one_character(
        setting = self.triclinic,
        integration_limit = fres.current_limit)

      results = trial["results"]
      obs = [item.get_obs(trial["spacegroup"]) for item in results]
      A.stats_mtz(trial,obs)
      #print A
    if 'trial' in vars().keys(): self.triclinic['integration'] = trial

  def integrate_one_character(self,setting,integration_limit):
    #from libtbx.development.timers import Profiler
    #P = Profiler("Preliminary")
    import copy
    local = copy.deepcopy(self.process_dictionary)
    local['cell']=cellstr(setting)

    print("Cell in setting",setting["counter"],local["cell"])

    frames = list(sorted(self.spotfinder_results.pd['osc_start'].keys()))

    local['maxcel']='0'
    local['xbeam']="%f"%setting['minimizer'].new['xbeam']
    local['ybeam']="%f"%setting['minimizer'].new['ybeam']
    local['distance']="%f"%setting['minimizer'].new['distance']
    local["resolution"]= "%f"%integration_limit

    from labelit.steps import primaries
    local['spacegroup'] = primaries[setting['bravais']]

    local['procstart'] = local['procend'] = "%d"%frames[0]

    self.pixel_size = float(local['pixel_size'])

    from labelit.dptbx import AutoIndexEngine, Parameters
    ai = AutoIndexEngine(local['endstation'])

    P = Parameters(xbeam=setting["refined x beam"],ybeam=setting["refined y beam"],
             distance=setting["refined distance"],twotheta=float(local["twotheta"]))
    ai.setBase(P)
    ai.setWavelength(float(local['wavelength']))
    ai.setMaxcell(float(local['ref_maxcel']))
    print("Deltaphi is",float(local['deltaphi']))
    ai.setDeltaphi(float(local['deltaphi'])*math.pi/180.)
    ai.setMosaicity(setting["mosaicity"])
    ai.setOrientation(setting["orient"])
    refimage = self.files.images[0]
    ai.set_active_areas(self.horizons_phil,
                        beam=(int(refimage.beamx/refimage.pixel_size),
                              int(refimage.beamy/refimage.pixel_size)))

    image_centers = [(math.pi/180.)*float(x) for x in local["osc_start"].values()]

    print("Limiting resolution",integration_limit)
    local["results"] = []
    for i in range(len(frames)):
      print("---------BEGIN Integrate one frame %d %s" % \
          (frames[i], os.path.split(self.files.filenames()[i])[-1]))
      #P = Profiler("worker")
      if self.horizons_phil.integration.combine_sym_constraints_and_3D_target and setting["counter"]>1:
        from rstbx.apps.stills.dials_refinement_preceding_integration import integrate_one_frame
        integrate_worker = integrate_one_frame(self.triclinic["integration"]["results"][0])
      else:
        from rstbx.apps.stills.deltapsi_refinement_preceding_integration import integrate_one_frame
        integrate_worker = integrate_one_frame()
      integrate_worker.inputai = ai

      integrate_worker.inputpd = dict(masks=local["masks"],
                                      size1=local["size1"],
                                      size2=local["size2"],
                                      symmetry=setting["best_subsym"])
        # carefully select only the data items needed for integrate_worker
        # avoid giving the whole process dictionary; reference to "local"
        # is a circular reference creating memory leak, while copying the
        # whole thing is a big performance hit.
      integrate_worker.frame_numbers = frames
      integrate_worker.imagefiles = self.files
      integrate_worker.spotfinder = self.spotfinder_results
      integrate_worker.image_centers = image_centers

      integrate_worker.limiting_resolution = integration_limit
      integrate_worker.setting_id = setting["counter"]
      integrate_worker.pixel_size = self.pixel_size
      integrate_worker.set_pixel_size(self.pixel_size)
      integrate_worker.set_detector_size(int(local["size1"]),int(local["size2"]))

      integrate_worker.set_detector_saturation(refimage.saturation)
      integrate_worker.set_up_mask_focus()
      integrate_worker.initialize_increments(i)
      integrate_worker.horizons_phil = self.horizons_phil
      if self.horizons_phil.indexing.verbose_cv:
        print("EFFECTIVE TILING"," ".join(
          ["%d"%z for z in refimage.get_tile_manager(self.horizons_phil).effective_tiling_as_flex_int()]))
      integrate_worker.integration_concept(image_number = i,
        cb_op_to_primitive = setting["cb_op_inp_best"].inverse(),
        verbose_cv = self.horizons_phil.indexing.verbose_cv,
        background_factor = self.horizons_phil.integration.background_factor,
        )
      #P = Profiler("proper")
      integrate_worker.integration_proper()
      local["results"].append(integrate_worker)
      local["r_xbeam"]=ai.xbeam()
      local["r_ybeam"]=ai.ybeam()
      local["r_distance"]=ai.distance()
      local["r_wavelength"]=ai.wavelength
      local["r_residual"]=integrate_worker.r_residual
      local["r_mosaicity"]=setting["mosaicity"]
      try:
        local["ewald_proximal_volume"]=integrate_worker.ewald_proximal_volume
      except Exception as e:
        local["ewald_proximal_volume"]=None

      if (self.horizons_phil.indexing.open_wx_viewer):
       if True: #use updated slip viewer
        try:
          import wx
          from rstbx.slip_viewer.frame import XrayFrame as SlipXrayFrame
          from rstbx.command_line.slip_viewer import master_str as slip_params
          from iotbx import phil
          from spotfinder import phil_str
          from spotfinder.command_line.signal_strength import additional_spotfinder_phil_defs

          work_phil = phil.process_command_line("",master_string=slip_params + phil_str + additional_spotfinder_phil_defs)
          work_params = work_phil.work.extract()

          app = wx.App(0)
          wx.SystemOptions.SetOption("osx.openfiledialog.always-show-types", "1")
          frame = SlipXrayFrame(None, -1, "X-ray image display", size=(800,720))
          frame.Show()

          # Update initial settings with values from the command line.  Needs
          # to be done before image is loaded (but after the frame is
          # instantiated).
          frame.inherited_params = integrate_worker.horizons_phil
          frame.params = work_params

          if (frame.pyslip is None):
            frame.init_pyslip()
          if (frame.settings_frame is None):
            frame.OnShowSettings(None)
          frame.Layout()

          frame.pyslip.tiles.user_requests_antialiasing = work_params.anti_aliasing
          frame.settings_frame.panel.center_ctrl.SetValue(True)
          frame.settings_frame.panel.integ_ctrl.SetValue(True)
          frame.settings_frame.panel.spots_ctrl.SetValue(False)
          frame.settings.show_effective_tiling = work_params.show_effective_tiling
          frame.settings_frame.panel.collect_values()
          paths = work_phil.remaining_args

          frame.user_callback = integrate_worker.slip_callback
          frame.load_image(self.files.filenames()[i])

          app.MainLoop()
          del app
        except Exception:
          pass # must use phenix.wxpython for wx display

       elif False : #original wx viewer
        try:
          from rstbx.viewer.frame import XrayFrame
          import wx
          from rstbx.viewer import display
          display.user_callback = integrate_worker.user_callback

          app = wx.App(0)
          frame = XrayFrame(None, -1, "X-ray image display", size=(1200,1080))
          frame.settings.show_spotfinder_spots = False
          frame.settings.show_integration = False
          #frame.settings.enable_collect_values = False
          frame.SetSize((1024,780))
          frame.load_image(self.files.filenames()[i])
          frame.Show()
          app.MainLoop()
          del app
        except Exception:
          pass # must use phenix.wxpython for wx display

      # for the wx image viewer
      filename = self.horizons_phil.indexing.indexing_pickle
      if filename != None:
        filename = "%s_%d_%d.pkl"%(filename,setting["counter"],keys[i])

        SIO = StringIO()
        table_raw = show_observations(integrate_worker.get_obs(
          local["spacegroup"]),out=SIO)
        limitobject = ResolutionAnalysisMetaClass(local, self.horizons_phil)
        info = dict(table = SIO.getvalue(),
          table_raw = table_raw,
          xbeam = setting["refined x beam"],
          ybeam = setting["refined y beam"],
          distance = setting["refined distance"],
          residual = integrate_worker.r_residual,
          resolution = limitobject.value, # FIXME not reliable?
          mosaicity = setting["mosaicity"],
          pointgroup = local["spacegroup"],
          hkllist = integrate_worker.hkllist,
          predictions = (1./integrate_worker.pixel_size)*integrate_worker.predicted,
          mapped_predictions = integrate_worker.detector_xy,
          integration_masks_xy = integrate_worker.integration_masks_as_xy_tuples(),
          background_masks_xy = integrate_worker.background_masks_as_xy_tuples()
        )
        assert info["predictions"].size() >= info["mapped_predictions"].size()
        assert info["predictions"].size() == info["hkllist"].size()
        G = open(filename,"wb")
        pickle.dump(info,G,pickle.HIGHEST_PROTOCOL)
      print("---------END Integrate one frame",frames[i])

    return local

  def find_best(self):
    self.best_counter = 1

    for index in self.M.best()[0:len(self.M.best())-1]:
      if 'status' in index and index['status'] in [
        'unlikely','very_unlikely']:continue
      if float(self.triclinic['integration']['resolution'])==0.0:
        raise Sorry("No signal detected in triclinic integration trial")

      if self.horizons_phil.known_setting == index['counter'] and \
         self.horizons_phil.integration.montecarlo_integration_limit is not None:
        index['integration'] = self.integrate_one_character(
        setting=index,
        integration_limit=self.horizons_phil.integration.montecarlo_integration_limit)
      elif self.horizons_phil.integration.greedy_integration_limit:
        index['integration'] = self.integrate_one_character(
        setting=index,
        integration_limit=float(self.triclinic['integration']['results'][0].limiting_resolution))
      else:
        index['integration'] = self.integrate_one_character(
        setting=index,
        integration_limit=float(self.triclinic['integration']['resolution']))

      A = ResolutionAnalysisMetaClass(index['integration'],self.horizons_phil)
      print(A)
      if (self.horizons_phil.known_cell!=None or
         self.horizons_phil.known_symmetry!=None):
        self.best_counter = index['counter']
        break

      if ( float(index['integration']['r_residual']) <
           self.horizons_phil.mosflm_rmsd_tolerance *
           float(self.triclinic['integration']['r_residual']) ):
        self.best_counter = index['counter']
        break

  def save_best(self):
    file = self.horizons_phil.indexing.completeness_pickle
    for index in self.M.best():
      if 'integration' in index:
        if index['counter']==self.best_counter:
          local = index["integration"]
          info = dict(
            xbeam = local["r_xbeam"],
            ybeam = local["r_ybeam"],
            distance = local["r_distance"],
            wavelength = float(local["r_wavelength"]),
            residual = local["r_residual"],
            mosaicity = local["r_mosaicity"],
            pointgroup = local["spacegroup"],
            observations = [a.get_obs(local["spacegroup"]) for a in local["results"]],
            mapped_predictions = [a.detector_xy for a in local["results"]],
            model_partialities = [getattr(a,"partialities",None) for a in local["results"]],
            sa_parameters = [getattr(a,"best_params","None") for a in local["results"]],
            max_signal = [getattr(a,"max_signal",None) for a in local["results"]],
            current_orientation = [getattr(a,"current_orientation",None) for a in local["results"]],
            current_cb_op_to_primitive = [getattr(a,"current_cb_op_to_primitive",None) for a in local["results"]],
            correction_vectors = [getattr(a, 'correction_vectors', None)
                                  for a in local['results']],
            effective_tiling = self.files.images[0].get_tile_manager(
              self.horizons_phil).effective_tiling_as_flex_int()
          )
          for correction_type in self.horizons_phil.integration.absorption_correction:
            if correction_type.apply and correction_type.algorithm=="fuller_kapton":
              info['fuller_kapton_absorption_correction'] = [a.fuller_kapton_absorption_correction for a in local["results"]]
          if self.horizons_phil.integration.model=="user_supplied":
            info['ML_half_mosaicity_deg'] = [getattr(a,"ML_half_mosaicity_deg",0) for a in local["results"]]
            info['ML_domain_size_ang'] = [getattr(a,"ML_domain_size_ang",0) for a in local["results"]]
            info['ewald_proximal_volume'] = [getattr(a,"ewald_proximal_volume",0) for a in local["results"]]
          info["identified_isoform"] = local["results"][0].__dict__.get("identified_isoform",None)
          if file is not None:
            G = open(file,"wb")
            pickle.dump(info,G,pickle.HIGHEST_PROTOCOL)
          return info

  def show(self):
    print()
    print("New Horizons Integration results:")
    print("Solution  SpaceGroup Beam x   y  distance  Resolution Mosaicity RMS")
    for index in self.M.best():
      if 'integration' in index:
        limitobject = ResolutionAnalysisMetaClass( index['integration'], self.horizons_phil )
        if index['counter']==self.best_counter:
          print(":)", end=' ')
          self.process_dictionary['best_integration']=index
        else: print("  ", end=' ')
        # only write out the triclinic integration results if there is
        #  an application for the data--future expansion
        if index['counter']==1 and len(self.M.best())>1:
          self.process_dictionary['triclinic']=index
        print("%3d"%index['counter'], end=' ')
        print("%12s"%index['integration']["spacegroup"], end=' ')
        print("%6.2f %6.2f"%(float(index['integration']["r_xbeam"]),
                             float(index['integration']["r_ybeam"])), end=' ')
        print("%7.2f   "%(float(index['integration']["r_distance"])), end=' ')
        if limitobject.value == 0.00:
          #Analysis of mtz file gives no resolution estimate; revert to limit detected by spotpicking
          index['integration']['r_resolution'] = float(self.process_dictionary["resolution_inspection"])
        else: index['integration']['r_resolution'] = limitobject.value
        print("%7.2f   "%index['integration']['r_resolution'], end=' ')
        print(index['integration']["r_mosaicity"], end=' ')
        print("  ", end=' ')
        print("%5.3f"%index['integration']["r_residual"])

        continue # the following code merely demonstrates the unpacking of partiality info
        if hasattr(index["integration"]["results"][0],"partialities"):
          hackobs = index["integration"]["results"][0].get_obs(index["integration"]["spacegroup"])
          hackpart = index["integration"]["results"][0].partialities["data"]
          hackhkl = list(index["integration"]["results"][0].partialities["indices"])
          from scitbx.array_family import flex
          xx = flex.double()
          yy = flex.double()
          for idx in range(hackobs.indices().size()):
            hkl = hackobs.indices()[idx]
            thisobs = hackobs.data()[idx]
            lookupidx = hackhkl.index(hkl)
            print(hkl,thisobs,hackhkl[lookupidx],hackpart[lookupidx])
            if thisobs>0.:
              resolution = hackobs.unit_cell().d(hkl)
              if resolution > 2.5 and resolution < 4.0:
                # correlation between partiality & Iobs only within resolution shells
                xx.append(math.log(thisobs))
                yy.append(hackpart[lookupidx])
          from matplotlib import pyplot as plt
          plt.plot(xx,yy,"r.")
          plt.show()


    #Sublattice analysis
    for index in self.M.best():
      if index['counter']==1 and 'integration' in index:
        results = index['integration']["results"]
        obs = [item.get_obs(index['integration']["spacegroup"]) for item in results]
        try:
          get_limits(params = index['integration'],
                   file = obs,
                   verbose = False,
                   sublattice_flag = True,
                   override_maximum_bins = 12,
                   horizons_phil = self.horizons_phil)
        except: # intentional
          #Numpy multiarray.error raises an object not derived from Exception
          print("Catch any problem with sublattice analysis & numpy masked arrays")
          return

class limits_fix_engine:
  def __init__(self):
    pass
  def rawprocess(self,first_pass):
    self.xbeam = float(first_pass.get('labelit_x',first_pass["xbeam"]))
    self.ybeam = float(first_pass.get('labelit_y',first_pass["ybeam"]))
    self.px = float(first_pass["pixel_size"])
    self.sz1 = float(first_pass.get("size1"))
    self.sz2 = float(first_pass.get("size2"))

  def corners(self):
    #coordinates of the four detector corners relative to the beam (mm)
    return ((-self.xbeam,                  -self.ybeam),
            (-self.xbeam,                  self.px*self.sz2-self.ybeam),
            (self.px*self.sz1-self.xbeam,-self.ybeam),
            (self.px*self.sz1-self.xbeam,self.px*self.sz2-self.ybeam))
  def edges(self):
    #coordinates of the four detector edges relative to the beam (mm)
    return ((-self.xbeam+(self.px*self.sz1/2.)  ,-self.ybeam),
            (-self.xbeam,                         (self.px*self.sz2/2.) - self.ybeam),
            (self.px*self.sz1-self.xbeam,       (self.px*self.sz2/2.)-self.ybeam),
            ((self.px*self.sz1/2.) - self.xbeam,self.px*self.sz2-self.ybeam))

class SafetyLimits:
  def __init__(self,INFO,AI):
    self.INFO = INFO
    self.AI = AI

  def determined_limit_from_screen(self):
    try:
      screen = float(self.INFO['best_integration']['integration']['r_resolution'])
    except Exception:
      # no previous integration success.  Rely on DISTL results.
      screen = float(self.INFO['resolution'])
    return screen

  def safety_limit(self,algorithm="corner"):
    epsilon = 0.0001
    # Originally, the limit at edge of detector on a MarCCD, otherwise lambda/2.
    # But due to comments from Ana Gonzalez on problems with MOSFLM integration
    # out to the corner (Gordon Conference, 2006), a) a limit will be placed
    # on all detectors based on detector geometry, and b) it will be a
    # configurable choice whether the limit is based on the farthest corner
    # or the farthest edge.
    XE = limits_fix_engine()
    XE.rawprocess(self.INFO)
    from labelit.mathsupport import length
    if algorithm=="edge":
      mpoint = max([length(point) for point in XE.edges()])
    elif algorithm=="corner":
      mpoint = max([length(point) for point in XE.corners()])
    theta = math.atan2(mpoint,self.AI.distance())/2.0
    return epsilon + self.AI.wavelength/(2.0*math.sin(theta))

  def determined_limit_with_safety(self):
    return max(self.determined_limit_from_screen(), self.safety_limit() )

class ResLimitControl:
  def __init__(self,pd,horizons_phil):
    #initialVolumeFactor = 1.8 # Legacy value=1.5 prior to 7/17/07
    initialVolumeFactor = horizons_phil.integration.initial_volume_factor
    self.subsequentVolumeFactor = 1.5 # this value not used by stats_mtz; see stats_mtz code
    if horizons_phil.mosflm_integration_reslimit_override!=None:
      self.initial_limit = horizons_phil.mosflm_integration_reslimit_override
      self.outer_limit = self.initial_limit
    else:
      #resolution_inspection from DISTL is a conservative estimate
      #increase the reciprocal space volume by factor of 1.5 for trial integration
      from math import pow
      trial_reslimit = float(pd["resolution_inspection"])/pow(initialVolumeFactor,1.0/3.0)

      from labelit.dptbx import AutoIndexEngine,Parameters
      ai = AutoIndexEngine(pd['endstation'],
           horizons_phil.model_refinement_minimum_N) #Just a container for a few parameters
      base = Parameters(xbeam = float(pd['xbeam']), ybeam = float(pd['ybeam']),
                        distance = float(pd['distance']), twotheta = 0.0 )
      ai.setBase(base)
      ai.setWavelength(float(pd['wavelength']))

      safety = SafetyLimits(pd,ai).safety_limit(
        algorithm = horizons_phil.mosflm_safety_algorithm)

      self.outer_limit = safety
      self.initial_limit = max((trial_reslimit,safety))
      """synopsis:
      LABELIT integrates to determine the best-estimate of
      dataset resolution.  After integration out to the image corner, an
      intensity log-plot is used to linearly extrapolate out to an I/sigma
      of 0.75, based on all partials and fulls.  Ana Gonzalez requested that
      this be changed to the image edge, because of MOSFLM problems in cases
      where the diffraction is of low quality or low resolution.  Unfortunately,
      using a cutoff at the edge is counter-productive in cases where diffraction
      extends past the corners.  In those cases, the log-linear extrapolation
      invariably leads to a resolution estimate that is more optimistic (lower
      Angstrom cutoff) than when all of the data is used out to the image corner.
      Edge cutoff is therefore not recommended.  However, it is now available as
      a configurable choice using the command line argument:
      mosflm_safety_algorithm=[corner|edge]  (default = corner).
      #print pd['file'][pd['file'].keys()[0]]
      #print "corner: %.2f"%(max((trial_reslimit,interface.SafetyLimits(pd,ai).safety_limit(algorithm="corner"))))
      #print "  edge: %.2f"%(max((trial_reslimit,interface.SafetyLimits(pd,ai).safety_limit(algorithm="edge"))))
      """
    self.current_limit = self.initial_limit
  def show(self):
    print("initial resolution limit",self.initial_limit)
    print("outer resolution limit",self.outer_limit)
  def limits(self):
    yield self.current_limit
    while True:
      self.current_limit /=  pow(self.subsequentVolumeFactor,1.0/3.0)
      if self.current_limit <= self.outer_limit: break
      yield self.current_limit

from labelit.diffraction.stats_mtz import get_limits
class ResolutionAnalysisMetaClass(get_limits):
  def __init__(self,integration_dict,horizons_phil,verbose=False):
    self.integration_dict = integration_dict
    self.horizons_phil = horizons_phil
    results = self.integration_dict["results"]
    obs = [item.get_obs(self.integration_dict["spacegroup"]) for item in results]
    if verbose:
      for item in results:
        show_observations(item.get_obs(self.integration_dict["spacegroup"]))
    self.stats_mtz(self.integration_dict,obs)
    self.integration_dict["resolution"] = self.value

  def stats_mtz(self,integration_dict,file):
    try:
      get_limits.__init__(self,params=integration_dict,file=file,horizons_phil=self.horizons_phil)
    except: # intentional
      #Numpy multiarray.error raises an object not derived from Exception
      print("Catch any problem with stats mtz & numpy masked arrays")

  def retest_required(self):
    return self.status.require_expanded_limit != None

  def target_resol(self):
    tmp = self.status.require_expanded_limit.split(" ")
    return float(tmp[3])

  def target_Isig(self):
    tmp = self.status.require_expanded_limit.split(" ")
    return float(tmp[9])

  def __getattr__(self,key):
    if key=='value':
      return self.status.value
    else:
      return self.__dict__[key]
