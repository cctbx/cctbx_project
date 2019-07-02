"""Specialization code for 2011 JCSG pilot experiment; not general use"""
from __future__ import absolute_import, division, print_function
from six.moves import range
from six.moves import cPickle as pickle

import math
from labelit.preferences import labelit_commands,labelit_phil
from labelit.command_line.default_param import establish_dict_for_refinement
from labelit.dptbx.autoindex import index_and_refine
from labelit.command_line.stats_index import best_character_to_IndexPrinter
from labelit.command_line.stats_index import AutoIndexOrganizer

class Empty:pass

class api:
  def __init__(self,*args):
    #application programmer interface for screen, using identical inputs
    # to the command line interface
    args = list(args)
    labelit_phil.merge_command_line(args)
    E = Empty()
    E.argv=['Empty']
    for x in range(len(args)):
      E.argv.append(args[x])
    E.horizons_phil=labelit_commands

    self.establish_dict_for_refinement = establish_dict_for_refinement
    self.index_and_refine = index_and_refine
    self.best_character_to_IndexPrinter = best_character_to_IndexPrinter

    self.Org = AutoIndexOrganizer(
      verbose=labelit_commands.distl.bins.verbose,argument_module=E)
    self.Org.setIndexingDelegate(self.index_and_integrate)
      # legacy:algorithm control could be excercised here
    self.horizons_phil = labelit_commands

  def __call__(self):
    return self.Org.process()

  def index_and_integrate(self,frames,files,spotfinder_results):
    self.frames = frames
    self.files = files
    self.spotfinder_results = spotfinder_results
    pd = self.establish_dict_for_refinement(frames,spotfinder_results)
    #------------------------------------------------------------
    ai,P = self.index_and_refine(pd,files,spotfinder_results,0)
    self.indexing_ai = ai
    #------------------------------------------------------------
    if labelit_commands.compatibility_allow==False:
      M = self.best_character_to_IndexPrinter(ai,P,pd,horizons_phil=self.horizons_phil)
    else:
      from labelit.diffraction.compatibility import \
        best_compatibility_to_IndexPrinter
      M = best_compatibility_to_IndexPrinter(ai,P,pd,files,spotfinder_results,horizons_phil=self.horizons_phil)
    #------------------------------------------------------------
    if "writer" in labelit_commands.__dict__:
      labelit_commands.writer.make_image_plots_detail(
        ai=ai,pd=pd,inframes=files,spotfinder_results=spotfinder_results)
    if not labelit_commands.index_only:
      IC = IntegrateCharacters(M,pd)
      IC.write_mosflm_matrices()
      IC.find_best()
      IC.show()
    return pd

  def create_case_only(self,frames,file,spotfinder_results):
    self.pd = self.establish_dict_for_refinement(frames)

  def analyze_one(self,solution):
    inputpd = self.Org.process()

    settings = pickle.load(open("LABELIT_possible","rb"))
    setting = [setting for setting in settings if setting["counter"]==solution][0]

    from labelit.preferences import labelit_commands as param

    pixel_size = float(inputpd['pixel_size'])
    self.pixel_size = pixel_size

    from labelit.dptbx import AutoIndexEngine, Parameters
    ai = AutoIndexEngine(inputpd['endstation'])

    P = Parameters(xbeam=setting["refined x beam"],ybeam=setting["refined y beam"],
             distance=setting["refined distance"],twotheta=float(inputpd["twotheta"]))
    ai.setBase(P)
    ai.setWavelength(float(inputpd['wavelength']))
    ai.setMaxcell(float(inputpd['ref_maxcel']))
    print("Deltaphi is",float(inputpd['deltaphi']))
    ai.setDeltaphi(float(inputpd['deltaphi'])*math.pi/180.)
    ai.setMosaicity(setting["mosaicity"])
    ai.setOrientation(setting["orient"])
    #why aren't hexagonal constraints applied here???
    print(inputpd["osc_start"])

    image_centers = [(math.pi/180.)*float(x) for x in inputpd["osc_start"].values()]

    limiting_resolution = param.distl_highres_limit
    print("Limiting resolution",limiting_resolution)

    #predict the spots
    spots = ai.predict_all(image_centers[0],limiting_resolution)
    pre2m = spots.vec3()
    self.pre2m = pre2m

    hkllist = spots.hkl()
    cell = ai.getOrientation().unit_cell()
    print(cell)
    for hkl in hkllist:
      #print "%25s %5.2f"%(str(hkl),cell.d(hkl))
      assert cell.d(hkl)>=limiting_resolution
    print("Number of hkls:",(hkllist).size(), end=' ')
    print("all inside the %4.2f Angstrom limiting sphere."%limiting_resolution)
    print("The unit cell is",cell)
    self.solution_setting_ai = ai
    self.solution_pd = inputpd
    self.image_centers = image_centers
    self.one_setting = setting
    return [ai.getOrientation().unit_cell(),hkllist]

  def parameters(self):
    image_centers = [(math.pi/180.)*float(x) for x in self.solution_pd["osc_start"].values()]
    P = dict( image_center_radians = image_centers,
              wavelength = float(self.solution_pd['wavelength']),
              deltaphi = float(self.solution_pd['deltaphi']),
              mosaicity_degrees = self.solution_setting_ai.getMosaicity(),
              orientation = self.solution_setting_ai.getOrientation(),
              )
    return P

from libtbx.development.timers import Timer
from rstbx.dials_core.integration_core import integration_core
from rstbx.apps.stills.simple_integration import IntegrationMetaProcedure as base_class
class IntegrationMetaProcedure(base_class):
  def __init__(self,inputs,backcompat_horizons_phil): # inputs is an instance of class api
    self.horizons_phil = backcompat_horizons_phil
    integration_core.__init__(self)
    self.inputai = inputs.solution_setting_ai #C++ autoindex engine
    self.indexing_ai = inputs.indexing_ai

    #Note...here we may be working with a non-primitive cell, so
    # must work with systematic absences...not implemented yet.
    if self.indexing_ai.getData().size() < 40: return # initial protection

    self.inputpd = inputs.solution_pd #parameter dictionary
    self.inputpd["symmetry"] = inputs.one_setting["best_subsym"]
    self.inputframes = inputs.frames
    self.imagefiles = inputs.files
    self.spotfinder = inputs.spotfinder_results
    self.frame_numbers = self.spotfinder.pd['osc_start'].keys()
    self.frame_numbers.sort()

    self.image_centers = inputs.image_centers
    print("IMAGE CENTERS",self.image_centers)

    # initial resolution from DISTL
    resolution_est = float(self.inputpd['resolution_inspection'])
    print("initial resolution from DISTL",resolution_est)

    # resolution limit of the strong spots used for indexing
    resolution_str = self.indexing_ai.high().reciprocal_spacing
    resolution = max(resolution_est,resolution_str)
    print("resolution limit of the strong spots used for indexing",resolution)
    self.limiting_resolution = resolution

    #print "resolution: %.2f %.2f"%(resolution_est,resolution_str)
    self.pixel_size = inputs.pixel_size
    self.set_pixel_size(inputs.pixel_size)
    self.set_detector_size(int(self.inputpd["size1"]),int(self.inputpd["size2"]))

    self.pre2m = inputs.pre2m
    self.set_up_mask_focus()
    self.initialize_increments()
    T = Timer("concept")
    from cctbx import sgtbx
    self.integration_concept(image_number = 0,
        cb_op_to_primitive = sgtbx.change_of_basis_op(), #identity; supports only primitive lattices
        verbose_cv = self.horizons_phil.indexing.verbose_cv,
        background_factor = self.horizons_phil.integration.background_factor)
    T = Timer("proper")
    self.integration_proper()
