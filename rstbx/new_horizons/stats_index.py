import os
from labelit.command_line.imagefiles import ImageFiles

class spotfinder_proxy:
  def __init__(self,old_spotfinder,phil,frames):
    self.frames = frames
    self.phil = phil
    self.old_pd = old_spotfinder.pd
    self.old_S = old_spotfinder
  def get_aitbx_inputs(self):
    pd = dict(xbeam = self.old_pd["xbeam"],
              ybeam = self.old_pd["ybeam"],
              osc_start = self.old_pd["osc_start"],
              binning = "1",
              size1 = self.old_pd["size1"],
              size2 = self.old_pd["size2"],
              pixel_size = self.old_pd["pixel_size"],
              distance = self.old_pd["distance"],
              wavelength = self.old_pd["wavelength"],
              deltaphi = self.old_pd["deltaphi"],
              indexing = self.old_S.get_aitbx_inputs()["indexing"],
              endstation = self.old_pd["endstation"],
              recommended_grid_sampling = self.old_S.get_aitbx_inputs()["recommended_grid_sampling"],
              twotheta = self.old_pd["twotheta"],
              resolution_inspection = self.old_pd["resolution_inspection"],
              smallest_spot_sep = self.old_S.get_aitbx_inputs()["smallest_spot_sep"],
              masks = self.old_pd["masks"], #see practical heuristics
              spot_convention = self.old_pd["spot_convention"],
              vendortype = self.old_pd["vendortype"],
              #characteristic_grid_sampling =  0.01845574110881109,
              #characteristic_resolution_mm = 43.390947190873604
             )
    pd["ref_maxcel"] = self.old_pd["ref_maxcel"] #post-get_aitbx_inputs
    self.images = self.old_S.images # sublattice average profile
    self.pd = pd

    from rstbx.new_horizons.speckfinder import speckfinder
    self.pd["indexing"]=[] # zero out the old spotfinder spots; use speckfinder spots instead
    for key in self.images.keys():
      self.specks = speckfinder(imgobj = self.frames.imageindex(key),
                       phil = self.phil,
                       inputpd = self.pd)
      self.pd["indexing"] += self.specks.get_active_data()
    return self.pd

class AutoIndexOrganizer:

  def __init__(self,verbose = 0,**kwargs):
    self.rundir = os.getcwd()
    self.verbose = verbose
    self.horizons_phil = kwargs["horizons_phil"]
    #self.horizons_phil.persist.show()
    assert kwargs.has_key('argument_module')
    self.setCommandInput(kwargs['argument_module'])
    if self.verbose: print "Process frames in directory:",self.Files.filenames.FN[0].cwd

    if kwargs.has_key('delegate'):
      self.setIndexingDelegate(kwargs['delegate'])
    self.exception_passthru = 0
    if kwargs.has_key('exception_passthru'):
      self.exception_passthru = kwargs['exception_passthru']
    print '\n'.join(self.Files.filenames())

  def setCommandInput(self,argument_module):
    self.Files = ImageFiles(argument_module,self.horizons_phil)
    self.frames = self.Files.frames()

  def printSpots(self):
    from labelit.procedure import spotfinder_and_pickle
    S = spotfinder_and_pickle(self.rundir,self.Files,
        spots_pickle = self.horizons_phil.spots_pickle,
        horizons_phil = self.horizons_phil)

    #print S.images

    NEW = spotfinder_proxy(S,self.horizons_phil,self.Files)
    NEW.images = {}
    NEW.overlapping = False
    NEW.phil_params = S.phil_params
    for frame in self.frames:
      NEW.images[frame]=dict(area=[1,]  # not actually used for new horizons
                            )

    self.S = NEW
    for frame in self.frames:
     if self.verbose:
      from labelit.command_line.stats_distl import pretty_image_stats,notes
      pretty_image_stats(S,frame)
      notes(S,self.frames[0])
    print

  def setIndexingDelegate(self,function):
    self.indexing_delegate = function

  def executeDelegate(self):
      self.info = self.indexing_delegate(self.frames,self.Files,self.S)

  def pickle_the_results(self):
      for key in ['best_integration','triclinic']:
        if self.info.has_key(key):
         if self.info[key].has_key('minimizer'): #not attained when best==tri
          del self.info[key]['minimizer'] # Must remove

         # temporary section pending an analysis of which data need to be persistent
         if self.info[key]["integration"].has_key('results'):
          #future options 1) make the whole object picklable--write test script
          #2) just pickle the data needed for the GUI
          del self.info[key]["integration"]['results']
      from labelit.dptbx.pickle_support import pickle_refinements
      pickle_refinements(self.info,self.horizons_phil.refinements_pickle)

  def process(self):
    self.printSpots()
    self.executeDelegate()
    if self.__dict__.has_key('info'): #if indexing worked
      self.pickle_the_results()
      return self.info
