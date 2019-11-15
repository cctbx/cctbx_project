from __future__ import absolute_import, division, print_function
import os
from spotfinder.applications.stats_distl import pretty_image_stats,notes

def spotfinder_factory(absrundir,frames,phil_params):

  local_frames=frames.frames()

  A = frames.images[0]
  #A.readHeader()--deprecate this because it squashes any overrides
  #                from dataset_preferences processed in imagefiles.py
  pd = {'directory':frames.filenames.FN[0].cwd,
        'template': frames.filenames.FN[0].template,
        'identifier':frames.filenames.FN[0].fileroot,
        'vendortype':A.vendortype,
        'binning':'%d'%A.bin,
        'distance':'%f'%A.distance,
        'wavelength':'%f'%A.wavelength,
        'deltaphi':'%f'%A.deltaphi,
        }

  #temp values for getting coordinate convention
  pd['pixel_size']='%f'%A.pixel_size
  pd['size1']='%f'%A.size1
  pd['size2']='%f'%A.size2
  pd['ybeam'] = '%f'%A.beamy
  pd['xbeam'] = '%f'%A.beamx
  try:
    pd['twotheta'] = '%f'%A.twotheta
  except Exception:
    pd['twotheta'] = '0.0'

  from spotfinder.applications.practical_heuristics import heuristics_base
  Spotfinder = heuristics_base(pd,phil_params)

  for framenumber in local_frames:
    try:
      assert framenumber in Spotfinder.images
    except Exception:
      Spotfinder.register_frames(framenumber,frames)
      if phil_params.spotfinder_verbose: Spotfinder.show()

  return Spotfinder

def dxtbx_spotfinder_factory(phil_params):

  import dxtbx.format.Registry
  reader = dxtbx.format.Registry.get_format_class_for_file(phil_params.distl.image[0])
  from spotfinder.dxtbx_toolbox.practical_heuristics import heuristics_base
  Spotfinder = heuristics_base(phil_params)
  return Spotfinder

class DistlOrganizer(object):

  def __init__(self,verbose = 0,**kwargs):
    self.rundir = os.getcwd()
    self.verbose = verbose
    self.phil_params = kwargs["phil_params"]
    if self.phil_params.distl.dxtbx:
      self.set_dxtbx_input()
      return
    if 'argument_module' in kwargs:
      # new interface
      self.setCommandInput(kwargs['argument_module'])

  def set_dxtbx_input(self):
    pass

  def setCommandInput(self,argument_module):
    from spotfinder.diffraction.imagefiles import spotfinder_image_files as ImageFiles
    self.Files = ImageFiles(argument_module,self.phil_params)
    self.frames = self.Files.frames()

  def update_spotfinder(self):
    # used by distl.image_viewer
    S = spotfinder_factory(self.rundir,self.Files,self.phil_params)
    self.S = S
    for frame in self.frames:
      if self.verbose:
        pretty_image_stats(S,frame)
        notes(S,self.frames[0])

  def printSpots(self):
    '''spotfinder and pickle implicitly assumes ADSC format'''
    if self.phil_params.distl.dxtbx:
      self.S = S = dxtbx_spotfinder_factory(self.phil_params)
    else:
      self.S = S = spotfinder_factory(self.rundir,self.Files,self.phil_params)
      for frame in self.frames:
        if self.verbose:
          pretty_image_stats(S,frame)
          notes(S,self.frames[0])
