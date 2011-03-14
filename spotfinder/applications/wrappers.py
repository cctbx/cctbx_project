import os
from labelit.command_line.imagefiles import ImageFiles
from labelit.preferences import procedure_preferences
from labelit import tnear2
from spotfinder.applications.stats_distl import pretty_image_stats,notes

def spotfinder_factory(absrundir,frames):

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
  except:
    pd['twotheta'] = '0.0'

  Spotfinder = tnear2.sf2(pd)

  for framenumber in local_frames:
    try:
      assert Spotfinder.images.has_key(framenumber)
    except:
      Spotfinder.register_frames(framenumber,frames)
      if procedure_preferences.spotfinder_verbose: Spotfinder.show()

  return Spotfinder

class DistlOrganizer:

  def __init__(self,verbose = 0,**kwargs):
    self.rundir = os.getcwd()
    self.verbose = verbose
    if kwargs.has_key('argument_module'):
      # new interface
      self.setCommandInput(kwargs['argument_module'])
    #print '\n'.join(self.Files.filenames())

  def setCommandInput(self,argument_module):
    # new way of doing things
    self.Files = ImageFiles(argument_module)
    self.frames = self.Files.frames()

  def printSpots(self):
    '''spotfinder and pickle implicitly assumes ADSC format'''
    S = spotfinder_factory(self.rundir,self.Files)
    self.S = S
    for frame in self.frames:
     if self.verbose:
      pretty_image_stats(S,frame)
      notes(S,self.frames[0])
