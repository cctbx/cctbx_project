from scitbx.array_family import flex

import boost.python
ext = boost.python.import_ext("iotbx_detectors_ext")
from iotbx_detectors_ext import *

import exceptions
from iotbx.detectors.adsc import ADSCImage
from iotbx.detectors.mar import MARImage
from iotbx.detectors.marIP import MARIPImage

class ImageException(exceptions.Exception):
  def __init__(self,string):
    self.message = string
  def __str__(self): return self.message

all_image_types = [ADSCImage,MARImage,MARIPImage]

def ImageFactory(filename):
  for itype in all_image_types:
    try:
      I = itype(filename)
      I.readHeader()
      return I
    except:
      pass
  raise ImageException(filename+" not recognized as any known detector image type")
