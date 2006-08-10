from scitbx.array_family import flex

import boost.python
ext = boost.python.import_ext("iotbx_detectors_ext")
from iotbx_detectors_ext import *
from iotbx_detectors_bruker_ext import Bruker_base

import exceptions
from iotbx.detectors.adsc import ADSCImage
from iotbx.detectors.mar import MARImage
from iotbx.detectors.marIP import MARIPImage
from iotbx.detectors.raxis import RAXISImage
from iotbx.detectors.raxis_nonsquare import NonSquareRAXISImage
from iotbx.detectors.macscience import DIPImage
from iotbx.detectors.saturn import SaturnImage
from iotbx.detectors.bruker import BrukerImage

class ImageException(exceptions.Exception):
  def __init__(self,string):
    self.message = string
  def __str__(self): return self.message

all_image_types = [SaturnImage,DIPImage,ADSCImage,
                  MARImage,MARIPImage,RAXISImage,
                  NonSquareRAXISImage,BrukerImage]

def ImageFactory(filename):
  for itype in all_image_types:
    try:
      I = itype(filename)
      I.readHeader()

      if itype==RAXISImage:
        assert I.head['sizeFast']==I.head['sizeSlow']
        assert 0.4 < I.head['wavelength'] < 10.0 #needed to disambiguate from Bruker
      if itype==ADSCImage:
        assert I.parameters.has_key('DETECTOR_SN')
      if itype==BrukerImage:
        assert I.distance > 0.0 #needed to disambiguate from RAXIS
      return I
    except:
      pass
  raise ImageException(filename+" not recognized as any known detector image type")
