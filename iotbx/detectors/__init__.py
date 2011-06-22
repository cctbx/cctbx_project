import os

import boost.python
ext = boost.python.import_ext("iotbx_detectors_ext")
from iotbx_detectors_ext import *
from iotbx_detectors_bruker_ext import Bruker_base # import dependency

import exceptions
class ImageException(exceptions.Exception):
  pass

from iotbx.detectors.adsc import ADSCImage
from iotbx.detectors.mar import MARImage
from iotbx.detectors.marIP import MARIPImage
from iotbx.detectors.cbf import CBFImage
from iotbx.detectors.dtrek import DTREKImage
from iotbx.detectors.raxis import RAXISImage
from iotbx.detectors.raxis_nonsquare import NonSquareRAXISImage
from iotbx.detectors.macscience import DIPImage
from iotbx.detectors.saturn import SaturnImage
from iotbx.detectors.bruker import BrukerImage
from iotbx.detectors.pilatus_minicbf import PilatusImage
from iotbx.detectors.edf import EDFImage
from iotbx.detectors.detectorbase import DetectorImageBase
from iotbx.detectors.pilatus_slice import pilatus_slice_from_file_url
from iotbx.detectors.pilatus_slice import pilatus_slice_from_http_url
from iotbx.detectors.adsc_module import ADSC_module_from_file_url

class EDFWrapper(EDFImage, DetectorImageBase):
  def __init__(self,filename):
    EDFImage.__init__(self,filename)
    self.vendortype = "Pilatus Single Module"
  def readHeader(self):
    EDFImage.readHeader(self)
    self.parameters['PIXEL_SIZE']=0.172
    self.parameters['SIZE1']=self.parameters['Dim_1']
    self.parameters['SIZE2']=self.parameters['Dim_2']
    self.parameters['BEAM_CENTER_X']=0.0 #Dummy argument
    self.parameters['BEAM_CENTER_Y']=0.0 #Dummy argument
    self.parameters['DISTANCE']=100.0 #Dummy argument

all_image_types = [EDFWrapper,SaturnImage,DIPImage,ADSCImage,
                  MARImage,MARIPImage,DTREKImage,RAXISImage,
                  NonSquareRAXISImage,PilatusImage,CBFImage,BrukerImage]

all_url_types = [pilatus_slice_from_file_url,pilatus_slice_from_http_url,
                 ADSC_module_from_file_url
]

names_and_types = { "ADSC"          : ADSCImage,
                    "Saturn"        : SaturnImage,
                    "DIP"           : DIPImage,
                    "MAR"           : MARImage,
                    "MARIP"         : MARIPImage,
                    "DTREK"         : DTREKImage,
                    "RAXIS"         : RAXISImage,
                    "NonSquareRAXIS": NonSquareRAXISImage,
                    "Pilatus"       : PilatusImage,
                    "CBF"           : CBFImage,
                    "Bruker"        : BrukerImage,
                    "EDF"           : EDFWrapper
                   }

def ImageFactory(filename):
  from iotbx.detectors import url_support
  if os.path.isfile(filename):
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
      except Exception:
        pass
  A = url_support.potential_url_request(filename)
  if A.is_url_request():
    for utype in all_url_types:
      try:
        I = utype(filename)
        I.readHeader()
        return I
      except Exception:
        pass
    raise ImageException(filename+" does not work as a functioning image URL")
  raise ImageException(filename+" not recognized as any known detector image type")

def TrySingleImageType(filename,image_type):
  I = names_and_types[ image_type ](filename)
  return I
