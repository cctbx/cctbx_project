import os
import re

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

def identify_dataset (path_name) :
  def get_file_name_components (file_name_) :
    base, ext = os.path.splitext(os.path.basename(file_name_))
    fields = base.split("_")
    suffix = fields[-1]
    common_base = "_".join(fields[:-1])
    return (common_base, suffix, ext)
  suffix = common_base = common_ext = None
  file_name = dir_name = None
  if (os.path.isfile(path_name)) :
    file_name = os.path.abspath(path_name)
  elif (os.path.isdir(path_name)) :
    dir_name = os.path.abspath(path_name)
  else :
    assert 0
  if (file_name is not None) :
    dir_name = os.path.dirname(file_name)
    (common_base, suffix, common_ext) = get_file_name_components(file_name)
  all_files = os.listdir(dir_name)
  stacks = {}
  suffixes = {}
  extensions = {}
  for fn in all_files :
    if (common_base is not None) :
      if (fn.startswith(common_base)) and (fn.endswith(common_ext)) :
        (base2, suffix2, ext2) = get_file_name_components(fn)
        if (common_base in stacks) :
          stacks[common_base].append(int(suffix2))
        else :
          suffixes[common_base] = "#" * len(suffix)
          stacks[common_base] = [ int(suffix2) ]
          extensions[common_base] = ext2
    else :
      (base2, suffix2, ext2) = get_file_name_components(fn)
      # FIXME probably not a comprehensive list...
      if (ext2 in [".img",".osc",".ccd",".mccd",".cbf"]) :
        if (base2 in stacks) :
          stacks[base2].append(int(suffix2))
        else :
          stacks[base2] = [ int(suffix2) ]
          suffixes[base2] = "#" * len(suffix2)
          extensions[base2] = ext2
  results = []
  for base in sorted(stacks.keys()) :
    ranges = []
    img_start = img_last = None
    for x in sorted(stacks[base]) :
      if (img_start is None) :
        img_start = x
      elif (img_last is not None) and (x > (img_last + 1)) :
        ranges.append((img_start, img_last))
        img_start = x
      img_last = x
    ranges.append((img_start, img_last))
    file_base = "%s_%s%s" % (base, suffixes[base], extensions[base])
    #print "%s %s" % (file_base,
    #  ", ".join([ "%d-%d" % (a,b) for (a, b) in ranges ]))
    dataset = dataset_info(
      base_name=os.path.join(dir_name, file_base),
      ranges=ranges)
    results.append(dataset)
  return results

class dataset_info (object) :
  def __init__ (self, base_name, ranges) :
    self.base_name = base_name
    self.ranges = ranges

  def format (self) :
    ranges_strs = []
    for start, end in self.ranges :
      if (start == end) :
        ranges_strs.append(str(start))
      else :
        ranges_strs.append("%d-%d" % (start, end))
    return "%s (%s)" % (self.base_name, ",".join(ranges_strs))

  def __str__ (self) :
    return self.format()

  def get_frame_path (self, frame) :
    assert isinstance(frame, int)
    serial_format = "%%0%dd" % (self.base_name.count("#"))
    format_str = re.sub("[#]{1,}", serial_format, self.base_name)
    return format_str % frame
