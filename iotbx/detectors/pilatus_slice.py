import copy,cPickle
from iotbx.detectors.pilatus_minicbf import PilatusImage
from iotbx.detectors import ImageException
from scitbx.array_family import flex

def pilatus_slice_from_http_url(url):
  #backward compatibility with Python 2.5
  try: from urlparse import parse_qs
  except Exception: from cgi import parse_qs

  from urlparse import urlparse
  parsed = urlparse(url)
  assert parsed.scheme in ["http","https"]
  from urllib2 import urlopen
  Response = urlopen(url)
  info = Response.info()
  #print info
  P = PilatusSlice()
  if "Image-slice" in info.keys():
    sliceindex = int(info["Image-slice"])
    assert 0 <= sliceindex <= 11 # slow slice section, possible index range
    P.sliceindex = sliceindex
  P.bin = int(info["Image-bin"])
  P.filename = info["Image-filename"]
  P.vendortype = "Pilatus-6M"
  P.parameters = {}
  for item in P.mandatory_keys:
      header_key = "Image-%s"%item.lower()
      P.parameters[item]=P.data_types[item](info[header_key])
  if info["Image-data_encryption"]=="pickle":
    P.slice_callback = P.slice_callback_with_portable_http_data
  elif info["Image-data_encryption"]=="byte":
    P.slice_callback = P.slice_callback_with_high_performance_http_data
  P.object = Response #hand over the data stream to the callback function
  return P

def pilatus_slice_from_file_url(url):
  #backward compatibility with Python 2.5
  try: from urlparse import parse_qs
  except Exception: from cgi import parse_qs

  from urlparse import urlparse
  parsed = urlparse(url)
  assert parsed.scheme == "file"
  file = parsed.path.split("?")[0]
  if file == parsed.path:
    return PilatusImage(file)
  qs = parse_qs(parsed.path.split("?")[1])
  sliceindex = int(qs["slice"][0])
  object = PilatusImage(file)
  object.readHeader()
  return pilatus_slice_from_object_and_slicenumber(object,sliceindex)

def pilatus_slice_from_object_and_slicenumber(object,sliceindex):
  P = PilatusSlice()
  assert 0 <= sliceindex <= 11 # slow slice section, possible index range
  P.sliceindex = sliceindex
  P.object = object
  P.copy_common_attributes_from_parent_instance(P.object)
  P.vendortype = "Pilatus-6M" # overrides the default copy
  P.parameters = P.slice_parameters(object.parameters)
  P.slice_callback = P.slice_callback_with_object_data
  return P

class PilatusSlice(PilatusImage):
  def __init__(self):
    self.already_read_data = False

  data_types = copy.copy(PilatusImage.data_types)
  data_types.update({"PIXEL_SIZE_UNITS":str, "DISTANCE_UNITS":str, "AXIS":str})

  def slice_parameters(self,param):
    #unchanged parameters first:
    result = {}
    for item in ['PIXEL_SIZE_UNITS', 'DISTANCE', 'PHI', 'WAVELENGTH',
    'SIZE2', 'TWOTHETA', 'DISTANCE_UNITS', 'OSC_RANGE', 'BEAM_CENTER_Y',
    'CCD_IMAGE_SATURATION', 'OSC_START', 'DETECTOR_SN', 'PIXEL_SIZE',
    'AXIS']:
      result[item]=copy.copy(param[item])
    result['SIZE1']=195
    #convert beam center X (slow) back to pixels, correct for slice, back to mm
    pixels = param["BEAM_CENTER_X"] / param["PIXEL_SIZE"]
    correction = pixels - self.sliceindex * (195+17)
    result["BEAM_CENTER_X"] = correction * result["PIXEL_SIZE"]
    return result

  def slice_callback_with_portable_http_data(self):
    linearintdata = cPickle.load(self.object)
    del self.object #once the data are copied, close the stream
    return linearintdata

  def slice_callback_with_high_performance_http_data(self):
    BYU = self.object.read()
    linearintdata = flex.int_from_byte_str(BYU)
    provisional_size = linearintdata.size()
    if provisional_size==480285:
      linearintdata.reshape(flex.grid((195,2463)))
    elif provisional_size==6224001:
      linearintdata.reshape(flex.grid((2527,2463)))
    else:
      raise ImageException("wrong number of pixels for Pilatus image")
    del self.object #once the data are copied, close the stream
    return linearintdata

  def slice_callback_with_object_data(self):
    self.object.read()
    start_index = 2463*((195+17)*self.sliceindex)
    stop_index = start_index + 2463*195
    linearintdata = flex.int(self.object.linearintdata[start_index:stop_index])
    linearintdata.reshape(flex.grid((195,2463)))
    del self.object #once the data are copied, no need to keep the original
    return linearintdata

  def read(self):
    if self.already_read_data: return
    self.bin_safe_set_data( self.slice_callback() )
    self.already_read_data = True

if __name__=='__main__':
  import sys
  full_path_to_file = sys.argv[1]
  a = PilatusImage(testing_file)
  a.read()
  print a
  print a.parameters
  print a.rawdata, len(a.rawdata), a.size1*a.size2
  for dataitem in ['bin', 'filename', 'header', 'headerlines', 'linearintdata', 'parameters', 'vendortype']:
    print dataitem,
    exec("print a.%s"%dataitem)
  print pilatus_slice_from_object_and_slicenumber(a,5)

  P = pilatus_slice_from_file_url(url="file://%s?slice=5"%full_path_to_file)
