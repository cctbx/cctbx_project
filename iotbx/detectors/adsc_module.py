import copy
from iotbx.detectors.adsc import ADSCImage
from scitbx.array_family import flex
from iotbx.detectors import image_divider

def vendor_specific_null_value(object):
    if object.vendortype.lower().find("pilatus")>=0:
      nullvalue = -1
    else:
      nullvalue = 0
    return nullvalue

def ADSC_module_from_file_url(url):
  #backward compatibility with Python 2.5
  try: from urlparse import parse_qs
  except Exception: from cgi import parse_qs

  from urlparse import urlparse
  parsed = urlparse(url)
  assert parsed.scheme == "file"
  file = parsed.path.split("?")[0]
  if file == parsed.path:
    return ADSCImage(file)
  qs = parse_qs(parsed.path.split("?")[1])
  sliceindex = int(qs["slice"][0])
  object = ADSCImage(file)
  object.readHeader()
  return ADSC_module_from_object_and_slicenumber(object,sliceindex)

def ADSC_module_from_object_and_slicenumber(object,moduleindex):
  P = ADSCModule()
  P.moduleindex = moduleindex
  P.object = object
  P.copy_common_attributes_from_parent_instance(P.object)
  P.vendortype = "ADSC module"
  P.parameters = P.module_parameters(P.object)
  P.slice_callback = P.slice_callback_with_object_data
  return P

class ADSCModule(ADSCImage):
  def __init__(self):
    self.already_read_data = False

  def module_parameters(self,object):
    param = object.parameters
    #unchanged parameters first:
    result = {}
    for item in ['DISTANCE', 'PHI', 'WAVELENGTH',
    'TWOTHETA', 'OSC_RANGE',
    'CCD_IMAGE_SATURATION', 'OSC_START', 'DETECTOR_SN', 'PIXEL_SIZE',
    ]:
      result[item]=copy.copy(param[item])

    #other parameters keep a record of parent settings; but child
    # settings are deduced on-the-fly once the data are read.
    for item in ['SIZE1','SIZE2','BEAM_CENTER_X','BEAM_CENTER_Y',
    ]:
      result["PARENT_"+item]=copy.copy(param[item])

    #the final parameters require a knowledge of the module boundaries.
    # In this implementation, these boundaries are determined on the fly
    # from the raw data, implying an up-front file read.
    # The only way to avoid this (not implemented here) is to encode
    # the vendor- & model-specific module boundaries.

    object.read()
    nullvalue = vendor_specific_null_value(object)
    ID = image_divider(data = object.linearintdata, nullvalue=nullvalue)
    print "module slow interval",ID.tile_slow_interval(self.moduleindex).first, ID.tile_slow_interval(self.moduleindex).last
    print "module fast interval",ID.tile_fast_interval(self.moduleindex).first, ID.tile_fast_interval(self.moduleindex).last
    result["SIZE1"] = ID.tile_slow_interval(self.moduleindex).size()
    result["SIZE2"] = ID.tile_fast_interval(self.moduleindex).size()
    print "size1",result["SIZE1"],"size2",result["SIZE2"]
    return result

  def set_beam_center_convention(self,beam_center_convention):

    #previously part of the module_parameters function
    #from iotbx.detectors.context.config_detector import beam_center_convention_from_image_object
    #beam_center_convention = beam_center_convention_from_image_object(object)
    print "CC",beam_center_convention
    assert self.object.beam_center_reference_frame == "instrument"
    nullvalue = vendor_specific_null_value(self.object)
    ID = image_divider(data = self.object.linearintdata, nullvalue=nullvalue)
    from iotbx.detectors.beam_center_convention import convert_beam_instrument_to_module
    self.parameters['BEAM_CENTER_X'],self.parameters['BEAM_CENTER_Y'] = convert_beam_instrument_to_module(
      self.object, ID, self.moduleindex, beam_center_convention)
    self.beam_center_reference_frame = "imageblock"
    self.beam_center_convention = beam_center_convention
    print "old beam center",self.object.parameters['BEAM_CENTER_X'],self.object.parameters['BEAM_CENTER_Y']
    print "for module",   self.moduleindex,"beam x,y is", self.parameters['BEAM_CENTER_X'],self.parameters['BEAM_CENTER_Y']


  def slice_callback_with_high_performance_http_data(self):
    BYU = self.object.read()
    linearintdata = flex.int_from_byte_str(BYU)
    provisional_size = linearintdata.size()
    assert provisional_size == self.size1*self.size2
    print "size assertion OK",provisional_size
    linearintdata.reshape(flex.grid((self.size1,self.size2)))
    del self.object #once the data are copied, close the stream
    return linearintdata

  def slice_callback_with_object_data(self):
    self.object.read()
    nullvalue = vendor_specific_null_value(self.object)
    ID = image_divider(data = self.object.linearintdata, nullvalue=nullvalue)
    assert 0 <= self.moduleindex < ID.module_count()
    new_data_array = ID.tile_data(self.moduleindex)
    if self.object.linearintdata.__dict__.get("bin2by2")==True:
      new_data_array.bin2by2=True
    del self.object #once the data are copied, no need to keep the original
    return new_data_array

  def read(self):
    if self.already_read_data: return
    self.bin_safe_set_data(self.slice_callback())
    self.already_read_data = True

if __name__=='__main__':
  import sys
  sliceidx=4
  full_path_to_file = sys.argv[1]
  a = ADSCImage(full_path_to_file)
  a.read()
  print a
  print a.parameters
  print a.rawdata, len(a.rawdata), a.size1*a.size2
  for dataitem in ['bin', 'filename', 'header', 'headerlines', 'linearintdata', 'parameters', 'vendortype']:
    print dataitem,
    exec("print a.%s"%dataitem)
  print ADSC_module_from_object_and_slicenumber(a,sliceidx)

  P = ADSC_module_from_file_url(url="file://%s?slice=%d"%(full_path_to_file,sliceidx))
  print "file://%s?slice=%d"%(full_path_to_file,sliceidx)
  print P
