import cStringIO as StringIO
from iotbx.detectors.adsc_module import ADSCModule
from iotbx.detectors.pilatus_slice import PilatusSlice

class attributes_from_dict:
  def __init__(self,parentdict):
    for key in parentdict.keys():
      self.__dict__[key]=parentdict[key]

def module_or_slice_from_http_request(request):
  # edit out data from list
  tags = {}
  for item in request.keys():
    assert len(request[item])==1 #list with one element
    tags[item]=request[item][0]
    if len(tags[item]) < 1000:  #short items could be converted to numbers
      try:
        tags[item] = float(tags[item])
        if tags[item]==round(tags[item],0):
          tags[item]= int(tags[item])
      except Exception: pass
      #print "posted data types",item, type(tags[item]), tags[item]
  # tags dictionary is now conditioned as to type
  tags_class = attributes_from_dict(tags)

  #Vendor selector.
  if tags["vendortype"].lower().find("adsc")>=0:
    P = ADSCModule()
  elif tags["vendortype"].lower().find("pilatus")>=0:
    P = PilatusSlice()

  P.moduleindex = tags["moduleindex"]
  P.object = StringIO.StringIO(tags["adsc_data"])
  P.copy_common_attributes_from_parent_instance(tags_class)
  P.parameters = {}
  for item in ['DISTANCE', 'PHI', 'WAVELENGTH',
    'TWOTHETA', 'OSC_RANGE',
    'CCD_IMAGE_SATURATION', 'OSC_START', 'DETECTOR_SN', 'PIXEL_SIZE',
    'SIZE1','SIZE2','BEAM_CENTER_X','BEAM_CENTER_Y'
    ]:
    P.parameters[item] = tags[item]
  from iotbx.detectors.beam_center_convention import convert_beam_instrument_to_imageblock
  convert_beam_instrument_to_imageblock(P,P.beam_center_convention)
  P.slice_callback = P.slice_callback_with_high_performance_http_data
  return P
