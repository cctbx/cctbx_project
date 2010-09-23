import re,types,struct
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from iotbx.detectors.detectorbase import DetectorImageBase
from iotbx.detectors import ReadDTrek

verbose=False

class DTREKImage(DetectorImageBase):
  """enforce dTREK Image Format v1.1, using Rigaku/MSC documentation."""

  def __init__(self,filename):
    DetectorImageBase.__init__(self,filename)
    self.vendortype = "RAXIS"

  def read_vendor_header(self):
    G = open(self.filename, "rb")
    assert G.read(14)=="{\nHEADER_BYTES"
    raw = [charac for charac in G.read(6) if charac.isdigit()]
    header_bytes = int("".join(raw))
    assert header_bytes%512==0
    G.seek(0)
    padded_header = G.read(header_bytes)
    G.close()
    unpadded_header = padded_header.rstrip()
    assert unpadded_header[0:2]=="{\n"
    assert unpadded_header[-2:]=="\n}"
    self.header = padded_header
    raw_key_values = unpadded_header[2:-3].split(";\n")
    self.headerlines = []; self.keys={}
    for i in raw_key_values:
      self.headerlines.append( i.split("=") )
      self.keys[self.headerlines[-1][0]]=self.headerlines[-1][1]

  def enforce_types(self):
    # mandate:=( regex_expression, python_type, length )
    # length:= 1: single variable
    #          0: unknown, could be multiple
    #         >1: definite multiple value
    self.enf= [("SIZE1",int,1),
               ("SIZE2",int,1),
               ("SATURATED_VALUE",int,1),
               ("DETECTOR_NAMES",types.StringType,0),
               (r"\n([0-9A-Za-z]+_DETECTOR_DIMENSIONS)",float,2),
               (r"\n([0-9A-Za-z]+_DETECTOR_SIZE)",float,2),
               (r"\n(ROTATION)=",float,10),
               (r"\n([0-9A-Za-z]+_GONIO_NAMES)",types.StringType,0),
               (r"\n([0-9A-Za-z]+_GONIO_UNITS)",types.StringType,0),
               (r"\n([0-9A-Za-z]+_GONIO_VALUES)",float,0),
               ("SOURCE_WAVELENGTH",float,2), #must be one wavelength
               (r"\n([0-9A-Za-z]+_SPATIAL_DISTORTION_TYPE)",types.StringType,1),
               (r"\n([0-9A-Za-z]+_SPATIAL_DISTORTION_INFO)",float,4),
               ("DETECTOR_TYPE",types.StringType,1),
               ("DATA_TYPE",types.StringType,1),
               ("BYTE_ORDER",types.StringType,1),
               ("RAXIS_COMPRESSION_RATIO",int,1),
               ("HEADER_BYTES",int,1),
              ]
    for mandate in self.enf:
      pattern = re.compile(mandate[0])
      matches = pattern.findall(self.header)
      for match in matches:
        if verbose: print match,
        if mandate[2]==1:
          self.keys[match] = mandate[1](self.keys[match])
        else:
          all_tokens = self.tokenize(self.keys[match])
          all_values = [mandate[1](az) for az in all_tokens]
          if mandate[2] > 1: assert len(all_values)==mandate[2]
          self.keys[match] = all_values
        if verbose:
          print self.keys[match],
          print

    for integer in ["SIZE1","SIZE2","SATURATED_VALUE"]:
      self.keys[integer]=int(self.keys[integer])

  def tokenize(self,string_):
    tokens = string_.split(" ")
    while "" in tokens: tokens.remove("")
    return tokens

  def generic_param_from_vendor_head(self):
      self.parameters={}
      if verbose:
       for i in self.headerlines:
        print "%29s"%i[0],self.keys[i[0]]
      # Note that SIZE1 is slow for ADSC/CBF but fast for RAXIS
      # Note that SIZE2 is fast for ADSC/CBF but slow for RAXIS
      self.parameters['SIZE1'] = self.keys["SIZE1"]
      self.parameters['SIZE2'] = self.keys["SIZE2"]
      self.parameters['CCD_IMAGE_SATURATION'] = self.keys["SATURATED_VALUE"]

      dname_prefix = self.keys["DETECTOR_NAMES"][0]

      sizes = flex.double(self.keys[dname_prefix+"DETECTOR_SIZE"])
      pixels = flex.double(self.keys[dname_prefix+"DETECTOR_DIMENSIONS"])
      pixel_sizes = sizes/pixels
      assert pixel_sizes[0]==pixel_sizes[1]
      self.parameters['PIXEL_SIZE'] = pixel_sizes[0]

      assert approx_equal(self.keys["ROTATION"][1]-self.keys["ROTATION"][0],
                          self.keys["ROTATION"][2])
      self.parameters['OSC_START'] = self.keys["ROTATION"][0]
      self.parameters['OSC_RANGE'] = self.keys["ROTATION"][2]

      distance_idx = self.keys[dname_prefix+"GONIO_NAMES"].index("Distance")
      assert self.keys[dname_prefix+"GONIO_UNITS"][distance_idx]=="mm"
      self.parameters['DISTANCE'] = self.keys[dname_prefix+"GONIO_VALUES"][distance_idx]

      assert self.keys["SOURCE_WAVELENGTH"][0]==1.0
      self.parameters['WAVELENGTH'] = self.keys["SOURCE_WAVELENGTH"][1]

      assert self.keys[dname_prefix+"SPATIAL_DISTORTION_TYPE"]=="Simple_spatial"
      check_pixel_sizes = flex.double(
        self.keys[dname_prefix+"SPATIAL_DISTORTION_INFO"][2:4])
      assert check_pixel_sizes == pixel_sizes
      beam_mm = flex.double(
        self.keys[dname_prefix+"SPATIAL_DISTORTION_INFO"][0:2])*pixel_sizes

      self.parameters['BEAM_CENTER_X'] = beam_mm[0]
      self.parameters['BEAM_CENTER_Y'] = beam_mm[1]
      tt_idx = self.keys[dname_prefix+"GONIO_NAMES"].index("2Theta")
      assert self.keys[dname_prefix+"GONIO_UNITS"][tt_idx]=="deg"
      self.parameters['TWOTHETA'] = self.keys[dname_prefix+"GONIO_VALUES"][tt_idx]
      self.parameters['DETECTOR_SN'] = self.keys["DETECTOR_TYPE"]


  def readHeader(self,):
    if not self.parameters:
      self.read_vendor_header()
      self.enforce_types()
      self.generic_param_from_vendor_head()

  def getEndian(self):
    self.readHeader()
    if self.keys['BYTE_ORDER'].lower().find('big')>=0:
      return 1 #big_endian
    else:
      return 0 #little_endian

  def read(self):
    G = open(self.filename, "rb")
    G.seek(self.keys["HEADER_BYTES"])

    endian_code = {'little_endian':'<','big_endian':'>'}[self.keys["BYTE_ORDER"]]
    type_code = {'signed char':'b',
                 'unsigned char':'B',
                 'short int':'h',
                 'long int':'i',
                 'unsigned short int':'H',
                 'unsigned long int':'I',
                 'float IEEE':'f',
                }[self.keys['Data_type']]
    type_size = {'b':1,'B':1,'h':2,'H':2,'i':4,'I':4,'f':4}[type_code]
    assert not type_code=="I" # for I, a flex.int() will exceed type limits
    array_size = self.parameters['SIZE1'] * self.parameters['SIZE2']
    rawdata = G.read(array_size * type_size)
    G.close()
      #Python prototype--
      #doesn't handle raxis uncompression & is 10x slower than C++ version
      #uncoded_data = struct.unpack(endian_code+type_code*array_size,rawdata)
      #provisional_data = flex.int(uncoded_data)
      #provisional_data.reshape(flex.grid((self.parameters['SIZE2'],
      #                                  self.parameters['SIZE1'])))
      #self.bin_safe_set_data(provisional_data)
    self.bin_safe_set_data(
                ReadDTrek(raw=rawdata,type_code=type_code,
                          slow=self.parameters['SIZE2'],
                          fast=self.parameters['SIZE1'],
                          swap=self.endian_swap_required(),
                          uncompress=self.keys.get("RAXIS_COMPRESSION_RATIO",1)
                          ))
