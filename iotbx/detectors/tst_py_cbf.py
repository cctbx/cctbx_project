import iotbx.cif
import sys, email.parser, copy, hashlib, base64
from cbflib_adaptbx import uncompress,assert_equal
from iotbx.detectors.detectorbase import DetectorImageBase

class cif_binary_section:
  endline = "\r\n"
  pattern = "--CIF-BINARY-FORMAT-SECTION--"
  init_boundary = endline + pattern + endline
  final_boundary = endline + pattern + "--" + endline
  binary_string_separator=endline + endline
  cbf_signature = chr(0x0C)+chr(0x1A)+chr(0x04)+chr(0xD5)

  def from_compressed_string(self,data,init,final):
    section_start = init + len(self.init_boundary)
    section_stop = final
    #offsets for section data, numbering is on the original file string
    self.init = init
    self.final = final

    divide = data.find(self.binary_string_separator,section_start,section_stop)
    assert divide > 0 # binary section consists of a header and data
    header = data[section_start:divide]

    self.header_dic = email.parser.Parser().parsestr(header)
    assert self.header_dic["Content-Type"].find("x-CBF_BYTE_OFFSET")>0

    bin_start = divide + len(self.binary_string_separator)
    assert data[bin_start:bin_start+4]==self.cbf_signature
    self.data = data[bin_start+4:section_stop]
    self.data_type = "compressed"

    m = hashlib.md5()
    m.update(self.data)
    derived_digest = base64.b64encode(m.digest())
    assert self.header_dic["Content-MD5"] == derived_digest

    self.size_fast = int(self.header_dic["X-Binary-Size-Fastest-Dimension"])
    self.size_slow = int(self.header_dic["X-Binary-Size-Second-Dimension"])
    total_elements = int(self.header_dic["X-Binary-Number-of-Elements"])
    compressed_size = int(self.header_dic["X-Binary-Size"])
    assert total_elements==self.size_fast*self.size_slow
    assert compressed_size==len(self.data)

    return self

  def uncompress_in_place(self):
    if self.data_type == "compressed":
      decompressed_data = uncompress(packed=self.data, fast=self.size_fast, slow=self.size_slow)
      self.data = decompressed_data
      self.data_type = "uncompressed"
    assert self.data_type=="uncompressed"
    return self.data

  def show(self):
    print self.header_dic

def get_binary_sections(raw):
  return_sections = []
  ptr = 0
  end = len(raw)
  oldinit = 0
  oldfinal = 0

  while ptr < end:
    init = raw.find(cif_binary_section.init_boundary, oldinit+1)
    if init==-1: break
    assert init > oldfinal

    final = raw.find(cif_binary_section.final_boundary, init)
    assert final > init

    return_sections.append(cif_binary_section().from_compressed_string(raw,init,final))
    oldinit = copy.copy(init)
    oldfinal = copy.copy(final)
    ptr = copy.copy(oldfinal)

  return return_sections

def get_header_sections(raw,binary_sections):
  provisional_slices = [(0,len(raw))]
  for i,section in enumerate(binary_sections):
    first = (provisional_slices[i][0],section.init)
    last = (section.final+len(cif_binary_section.final_boundary),provisional_slices[i][1])
    provisional_slices[i]=first
    provisional_slices.append(last)
  return cif_binary_section.endline.join([raw[s[0]:s[1]] for s in provisional_slices])

class Goniometer:
  def __init__(self,model):
    #problem to report back to Richard Gildea; not possible to print model["_diffrn_measurement"]
    axis_id = model["_diffrn_measurement_axis.axis_id"]

    #get the PHI axis setting
    row_idx = [str(a) for a in model["_diffrn_scan_frame_axis.axis_id"]].index(axis_id)
    self.osc_start = float(model["_diffrn_scan_frame_axis.angle"][row_idx])

    row_idx = [str(a) for a in model["_diffrn_scan_axis.axis_id"]].index(axis_id)
    self.osc_range = float(model["_diffrn_scan_axis.angle_increment"][row_idx])

def get_ad_hoc_beam(model):
    row_idx = [str(a) for a in model["_axis.id"]].index("ELEMENT_X")
    beamx = float(model["_axis.offset[2]"][row_idx])
    beamy = -float(model["_axis.offset[1]"][row_idx])
    return beamx,beamy

from iotbx.detectors.cbf import CBFImage
class pyCBFImage(CBFImage):
  def __init__(self, file_name):
    DetectorImageBase.__init__(self, file_name)
    raw = open(file_name, "rb").read()
    self.binary_sections = get_binary_sections(raw)
    self.header_sections = get_header_sections(raw, self.binary_sections)

    assert len(self.binary_sections)==1

    cif = iotbx.cif.fast_reader(input_string=self.header_sections)
    self.cif_model = cif.model()
    im1 = self.cif_model["image_1"]

    self.vendortype = "CBF"
    self.readHeader()

  def readHeader(self):
    model = self.cif_model["image_1"]
    binaries = self.binary_sections
    goniometer = Goniometer(model)
    ad_hoc_beam = get_ad_hoc_beam(model)

    self.parameters = {'SIZE1': binaries[0].size_fast, #not sure about precedence; implement later
                       'SIZE2': binaries[0].size_slow,
                       'CCD_IMAGE_SATURATION':int(model["_array_intensities.overload"]),
                       'PIXEL_SIZE':self.cbf_simple_py_get_pixel_size(model),
                       'OSC_START':goniometer.osc_start,
                       'DISTANCE':float(model["_diffrn_measurement.sample_detector_distance"]),
                       'WAVELENGTH':self.cbf_simple_py_get_wavelength(model),

                       #Ad-hoc implementation to get a quick beam center.  This is valid
                       #only for winter/diamond dataset; a general implementation will be needed.
                       'BEAM_CENTER_X':ad_hoc_beam[0],
                       'BEAM_CENTER_Y':ad_hoc_beam[1],
                       'OSC_RANGE':goniometer.osc_range,
                       'TWOTHETA':0.0,  #non-zero twotheta not supported (yet)!
                       'DETECTOR_SN':0
                       }
    self.binaries = binaries

  def read(self):
    assert len(self.binaries)==1
    self.linearintdata = self.binaries[0].uncompress_in_place()

  def cbf_simple_py_get_pixel_size(self,model):
    element_number = 0 # used in C version; not implemented in py version (yet)
    axis_number = 1
    array_id = model["_diffrn_data_frame.array_id"]

    #Given the input axis_number, look in the array_structure_list
    #  Table to find the index of the axis whose precedence==axis_number,
    #  for the correct array_id
    array_mask = model["_array_structure_list.array_id"]==array_id
    precedence = [int(a) for a in model["_array_structure_list.precedence"]]
    for i in xrange(len(array_mask)):
      if not array_mask[i]: precedence[i]=0
    idx = precedence.index(axis_number)
    axis_index = int( model["_array_structure_list.index"][idx] )
    assert axis_index > 0

    #Now find the array element size for the given axis_index
    array_mask = model["_array_element_size.array_id"]==array_id
    index_array = [int(a) for a in model["_array_element_size.index"]]
    for i in xrange(len(array_mask)):
      if not array_mask[i]: index_array[i]=0
    size_index = index_array.index(axis_index)

    pixel_size = 1000. * float(model["_array_element_size.size"][size_index])
    return pixel_size

  def cbf_simple_py_get_wavelength(self,model):
    return float(model["_diffrn_radiation_wavelength.wavelength"])

def run(file_name):
  from libtbx.test_utils import approx_equal
  py_image_obj = pyCBFImage(file_name)
  py_image_obj.read()
  c_image_obj = CBFImage(file_name)
  c_image_obj.read()
  assert_equal(py_image_obj.linearintdata, c_image_obj.linearintdata)
  assert py_image_obj.size1 == c_image_obj.size1
  assert py_image_obj.size2 == c_image_obj.size2
  assert approx_equal(py_image_obj.saturation, c_image_obj.saturation)
  assert approx_equal(py_image_obj.pixel_size, c_image_obj.pixel_size)
  assert approx_equal(py_image_obj.osc_start, c_image_obj.osc_start)
  assert approx_equal(py_image_obj.deltaphi, c_image_obj.deltaphi)
  assert approx_equal(py_image_obj.wavelength, c_image_obj.wavelength)
  assert approx_equal(py_image_obj.distance, c_image_obj.distance)
  assert approx_equal(py_image_obj.beamx, c_image_obj.beamx)
  assert approx_equal(py_image_obj.beamy, c_image_obj.beamy)

if (__name__ == "__main__"):
  args = sys.argv[1:]
  for file_name in args:
    #print file_name
    run(file_name)
  print "OK"
