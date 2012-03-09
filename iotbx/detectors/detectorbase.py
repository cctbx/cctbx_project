import copy
from iotbx.detectors import ReadADSC
import sys

class DetectorImageBase(object):
  def __init__(self,filename):
    self.filename=filename
    self.parameters=None
    self.linearintdata=None
    self.bin=1
    self.vendortype = "baseclass"
    self.beam_center_reference_frame = "instrument"#cf beam_center_convention.py
    self.beam_center_convention = None

  def copy_common_attributes_from_parent_instance(self, parentobject):
    self.filename = copy.copy(parentobject.filename)
    self.bin = copy.copy(parentobject.bin)
    self.vendortype = copy.copy(parentobject.vendortype)
    self.beam_center_reference_frame = copy.copy(parentobject.beam_center_reference_frame)
    self.beam_center_convention = copy.copy(parentobject.beam_center_convention)
    self.header = copy.copy(parentobject.header)
    self.headerlines = copy.copy(parentobject.headerlines)

  def setBin(self,bin): #software binning.
                        # the only bin values supported are 1 & 2
    if self.bin!=1 or bin!=2: return
    if self.size1%bin!=0: return
    self.parameters['SIZE1']=self.parameters['SIZE1']//bin
    self.parameters['SIZE2']=self.parameters['SIZE2']//bin
    if self.parameters.has_key('CCD_IMAGE_SATURATION'):
      self.parameters['CCD_IMAGE_SATURATION']=self.parameters['CCD_IMAGE_SATURATION']*bin*bin
    self.parameters['PIXEL_SIZE']=self.parameters['PIXEL_SIZE']*bin
    self.bin = bin
    self.bin_safe_set_data(self.linearintdata)

  def set_beam_center_convention(self,beam_center_convention):
    from iotbx.detectors.beam_center_convention import convert_beam_instrument_to_imageblock
    convert_beam_instrument_to_imageblock(self,beam_center_convention)

  def fileLength(self):
    self.readHeader()
    return self.dataoffset()+self.size1*self.size2*self.integerdepth()
    # dataoffset() and integerdepth() must be defined in derived class
    # pure supposition:
    #  size1 corresponds to number of rows.  Columns are slow.
    #  size2 corresponds to number of columns.  Rows are fast.

  def getEndian(self): pass
    # must be defined in derived class

  def endian_swap_required(self):
    data_is_big_endian = self.getEndian()
    import struct
    platform_is_big_endian = (
      struct.unpack('i',struct.pack('>i',3000))[0] == 3000
    )
    return data_is_big_endian != platform_is_big_endian

  def read(self):
    self.fileLength()
    self.bin_safe_set_data(
         ReadADSC(self.filename,self.dataoffset(),
         self.size1*self.bin,self.size2*self.bin,self.getEndian())
         )

  def bin_safe_set_data(self, new_data_array):
    #private interface for software binning 2 X 2.
    #  Any setting of linearintdata must be through this function
    #  self.bin==2: when data are read lazily, they must be binned
    #  new_data_array.bin2by2==True: the data have been binned
    if self.bin==2 and \
       new_data_array != None and\
       new_data_array.__dict__.get("bin2by2")!=True:
      from iotbx.detectors import Bin2_by_2
      self.linearintdata = Bin2_by_2(new_data_array)
      self.linearintdata.bin2by2 = True
    else:
      self.linearintdata = new_data_array

  def get_data_type(self):
    typehash = str(self.linearintdata.__class__)
    if typehash.find("int")>=0: return "int"
    elif typehash.find("double")>=0: return "double"

  def get_flex_image(self,binning=1,brightness=1.0):
    datatype = self.get_data_type()
    if datatype=="int":
      from iotbx.detectors import FlexImage
    elif datatype=="double":
      from iotbx.detectors import FlexImage_d as FlexImage
    return FlexImage(
      rawdata=self.linearintdata,
      binning=binning,
      vendortype=self.vendortype,
      brightness=brightness,
      saturation=int(getattr(self, "saturation", 65535)))

  data_types = dict( SIZE1=int, SIZE2=int, PIXEL_SIZE=float,
                     DISTANCE=float, TWOTHETA=float, OSC_RANGE=float,
                     OSC_START=float, PHI=float, WAVELENGTH=float,
                     BEAM_CENTER_X=float, BEAM_CENTER_Y=float,
                     CCD_IMAGE_SATURATION=int, DETECTOR_SN=str )

  def get_spotfinder(self,distl_params): #following heuristics_base.register_frames() example
    #application-specific adjustments to parameters
    #XXX this should probably be a deep copy of parameters.
    if distl_params.distl.res.inner!=None:
      distl_params.distl_lowres_limit = distl_params.distl.res.inner
    if distl_params.distl.res.outer!=None:
      distl_params.force_method2_resolution_limit = distl_params.distl.res.outer
      distl_params.distl_highres_limit = distl_params.distl.res.outer

    distl_params.distl_force_binning = False
    distl_params.distl_permit_binning = False
    distl_params.wedgelimit = 1
    distl_params.spotfinder_header_tests = False

    #unusual location for min spot area tests...
    from iotbx.detectors.context.config_detector import beam_center_convention_from_image_object
    beam_center_convention_from_image_object(self,distl_params)
    # end special min spot area treatment

    from spotfinder.applications.practical_heuristics import heuristics_base
    from spotfinder.diffraction.imagefiles import file_names
    class empty:pass
    E = empty()
    E.argv = ["Empty",self.filename]
    names = file_names(E)
    this_frame = names.frames()[0]
    process_dictionary = dict(twotheta = "%f"%self.twotheta,
       ybeam = "%f"%self.beamy,
       xbeam = "%f"%self.beamx,
       distance = "%f"%self.distance,
       wavelength = "%f"%self.wavelength,
       template = [item.template for item in names.FN if item.number==this_frame][0],
                              )
    Spotfinder = heuristics_base(process_dictionary,distl_params)
    Spotfinder.images[this_frame] = Spotfinder.oneImage(this_frame,
      Spotfinder.pd, self)
    Spotfinder.determine_maxcell(this_frame,Spotfinder.pd)
    Spotfinder.images[this_frame]['spotoutput']['relpath']=self.filename
    from spotfinder.applications.stats_distl import pretty_image_stats
    pretty_image_stats(Spotfinder,this_frame)
    return Spotfinder,this_frame

  def debug_write(self,fileout,mod_data=None):
    if not self.parameters.has_key("TWOTHETA"):
      self.parameters["TWOTHETA"]=0.0
    if self.getEndian()==1:
      self.parameters["BYTE_ORDER"]="big_endian"
    else:
      self.parameters["BYTE_ORDER"]="little_endian"
    info = """{
HEADER_BYTES= 1024;
DIM=2;
BYTE_ORDER=%(BYTE_ORDER)s;
TYPE=unsigned_short;
SIZE1=%(SIZE1)4d;
SIZE2=%(SIZE1)4d;
PIXEL_SIZE=%(PIXEL_SIZE)8.6f;
TIME=0.000000;
DISTANCE=%(DISTANCE).2f;
TWOTHETA=%(TWOTHETA).2f;
PHI=%(OSC_START).3f;
OSC_START=%(OSC_START).3f;
OSC_RANGE=%(OSC_RANGE).3f;
WAVELENGTH=%(WAVELENGTH).6f;
BEAM_CENTER_X=%(BEAM_CENTER_X).2f;
BEAM_CENTER_Y=%(BEAM_CENTER_Y).2f;
CCD_IMAGE_SATURATION=65535;
}\f"""%self.parameters
    F = open(fileout,"wb")
    F.write(info)
    len_null=1024-len(info)
    F.write('\0'*len_null)
    F.close()
    from iotbx.detectors import WriteADSC
    if mod_data==None: mod_data=self.linearintdata
    WriteADSC(fileout,mod_data,self.size1,self.size2,self.getEndian())

  def __getattr__(self, attr):
    if   attr=='size1' : return self.parameters['SIZE1']
    elif attr=='size2' : return self.parameters['SIZE2']
    elif attr=='npixels' : return self.parameters['SIZE1'] * self.parameters['SIZE2']
    elif attr=='saturation' : return self.parameters.get('CCD_IMAGE_SATURATION',65535)
    elif attr=='rawdata' : return self.linearintdata
    elif attr=='pixel_size' : return self.parameters['PIXEL_SIZE']
    elif attr=='osc_start' : return self.parameters['OSC_START']
    elif attr=='distance' : return self.parameters['DISTANCE']
    elif attr=='wavelength' : return self.parameters['WAVELENGTH']
    elif attr=='beamx' : return self.parameters['BEAM_CENTER_X']
    elif attr=='beamy' : return self.parameters['BEAM_CENTER_Y']
    elif attr=='deltaphi' : return self.parameters['OSC_RANGE']
    elif attr=='twotheta' : return self.parameters.get('TWOTHETA',0.0)
    elif attr=='serial_number' : return self.parameters['DETECTOR_SN']

  def show_header(self, out=None):
    if (out is None) :
      out = sys.stdout
    print >> out, "File:",self.filename
    print >> out, "Number of pixels: slow=%d fast=%d"%(self.size1,self.size2)
    print >> out, "Pixel size: %f mm"%self.pixel_size
    print >> out, "Saturation: %.0f"%self.saturation
    print >> out, "Detector distance: %.2f mm"%self.distance
    print >> out, "Detector 2theta swing: %.2f deg."%self.twotheta
    print >> out, "Rotation start: %.2f deg."%self.osc_start
    print >> out, "Rotation width: %.2f deg."%self.deltaphi
    print >> out, "Beam center: x=%.2f mm  y=%.2f mm"%(self.beamx,self.beamy)
    print >> out, "Wavelength: %f Ang."%self.wavelength

  # code developed for the image viewer. phil_parameters is a scope extract
  def initialize_viewer_properties(self,phil_parameters):

    self._invert_beam_center = False
    from iotbx.detectors.context.config_detector import \
      beam_center_convention_from_image_object
    bc = beam_center_convention_from_image_object(self,phil_parameters)
    print "beam center convention: %d" % bc
    # FIXME what about 2-4 & 6-7?
    if (bc == 0) :
      self._invert_beam_center = True
      self._invert_y = True
    elif (bc == 1) :
      self._invert_y = False
    elif (bc == 5) :
      self._invert_y = True

    self.image_size_fast = self.size2 # width
    self.image_size_slow = self.size1 # height
    self.pixel_resolution = self.pixel_size

  def detector_coords_as_image_coords_float (self, x, y) :
    """
    Convert absolute detector position (in mm) to floating-value image pixel coordinates.
    """
    dw = self.image_size_fast * self.pixel_resolution
    dh = self.image_size_slow * self.pixel_resolution
    x_frac = x / dw
    if (self._invert_y) :
      y_frac = - ((y / dh) - 1.0)
    else :
      y_frac = y / dh
    return x_frac * self.image_size_fast, \
           y_frac * self.image_size_slow

  def detector_coords_as_image_coords (self, x, y) :
    """
    Convert absolute detector position (in mm) to integer-value image pixel coordinates.
    """
    x_point,y_point = self.detector_coords_as_image_coords_float(x,y)
    return (int(x_point), int(y_point))

  def image_coords_as_detector_coords (self, x, y) :
    """
    Convert image pixel coordinates to absolute position on the detector
    (in mm).
    """
    dw = self.image_size_fast * self.pixel_resolution
    dh = self.image_size_slow * self.pixel_resolution
    x_frac = x / self.image_size_fast
    y_frac = y / self.image_size_slow
    x_detector = x_frac * dw
    if (self._invert_y) :
      y_detector = (1.0 - y_frac) * dh
    else :
      y_detector = y_frac * dh
    return x_detector, y_detector

  def get_beam_center_mm (self) :
    # FIXME Pilatus and ADSC images appear to have different conventions???
    if (self._invert_beam_center) :
      center_x = self.beamy
      center_y = self.beamx
    else :
      center_x = self.beamx
      center_y = self.beamy
    return center_x, center_y

  def get_beam_center_pixels_fast_slow(self):
    center_x, center_y = self.get_beam_center_mm()
    return self.detector_coords_as_image_coords_float(center_x, center_y)

  def get_pixel_intensity(self,coords):
    try:
      return self.linearintdata[(int(coords[0]), int(coords[1]))]
    except IndexError:
      return None
