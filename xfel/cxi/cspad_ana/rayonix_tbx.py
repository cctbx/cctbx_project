from __future__ import absolute_import, division, print_function
# Utility functions for the Rayonix Detector.
from scitbx.array_family import flex

# given value of rayonix detector saturation xppi6115
rayonix_saturated_value = 2**16 -1

# minimum value for rayonix data (actually this number is not trusted, numbers above it are)
rayonix_min_trusted_value = -1

def get_rayonix_pixel_size(bin_size):
  ''' Given a bin size determine a pixel size.

Michael Blum from Rayonix said The pixel size is recorded in the header,
but can be derived trivially from the overall dimension of the corrected imaging
area (170mm) and the number of pixels. (3840 unbinned). The corrected image is
forced to this size.

unbinned 170/3840  = 0.04427

I believe the accuracy of the MEAN pixel size to be at least as good as 0.1%
which is the limit to which I can measure our calibration plate and exceeds the
 parallax error in our calibration station.

Note, the Rayonix MX340 has the same pixel size as the MX170:

unbinned 340/7680  = 0.04427

  @param bin_size rayonix bin size as an integer
  '''
  pixel_size=bin_size*170/3840
  return pixel_size

def get_rayonix_detector_dimensions(env):
  ''' Given a psana env object, find the detector dimensions
      @param env psana environment object
  '''
  import psana
  cfgs = env.configStore()
  rayonix_cfg = cfgs.get(psana.Rayonix.ConfigV2, psana.Source('Rayonix'))
  if not rayonix_cfg: return None, None
  return rayonix_cfg.width(), rayonix_cfg.height()

def get_rayonix_cbf_handle(tiles, metro, timestamp, cbf_root, wavelength, distance, bin_size, detector_size, verbose = True, header_only = False):
  # set up the metrology dictionary to include axis names, pixel sizes, and so forth
  dserial = None
  dname = None
  detector_axes_names = [] # save these for later
  pixel_size = get_rayonix_pixel_size(bin_size)
  for key in sorted(metro):
    basis = metro[key]
    if len(key) == 1:
      assert dserial is None # only one detector allowed for now
      dserial = key[0]

      dname = "AXIS_D%d"%dserial #XXX check if DS1 is here
      for a in ["_X","_Y","_Z","_R"]: detector_axes_names.append(dname+a)
      basis.equipment_component = "detector_arm"
      basis.depends_on = dname+"_X"
      basis.include_translation = False # don't include the translation in the rotation axis offset below, instead it will be
                                        # included in the axis_settings table below
      basis.pixel_size = (pixel_size,pixel_size)
      basis.dimension = detector_size
      from xfel.cxi.cspad_ana.rayonix_tbx import rayonix_saturated_value
      basis.trusted_range = (rayonix_min_trusted_value, rayonix_saturated_value)
    else:
      assert False # shouldn't be reached as it would indicate a hierarchy for this detector
    basis.axis_name = detector_axes_names[-1]

  # the data block is the root cbf node
  from xfel.cftbx.detector.cspad_cbf_tbx import cbf_wrapper
  import os
  cbf=cbf_wrapper()
  cbf.new_datablock(os.path.splitext(os.path.basename(cbf_root))[0].encode())

  # Each category listed here is preceded by the imageCIF description taken from here:
  # http://www.iucr.org/__data/iucr/cifdic_html/2/cif_img.dic/index.html

  """Data items in the DIFFRN category record details about the
   diffraction data and their measurement."""
  cbf.add_category("diffrn",["id"])
  cbf.add_row(["Rayonix"])

  """Data items in the DIFFRN_SOURCE category record details of
  the source of radiation used in the diffraction experiment."""
  cbf.add_category("diffrn_source", ["diffrn_id","source","type"])
  cbf.add_row(["Rayonix","xfel","LCLS Endstation"])

  """Data items in the DIFFRN_DETECTOR category describe the
   detector used to measure the scattered radiation, including
   any analyser and post-sample collimation."""
  cbf.add_category("diffrn_detector", ["diffrn_id","id","type","details","number_of_axes"])
  cbf.add_row(["Rayonix","MFX_Rayonix","Rayonix",".",str(len(detector_axes_names))])

  """Data items in the DIFFRN_DETECTOR_AXIS category associate
     axes with detectors."""
  # Note, does not include the fast and the slow axes
  cbf.add_category("diffrn_detector_axis",["detector_id","axis_id"])
  for name in detector_axes_names:
    cbf.add_row(["MFX_Rayonix",name])

  """Data items in the DIFFRN_DETECTOR_ELEMENT category record
   the details about spatial layout and other characteristics
   of each element of a detector which may have multiple elements."""
  cbf.add_category("diffrn_detector_element",["id","detector_id"])

  cbf.add_row(["ELE_Rayonix", "MFX_Rayonix"])

  """Data items in the DIFFRN_DATA_FRAME category record
   the details about each frame of data."""
  cbf.add_category("diffrn_data_frame",["id","detector_element_id","array_id","binary_id"])

  cbf.add_row(["FRAME1","ELE_Rayonix","ARRAY_Rayonix","1"])

  if not header_only:
    add_frame_specific_cbf_tables(cbf, wavelength, timestamp, [metro[k].trusted_range for k in metro.keys()])

  """Data items in the AXIS category record the information required
     to describe the various goniometer, detector, source and other
     axes needed to specify a data collection.  The location of each
     axis is specified by two vectors: the axis itself, given as a unit
     vector, and an offset to the base of the unit vector.  These vectors
     are referenced to a right-handed laboratory coordinate system with
     its origin in the sample or specimen"""
  # More detail here: http://www.iucr.org/__data/iucr/cifdic_html/2/cif_img.dic/Caxis.html
  # Note we also use two new columns not in the latest imageCIF dictionary: rotation and rotation_axis.
  # We use them to specify an translation and a rotation in a single axis to describe a change in setting
  # when moving from one frame (say a quadrant) to another (say a sensor)
  cbf.add_category("axis",["id","type","equipment","depends_on","vector[1]","vector[2]","vector[3]",
                                                                "offset[1]","offset[2]","offset[3]",
                                                                "equipment_component"])

  # Keep a list of rows to add to the scan frame axis table
  axis_settings = []
  # keep a list of rows to add to the scan axis table
  axis_names = []

  # Create a series of axis describing frame shifts from each level of the detector to the next
  cbf.add_row( "AXIS_SOURCE  general     source   .        0  0  1 . . . .".split())                           ; axis_names.append("AXIS_SOURCE")
  cbf.add_row( "AXIS_GRAVITY general     gravity  .        0 -1  0 . . . .".split())                           ; axis_names.append("AXIS_GRAVITY")
  cbf.add_row(("%s_Z         translation detector .        0  0  1 . . . detector_arm"%(dname)).split())       ; axis_names.append("%s_Z"%dname)
  cbf.add_row(("%s_Y         translation detector %s_Z     0  1  0 . . . detector_arm"%(dname,dname)).split()) ; axis_names.append("%s_Y"%dname)
  cbf.add_row(("%s_X         translation detector %s_Y     1  0  0 . . . detector_arm"%(dname,dname)).split()) ; axis_names.append("%s_X"%dname)

  root_basis = metro[(0,)]

  axis_settings.append(["AXIS_SOURCE" ,"FRAME1","0","0"])
  axis_settings.append(["AXIS_GRAVITY","FRAME1","0","0"])
  axis_settings.append([dname+"_X"    ,"FRAME1","0",str(root_basis.translation[0])])
  axis_settings.append([dname+"_Y"    ,"FRAME1","0",str(root_basis.translation[1])])
  axis_settings.append([dname+"_Z"    ,"FRAME1","0",str(root_basis.translation[2])])

  for key in sorted(metro):
    basis = metro[key]
    cbf.add_frame_shift(basis, axis_settings)
    axis_names.append(basis.axis_name)

    dim_pixel = basis.pixel_size
    dim_readout = basis.dimension

    # Add the two vectors for each asic that describe the fast and slow directions pixels should be laid out in real space
    offset_fast = -dim_pixel[0]*((dim_readout[0]) / 2)
    offset_slow = +dim_pixel[1]*((dim_readout[1]) / 2)

    aname = "D%d"%key

    cbf.add_row(["AXIS_"+ aname + "_S", "translation","detector",basis.axis_name,"0", "-1","0",str(offset_fast),str(offset_slow),"0.0", "detector_asic"])
    cbf.add_row(["AXIS_"+ aname + "_F", "translation","detector","AXIS_"+aname +"_S","1","0","0","0","0","0.0", "detector_asic"])
    axis_names.append("AXIS_"+ aname + "_F"); axis_names.append("AXIS_"+ aname + "_S")
    axis_settings.append(["AXIS_"+ aname + "_F","FRAME1","0","0"])
    axis_settings.append(["AXIS_"+ aname + "_S","FRAME1","0","0"])

  """Data items in the DIFFRN_SCAN_AXIS category describe the settings of
     axes for particular scans.  Unspecified axes are assumed to be at
     their zero points."""
  # leave all the settings zero. Levels with settings are set below.
  cbf.add_category("diffrn_scan_axis",["axis_id","scan_id","angle_start","angle_range","angle_increment",
                                       "displacement_start","displacement_range","displacement_increment"])
  for name in axis_names:
    cbf.add_row([name,"SCAN1","0","0","0","0","0","0"])

  """Data items in the DIFFRN_SCAN_FRAME_AXIS category describe the
     settings of axes for particular frames.  Unspecified axes are
     assumed to be at their zero points."""
  cbf.add_category("diffrn_scan_frame_axis",["axis_id","frame_id","angle","displacement"])
  for row in axis_settings:
    cbf.add_row(row)

  """Data items in the ARRAY_STRUCTURE_LIST category record the size
     and organization of each array dimension.
     The relationship to physical axes may be given."""

  # find the asic sizes
  for tilekey in metro.keys():
    b = metro[tilekey]
    if not "x_dim" in locals():
      x_dim = b.dimension[0]
      y_dim = b.dimension[1]
    else:
      assert x_dim == b.dimension[0] and y_dim == b.dimension[1]

  z_dim = len(metro)
  root_key = (0,); assert root_key in metro; assert z_dim == 1 # single panel monolithic

  cbf.add_category("array_structure_list",["array_id","index","dimension","precedence","direction","axis_set_id"])
  cbf.add_row(["ARRAY_Rayonix","1","%d"%(x_dim),"1","increasing",dname+"_F"])
  cbf.add_row(["ARRAY_Rayonix","2","%d"%y_dim,"2","increasing",dname+"_S"])

  """Data items in the ARRAY_STRUCTURE_LIST_SECTION category identify
     the dimension-by-dimension start, end and stride of each section of an
     array that is to be referenced."""
  # no array sections in monolithic rayonix
  #cbf.add_category("array_structure_list_section",["array_id","id","index","start","end"])

  """Data items in the ARRAY_STRUCTURE_LIST_AXIS category describe
     the physical settings of sets of axes for the centres of pixels that
     correspond to data points described in the
     ARRAY_STRUCTURE_LIST category."""
  cbf.add_category("array_structure_list_axis",["axis_set_id","axis_id","displacement","displacement_increment"])
  cbf.add_row([dname+"_F",dname+"_F","0.0",str(metro[root_key].pixel_size[0])])
  cbf.add_row([dname+"_S",dname+"_S","0.0",str(metro[root_key].pixel_size[1])])

  if not header_only:
    add_data_to_cbf(cbf, tiles)

  return cbf

from dxtbx.format.FormatCBFCspad import FormatCBFFullStillInMemory
class FormatCBFRayonixInMemory(FormatCBFFullStillInMemory):
  """Mixin class for Rayonix in memory"""

def get_dxtbx_from_params(params, detector_size):
  """ Build a dxtbx format object for the Rayonix based on input paramters (beam center and binning) """
  from xfel.cftbx.detector.cspad_cbf_tbx import basis
  from scitbx.matrix import col
  fake_distance = 100
  null_ori = col((0,0,1)).axis_and_angle_as_unit_quaternion(0, deg=True)
  if params.override_beam_x is None and params.override_beam_y is None:
    metro = {(0,): basis(null_ori, col((0, 0, 0)))}
  elif params.override_beam_x is not None and params.override_beam_y is not None:
    # compute the offset from the origin given the provided beam center override
    pixel_size = get_rayonix_pixel_size(params.bin_size)
    image_center = col(detector_size)/2
    override_center = col((params.override_beam_x, params.override_beam_y))
    delta = (image_center-override_center)*pixel_size
    metro = {(0,): basis(null_ori, col((delta[0], -delta[1], 0)))} # note the -Y
  else:
    assert False, "Please provide both override_beam_x and override_beam_y or provide neither"
  cbf = get_rayonix_cbf_handle(None, metro, None, "test", None, fake_distance, params.bin_size, detector_size, verbose = True, header_only = True)
  base_dxtbx = FormatCBFRayonixInMemory(cbf)
  return base_dxtbx

def get_data_from_psana_event(evt, address):
  """ Read the pixel data for a Rayonix image from an event
  @param psana event object
  @param address old style psana detector address
  @return numpy array with raw data"""
  from psana import Source, Camera
  from xfel.cxi.cspad_ana import cspad_tbx
  import numpy as np
  address = cspad_tbx.old_address_to_new_address(address)
  src=Source('DetInfo(%s)'%address)
  data = evt.get(Camera.FrameV1,src)
  if data is not None:
    data = data.data16().astype(np.float64)
  return data

def format_object_from_data(base_dxtbx, data, distance, wavelength, timestamp, address, round_to_int=True):
  """
  Given a preloaded dxtbx format object and raw data, assemble the tiles
  and set the distance.
  @param base_dxtbx A header only dxtbx format object
  @param data rayonix byte array from XTC stream
  @param distance Detector distance (mm)
  @param wavelength Shot wavelength (angstroms)
  @param timestamp Human readable timestamp
  @param address Detector address, put in CBF header
  """
  import numpy as np
  from xfel.cftbx.detector import cspad_cbf_tbx
  cbf = cspad_cbf_tbx.copy_cbf_header(base_dxtbx._cbf_handle, skip_sections=True)
  rayonix_img = FormatCBFRayonixInMemory(cbf)
  cbf.set_datablockname((address + "_" + timestamp).encode())

  if round_to_int:
    data = flex.double(data.astype(np.float64)).iround()
  else:
    data = flex.double(data.astype(np.float64))

  n_asics = data.focus()[0] * data.focus()[1]
  cspad_cbf_tbx.add_frame_specific_cbf_tables(cbf, wavelength,timestamp,
    [(rayonix_min_trusted_value, rayonix_saturated_value)]*n_asics)

  # Set the distance, I.E., the length translated along the Z axis
  cbf.find_category(b"diffrn_scan_frame_axis")
  cbf.find_column(b"axis_id")
  cbf.find_row(b"AXIS_D0_Z") # XXX discover the Z axis somehow, don't use D0 here
  cbf.find_column(b"displacement")
  cbf.set_value(b"%f"%(-distance))

  # Explicitly reset the detector object now that the distance is set correctly
  rayonix_img._detector_instance = rayonix_img._detector()

  # Explicitly set up the beam object now that the tables are all loaded correctly
  rayonix_img._beam_instance = rayonix_img._beam()

  # Get the data and add it to the cbf handle.
  add_data_to_cbf(cbf,data)

  return rayonix_img

def add_data_to_cbf(cbf, data, verbose = False):
  """
  Given a cbf handle, add the raw data and the necessary tables to support it
  """
  import pycbf

  cbf.find_category(b"diffrn_data_frame")
  while True:
    try:
      cbf.find_column(b"array_id")
      array_name = cbf.get_value().decode()
      cbf.next_row()
    except Exception as e:
      assert "CBF_NOTFOUND" in str(e)
      break

  assert len(data.focus()) == 2
  if isinstance(data,flex.int):
    data_is_int = True
  elif isinstance(data,flex.double):
    data_is_int = False
  else:
    raise TypeError("Ints or doubles are required")

  """ Data items in the ARRAY_STRUCTURE category record the organization and
  encoding of array data in the ARRAY_DATA category."""
  cbf.add_category("array_structure",["id","encoding_type","compression_type","byte_order"])
  if data_is_int:
    cbf.add_row([array_name,"signed 32-bit integer","packed","little_endian"])
  else:
    cbf.add_row([array_name,"signed 64-bit real IEEE","packed","little_endian"])

  """ Data items in the ARRAY_DATA category are the containers for the array data
  items described in the category ARRAY_STRUCTURE. """
  cbf.add_category("array_data",["array_id","binary_id","data"])

  if verbose:
    print("Compressing tiles...", end=' ')

  focus = data.focus()

  cbf.add_row([array_name,"1"])

  binary_id = 1
  elements = len(data)
  data = data.copy_to_byte_str()
  byteorder = b"little_endian"
  dimfast = focus[1]
  dimmid = focus[0]
  dimslow = 1
  padding = 0

  if data_is_int:
    elsize = 4
    elsigned = 1

    cbf.set_integerarray_wdims_fs(\
      pycbf.CBF_PACKED,
      binary_id,
      data,
      elsize,
      elsigned,
      elements,
      byteorder,
      dimfast,
      dimmid,
      dimslow,
      padding)
  else:
    elsize = 8

    cbf.set_realarray_wdims_fs(\
      pycbf.CBF_CANONICAL,
      binary_id,
      data,
      elsize,
      elements,
      byteorder,
      dimfast,
      dimmid,
      dimslow,
      padding)
