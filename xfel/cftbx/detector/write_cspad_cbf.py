from __future__ import division
#-*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id: write_cspad_cbf.py
#

import pycbf, os

class cbf_wrapper(pycbf.cbf_handle_struct):
  """ Wrapper class that provids convience functions for working with cbflib"""

  def add_category(self, name, columns):
    """ Create a new category and populate it with column names """
    self.new_category(name)
    for column in columns:
      self.new_column(column)

  def add_row(self, data):
    """ Add a row to the current category.  If data contains more entries than
      there are columns in this category, then the remainder is truncated
      Use '.' for an empty value in a row. """
    self.new_row()
    self.rewind_column()
    for item in data:
      self.set_value(item)
      if item == '.':
        self.set_typeofvalue("null")
      try:
        self.next_column()
      except Exception:
        break

  def add_frame_shift(self, name, parent, params, axis_settings, equipment_component):
    """Add an axis representing a frame shift (a rotation axis with an offset)"""
    angle, axis = angle_and_axis(params)

    if angle == 0:
      axis = (0,0,1)

    self.add_row([name,"rotation","detector",parent,
                  str(axis[0]),str(axis[1]),str(axis[2]),
                  str(params.translation[0]*1000),
                  str(params.translation[1]*1000),
                  str(params.translation[2]*1000),
                  equipment_component])

    axis_settings.append([name, "FRAME1", str(angle), "0"])

def get_tile_metrology(detector, key):
  """ Searches a metrology phil object for an asic that matches the given key.
  @param  detector  metrology phil object
  @param  key       4-tuple in the form (detector, panel, sensor, key)"""
  kdetector,kpanel,ksensor,kasic = key

  assert kdetector == detector.serial

  for panel in detector.panel:
    if panel.serial == kpanel:
      for sensor in panel.sensor:
        if sensor.serial == ksensor:
          for asic in sensor.asic:
            if asic.serial == kasic:
              return asic

  return None

from scitbx import matrix

def angle_and_axis(params):
  """Normalize a quarternion and return the angle and axis
  @param params metrology object"""
  q = matrix.col(params.orientation).normalize()
  return q.unit_quaternion_as_axis_and_angle(deg=True)

def write_cspad_cbf(tiles, metro, timestamp, destpath):

  # the metrology object contains translations and rotations for each asic
  for o in metro.objects:
      if o.is_scope:
        #one of the subkeys of the root object will be the detector phil. it will be the only one not extracted.
        detector_phil = o.extract()
        break
  metro = metro.extract()

  cbf=cbf_wrapper()

  # the data block is the root cbf node
  cbf.new_datablock(os.path.splitext(os.path.basename(destpath))[0])

  # Each category listed here is preceded by the imageCIF description taken from here:
  # http://www.iucr.org/__data/iucr/cifdic_html/2/cif_img.dic/index.html

  """Data items in the DIFFRN category record details about the
   diffraction data and their measurement."""
  cbf.add_category("diffrn",["id"])
  cbf.add_row(["DS1"])

  """Data items in the DIFFRN_SOURCE category record details of
  the source of radiation used in the diffraction experiment."""
  cbf.add_category("diffrn_source", ["diffrn_id","source","type"])
  cbf.add_row(["DS1","xfel","LCLS CXI Endstation"])

  """Data items in the DIFFRN_RADIATION category describe
   the radiation used for measuring diffraction intensities,
   its collimation and monochromatization before the sample.

   Post-sample treatment of the beam is described by data
   items in the DIFFRN_DETECTOR category."""
  cbf.add_category("diffrn_radiation", ["diffrn_id","wavelength_id","probe"])
  cbf.add_row(["DS1","WAVELENGTH1","x-ray"])

  """ Data items in the DIFFRN_RADIATION_WAVELENGTH category describe
   the wavelength of the radiation used in measuring the diffraction
   intensities. Items may be looped to identify and assign weights
   to distinct wavelength components from a polychromatic beam."""
  cbf.add_category("diffrn_radiation_wavelength", ["id","wavelength","wt"])
  cbf.add_row(["WAVELENGTH1","%f"%metro.wavelength,"1.0"])

  """Data items in the DIFFRN_DETECTOR category describe the
   detector used to measure the scattered radiation, including
   any analyser and post-sample collimation."""
  cbf.add_category("diffrn_detector", ["diffrn_id","id","type","details","number_of_axes"])
  # figure out how many axes are in this detector
  dname = "AXIS_D%d"%detector_phil.serial #XXX check if DS1 is here
  detector_axes_names = [dname+a for a in ["_X","_Y","_Z","_R"]]

  for p in detector_phil.panel:
    detector_axes_names.append("FS_D%dQ%d"%(detector_phil.serial,p.serial))
    for s in p.sensor:
      detector_axes_names.append("FS_D%dQ%dS%d"%(detector_phil.serial,p.serial,s.serial))
      for a in s.asic:
        detector_axes_names.append("FS_D%dQ%dS%dA%d"%(detector_phil.serial,p.serial,s.serial,a.serial))

  cbf.add_row(["DS1","CSPAD_FRONT","CS PAD",".",str(len(detector_axes_names))])

  """Data items in the DIFFRN_DETECTOR_AXIS category associate
     axes with detectors."""
  # Note, does not include the fast and the slow axes
  cbf.add_category("diffrn_detector_axis",["detector_id","axis_id"])
  for name in detector_axes_names:
    cbf.add_row(["CSPAD_FRONT",name])

  tilestrs = []
  tilekeys = sorted(tiles)

  # create a series of strings representing each asic.  Q here is for quadrant instead of panel.
  for tilekey in tilekeys:
    tilestrs.append("D%dQ%dS%dA%d"%tilekey)

  """Data items in the DIFFRN_DETECTOR_ELEMENT category record
   the details about spatial layout and other characteristics
   of each element of a detector which may have multiple elements."""
  cbf.add_category("diffrn_detector_element",["id","detector_id"])

  for tilename in tilestrs:
    cbf.add_row(["ELE_" + tilename, "CSPAD_FRONT"])

  """Data items in the DIFFRN_DATA_FRAME category record
   the details about each frame of data."""
  cbf.add_category("diffrn_data_frame",["id","detector_element_id","array_id","binary_id"])

  for i, tilename in enumerate(tilestrs):
    cbf.add_row(["FRAME1","ELE_"+tilename,"ARRAY_"+tilename,"%d"%(i+1)])

  """Data items in the DIFFRN_MEASUREMENT category record details
   about the device used to orient and/or position the crystal
   during data measurement and the manner in which the
   diffraction data were measured."""
  cbf.add_category("diffrn_measurement",["diffrn_id","id","number_of_axes","method","details"])
  cbf.add_row(["DS1","INJECTION","0","electrospray","crystals injected by electrospray"])

  """ Data items in the DIFFRN_SCAN category describe the parameters of one
     or more scans, relating axis positions to frames."""
  cbf.add_category("diffrn_scan",["id","frame_id_start","frame_id_end","frames"])
  cbf.add_row(["SCAN1","FRAME1","FRAME1","1"])

  """Data items in the DIFFRN_SCAN_FRAME category describe
   the relationships of particular frames to scans."""
  cbf.add_category("diffrn_scan_frame",["frame_id","frame_number","integration_time","scan_id","date"])
  cbf.add_row(["FRAME1","1","0.0","SCAN1",timestamp])

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

  axis_settings.append(["AXIS_SOURCE" ,"FRAME1","0","0"])
  axis_settings.append(["AXIS_GRAVITY","FRAME1","0","0"])
  axis_settings.append([dname+"_X"    ,"FRAME1","0","0"])
  axis_settings.append([dname+"_Y"    ,"FRAME1","0","0"])
  axis_settings.append([dname+"_Z"    ,"FRAME1","0",str(-metro.distance)])

  cbf.add_frame_shift(dname+"_R",dname+"_X",detector_phil,axis_settings, "detector_arm")
  axis_names.append(dname+"_R")

  dname = "D%d"%detector_phil.serial

  for panel_phil in detector_phil.panel:
    pname = dname+"Q%d"%panel_phil.serial
    cbf.add_frame_shift("FS_"+pname,"AXIS_"+dname+"_R",panel_phil,axis_settings, "detector_quadrant")
    axis_names.append("FS_"+pname)
    for sensor_phil in panel_phil.sensor:
      sname = pname+"S%d"%sensor_phil.serial
      cbf.add_frame_shift("FS_"+sname,"FS_"+pname,sensor_phil,axis_settings, "detector_sensor")
      axis_names.append("FS_"+sname)
      for asic_phil in sensor_phil.asic:
        aname = sname+"A%d"%asic_phil.serial
        cbf.add_frame_shift("FS_"+aname,"FS_"+sname,asic_phil,axis_settings, "detector_asic")
        axis_names.append("FS_"+aname)

        dim_pixel = asic_phil.pixel_size
        dim_readout = asic_phil.dimension

        # Add the two vectors for each asic that describe the fast and slow dirctions pixels should be laid out in real space
        offset_fast = -dim_pixel[0]*1e3*((dim_readout[0] - 1) / 2)
        offset_slow = +dim_pixel[1]*1e3*((dim_readout[1] - 1) / 2)

        cbf.add_row(["AXIS_"+ aname + "_F", "translation","detector","FS_"  + aname      ,"1","0","0","%f"%offset_fast,"%f"%offset_slow,"0.0", "detector_asic"])
        cbf.add_row(["AXIS_"+ aname + "_S", "translation","detector","AXIS_"+ aname +"_F","0","-1","0","0","0","0.0", "detector_asic"])
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
  cbf.add_category("array_structure_list",["array_id","index","dimension","precedence","direction","axis_set_id"])
  for tilename,tilekey in zip(tilestrs,tilekeys):
    data = tiles[tilekey]
    cbf.add_row(["ARRAY_"+tilename,"1","%d"%data.focus()[1],"2","increasing","AXIS_"+tilename+"_F"])
    cbf.add_row(["ARRAY_"+tilename,"2","%d"%data.focus()[0],"1","increasing","AXIS_"+tilename+"_S"])

  """Data items in the ARRAY_STRUCTURE_LIST_AXIS category describe
     the physical settings of sets of axes for the centres of pixels that
     correspond to data points described in the
     ARRAY_STRUCTURE_LIST category."""
  cbf.add_category("array_structure_list_axis",["axis_set_id","axis_id","displacement","displacement_increment"])
  for tilename,tilekey in zip(tilestrs,tilekeys):
    tm = get_tile_metrology(detector_phil,tilekey)

    cbf.add_row(["AXIS_"+tilename+"_F","AXIS_"+tilename+"_F","0.0","%f"%(tm.pixel_size[0]*1000)])
    cbf.add_row(["AXIS_"+tilename+"_S","AXIS_"+tilename+"_S","0.0","%f"%(tm.pixel_size[1]*1000)])

  """ Data items in the ARRAY_INTENSITIES category record the
   information required to recover the intensity data from
   the set of data values stored in the ARRAY_DATA category."""
  # More detail here: http://www.iucr.org/__data/iucr/cifdic_html/2/cif_img.dic/Carray_intensities.html
  cbf.add_category("array_intensities",["array_id","binary_id","linearity","gain","gain_esd","overload","undefined_value"])
  for i, (tilename, tilekey) in enumerate(zip(tilestrs,tilekeys)):
    tm = get_tile_metrology(detector_phil,tilekey)
    cbf.add_row(["ARRAY_"+ tilename,str(i+1),"linear","1.0","0.1",str(tm.saturation),"0.0"])

  """ Data items in the ARRAY_STRUCTURE category record the organization and
     encoding of array data in the ARRAY_DATA category."""
  cbf.add_category("array_structure",["id","encoding_type","compression_type","byte_order"])

  for tilename in tilestrs:
    cbf.add_row(["ARRAY_"+tilename,"signed 64-bit real IEEE","canonical","little_endian"])

  """ Data items in the ARRAY_STRUCTURE category record the organization and
     encoding of array data in the ARRAY_DATA category."""
  cbf.add_category("array_data",["array_id","binary_id","data"])

  for i, tilekey in enumerate(tilekeys):
    detector, quadrant, sensor, asic = tilekey
    focus = tiles[tilekey].focus()

    print i, "Compressing tile D %s Q %s S %s A %s"%tilekey

    cbf.add_row(["ARRAY_" + tilestrs[i],str(i+1)])

    binary_id = i+1
    data = tiles[tilekey].copy_to_byte_str()
    elsize = 8
    elements = len(tiles[tilekey])
    byteorder = "little_endian"
    dimfast = focus[1]
    dimmid = focus[0]
    dimslow = 1
    padding = 0

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


  cbf.write_widefile(destpath,pycbf.CBF,\
      pycbf.MIME_HEADERS|pycbf.MSG_DIGEST|pycbf.PAD_4K,0)
