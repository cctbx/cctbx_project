from __future__ import absolute_import, division, print_function
import math,time
from iotbx.detectors.context.camera_convention import Cameras
from iotbx.detectors.context.config_detector import ADSC910_at_BioCARS

class EndStation:
  def __init__(self):
    #set defaults
    self.mos = {}
    self.mos['mosflm_detector']="""#detector-take defaults\n"""
    self.mos['mosflm_beamline']=beamlines['ALS']
    self.set_camera_convention(1)
    self.set_rotation_axis("ROTATION HORIZ ANTI")

  def set_camera_convention(self,number):
    self.cam_con = Cameras(number)

  def camera_convention(self,):
    return self.cam_con

  def set_rotation_axis(self,axis):
    if type(axis) == type("string"):
      self.rot_axi = rotation_lookup[axis]
      self.rot_axi_string = axis
    else: #tuple type of 3-component vector elements of normalized direction
      self.rot_axi= axis

  def rotation_axis(self):
    return self.rot_axi

  def mosflm(self): return self.mos

beamlines = {
"ALS":
"""#beam
SYNCHROTRON POLARIZATION 0.9
DIVERGENCE 0.100 0.020
DISPERSION 0.0001
""",
"Australian":
"""#beam
SYNCHROTRON POLAR 0.9
DIVERGENCE 0.11 0.001
DISPER 0.001
""",
"CHESS":
"""#beam
SYNCHROTRON POLAR 0.89
DIVE 0.030 0.010
DISPER 0.0025
""",
"RAXIS":
"""#beam
""",
}

rotation_lookup = {
"ROTATION HORIZ ANTI":(0,1,0),
"ROTATION HORIZ CLOCK":(0,-1,0), #reverse phi
"ROTATION VERT ANTI":(-1,0,0),
"ROTATION VERT CLOCK":(1,0,0),
}
#ALS Quantum 4:
#omega 90 rotation horiz anti fast horiz origin ur rectang type adsc
#rmin 5 rmax 188.0 ymax 188.0 xmax 188.0 xscan 188.0 yscan 188.0

def EndStation_from_ImageObject(imageobject,phil_params):
  endstation = EndStation()
  endstation.set_camera_convention(1)
  endstation.set_rotation_axis("ROTATION HORIZ ANTI")

  import six
  if isinstance(imageobject.parameters["DETECTOR_SN"],six.string_types) and \
    "S/N E-32-0105" in imageobject.parameters["DETECTOR_SN"]:
     # vertical goniometer axis at Max-IV
     endstation.set_rotation_axis("ROTATION VERT ANTI")
     print("MAX-IV Eiger 16M")
     endstation.mos['mosflm_detector'] = """
# Specific implementation for Max IV BioMAX
DETECTOR EIGE OMEGA 270
"""

  if imageobject.vendortype == "Bruker Proteus CCD":
     endstation.set_rotation_axis("ROTATION VERT ANTI")
     print("BRUKER rotation", endstation.rot_axi)

  if imageobject.vendortype == "RAXIS":
     endstation.set_rotation_axis("ROTATION VERT CLOCK")

  #clockwise horizontal phi at most CHESS beamlines
  #also at Australian Synchrotron
  if imageobject.vendortype == "ADSC" and \
     imageobject.serial_number in [406,414,441,448,457,471,924,928]:
     endstation.set_rotation_axis("ROTATION HORIZ CLOCK")

  if imageobject.vendortype == "ADSC" and imageobject.serial_number == 910:
    if ADSC910_at_BioCARS(imageobject):
      endstation.set_rotation_axis("ROTATION HORIZ CLOCK")
    else:
      endstation.set_rotation_axis("ROTATION VERT ANTI") # just a hypothesis
  #vertical phi at APS 19ID
  if imageobject.vendortype == "ADSC" and \
     imageobject.serial_number in [914]:
     endstation.set_rotation_axis("ROTATION HORIZ CLOCK")

  if imageobject.vendortype == "RigakuSaturn":
     endstation.set_rotation_axis("ROTATION VERT ANTI")

  #change in phi axis rotation at CHESS A1
  if imageobject.vendortype == "ADSC" and \
     imageobject.serial_number in [441]:
       record_date = imageobject.parameters["DATE"]
       record_tse = time.mktime(time.strptime(record_date))
       cutoff_441 = time.mktime(time.strptime("Fri Oct 01 00:00:00 2004"))
       if record_tse < cutoff_441:
         endstation.set_rotation_axis("ROTATION HORIZ ANTI")

  # Special cases:  Pringle-Shen goniometer at CHESS F3
  if imageobject.vendortype == "ADSC" and \
     imageobject.serial_number in [414] and \
     phil_params.goniometer_rotation.lower().find(
       'pringle-shen')>=0 and \
     'AXIS' in imageobject.parameters and \
     imageobject.parameters['AXIS']=='phi':
       endstation.set_rotation_axis("ROTATION VERT ANTI")
       if 'OMEGA' in imageobject.parameters and \
          imageobject.parameters['OMEGA'] != 0.0 :
          omega = imageobject.parameters['OMEGA']
          from iotbx.detectors import rotate_vector_around
          endstation.set_rotation_axis(
             rotate_vector_around(endstation.rotation_axis(),
                        (0,1,0),-omega*math.pi/180.)
          )
  '''
  #tested for CHESS F1 s/n 406:
  SCANNER ROTATION HORIZ CLOCK FAST horizontal ORIGIN UR RECT TYPE ADSC
  LIMITS RMIN 5 RMAX 137.2 XMAX 96.5 YMAX 97.5 XSCAN 94.0 YSCAN 94.0
  BACKSTOP RADIUS 7.00 CENTRE 90.100 91.500
  ADCOFFSET 20
  NULLPIX 0
  GAIN 0.300
  #tested for CHESS F3 Pringle-Shen
  SCANNER ROTATION VERT ANTI FAST horizontal ORIGIN UR RECT TYPE ADSC
  LIMITS RMIN 5 RMAX 143.3 XMAX 106.8 YMAX 95.7 XSCAN 94.0 YSCAN 94.0
  BACKSTOP RADIUS 4.00 CENTRE 106.600 94.200
  GAIN 0.500
  BIAS 5
  '''
  if imageobject.vendortype == "MacScience":
    if imageobject.size1==3000:
      endstation.mos['mosflm_detector'] = """DETECTOR DIP2030\n"""
      endstation.set_rotation_axis("ROTATION HORIZ CLOCK")

  if imageobject.vendortype == "MARCCD":
   if imageobject.size1*imageobject.bin>=4096:
     parameters = {'sz' : imageobject.size1*imageobject.bin,
                   'pix': imageobject.pixel_size / imageobject.bin}
     parameters['scan']=parameters['pix']*parameters['sz']/2

     endstation.mos['mosflm_detector'] = """#MARCCD detector
LIMITS XMIN 0 XMAX xmax_tag YMIN 0 YMAX ymax_tag xscan %(scan)d yscan %(scan)d
SIZE %(sz)d %(sz)d HEADER 1 byte 4096
PIXEL %(pix)f
NULLPIX 0
#4Kx4K MarCCD format is unknown to MOSFLM, which apparently defaults
#the nullpix to 10.  This is a problem for weak-background images.
"""%parameters
     '''the correct xmax_tag and ymax_tag are added later in the interface module'''

   #identification of specific beamlines at Spring8 with Reversephi:
     # BL41XU Mar MX225HE--Serial number 40
     # BL32XU Mar MX225HE--Serial number 31
   # rely on detector serial number, uncertain how to decode geometric description
   # of rotation axis within the header.

   if imageobject.parameters["DETECTOR_SN"] in [7]:
     endstation.set_rotation_axis("ROTATION HORIZ CLOCK")
     endstation.mos['mosflm_detector'] = """
# Specific implementation for APS SER-CAT BM22, chi=180 setting
DETECTOR MARCCD
DETECTOR REVERSEPHI
SIZE 4096 4096
PIXEL 0.07324 0.07324
"""
     endstation.mos['mosflm_beamline'] = """GAIN 0.37
POLARISATION 0.99
DIVE 0.0001 0.00001
DISPER 0.0001
"""

   if imageobject.parameters["DETECTOR_SN"] in [31,40]:
     endstation.set_rotation_axis("ROTATION HORIZ CLOCK")
     endstation.mos['mosflm_detector'] = """
# Specific implementation for Spring8 BL41XU Mar MX225HE
DETECTOR MARCCD
DETECTOR REVERSEPHI
SIZE 3072 3072
PIXEL 0.07324 0.07324
LIMITS RMIN 2 RMAX 159.1 XMAX 159.1 YMAX 112.5 XSCAN 159.1 YSCAN 159.1
"""
     endstation.mos['mosflm_beamline'] = """GAIN 0.37
POLARISATION 0.99
DIVE 0.0001 0.00001
DISPER 0.0001
!offset 0.0 0.0
"""

  if imageobject.vendortype == "ADSC" and \
     endstation.rot_axi_string!="ROTATION HORIZ ANTI":

     # Rough idea of the MOSFLM LIMITS
     mmsize = imageobject.size1 * imageobject.pixel_size
     maxx = max( abs(imageobject.beamx-mmsize), abs(imageobject.beamx) )
     maxy = max( abs(imageobject.beamy-mmsize), abs(imageobject.beamy) )
     parameters = {'rotation': endstation.rot_axi_string,
                   'rmax': math.sqrt(maxx*maxx + maxy*maxy),
                   'xmax': maxx,
                   'ymax': maxy,
                   'xscan': mmsize/2.,
                   'yscan': mmsize/2.}
     if imageobject.serial_number in [457,928]: parameters['gain']=0.32
     else: parameters['gain']=0.30

     endstation.mos['mosflm_detector'] = """#detector
SCANNER %(rotation)s FAST horizontal ORIGIN UR RECT TYPE ADSC
LIMITS RMIN 5 RMAX %(rmax).1f XMAX %(xmax).1f YMAX %(ymax).1f XSCAN %(xscan).1f YSCAN %(yscan).1f
GAIN %(gain).3f
"""%parameters

     if imageobject.serial_number in [457,928]:
       endstation.mos['mosflm_beamline'] = beamlines['Australian']
     else:
       endstation.mos['mosflm_beamline'] = beamlines['CHESS']

  if imageobject.vendortype == "ADSC" and imageobject.serial_number == 910:
    if ADSC910_at_BioCARS(imageobject):
      endstation.mos['mosflm_detector']=endstation.mos['mosflm_detector']+\
       "#BIOCARS 14-BM-C S/N=910"

  if imageobject.vendortype in [ "RAXIS" ]:
     endstation.mos['mosflm_detector'] = """#detector
ADCOFFSet 5
"""
     if imageobject.serial_number.lower().find('dr. r-axis iv')==0:
       endstation.mos['mosflm_detector']=endstation.mos['mosflm_detector']+\
                                         'DETECTOR RAXISIV'
       #At least for MOSFLM 6.2.4, this seems to be important because it
       #  allows the program to accept large standard deviations for spots
       #  (above 32767).  Otherwise MOSFLM can crash with SERIOUS ERROR message.

     # Both DTrek and Raxis formats have "RAXIS" vendortype but only Raxis has "head" attribute
     if "head" in imageobject.__dict__ and \
        imageobject.head['Device'].lower().find('r-axis2')==0:
       endstation.mos['mosflm_detector']=endstation.mos['mosflm_detector']+\
                                         'detector raxis'
     endstation.mos['mosflm_beamline'] = beamlines['RAXIS']

  if imageobject.vendortype in [ "CBF" ]:
     endstation.mos['mosflm_detector'] = """#detector cbf"""

  if imageobject.vendortype in [ "Pilatus-6M" ]:
     endstation.mos['mosflm_detector'] = """#detector Pilatus-6M"""
     #flags mosflm interface to include the start & angle fix, mosflm 7.0.3 & below

  if imageobject.vendortype in [ "MARIP" ]:
     endstation.mos['mosflm_detector'] = """#detector
ADCOFFSet 5
"""#Just a guess: image plates require an offset to insure stability
   # over long integration sweeps.  Example: TM0064/11_20_01/1c3p3

  #additional information for universally specifying image size
  endstation.mos['mosflm_detector'] = endstation.mos['mosflm_detector'] + \
"\n#UIS_PIXEL %(pix)f\n#UIS_SIZE %(sz)d"%{
'sz' : imageobject.size1*imageobject.bin,
'pix': imageobject.pixel_size / imageobject.bin}
  return endstation
