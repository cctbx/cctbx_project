from __future__ import absolute_import, division, print_function
from iotbx.detectors.adsc import ADSCImage
from iotbx.detectors.detectorbase import DetectorImageBase

class HamamatsuImage(ADSCImage):
  """Hamamatsu CMOS Detector
  Beamline BL32XU, SPring-8
  RIKEN/SPring-8 Center
  Research Infrastructure Group,
  SR Life Science Instrumentation Unit

  1-1-1 Kouto Sayo-cho Sayo-gun
  Hyogo, 679-5148 JAPAN
  Contact: Kunio Hirata
  """
  def __init__(self,filename):
    DetectorImageBase.__init__(self,filename)
    self.vendortype = "Hamamatsu"

  def readHeader(self):
    ADSCImage.readHeader(self,external_keys=[("DETECTOR_NAME","DETECTOR_NAME",str)])
    assert self.parameters["DETECTOR_NAME"].lower().find("hamamatsu")>=0

    #above code validates the Hamamatsu signature, as in
    """HEADER_BYTES=512;
DIM=2;
BYTE_ORDER=little_endian;
TYPE=unsigned_short;
SIZE1=2352;
SIZE2=2352;
PIXEL_SIZE=0.050000;
BIN=1x1;
DETECTOR_NAME=Hamamatsu C10158DK;
DATE=Tue Jun 14 15:31:08 2011;
TIME=0.33;
DISTANCE=144.000;
OSC_RANGE=1.000;
OMEGA=1.000;
OSC_START=0.000;
TWOTHETA=0.000;
AXIS=Omega;
WAVELENGTH=1.00000;
BEAM_CENTER_X=58.425;
BEAM_CENTER_Y=58.775;
"""
