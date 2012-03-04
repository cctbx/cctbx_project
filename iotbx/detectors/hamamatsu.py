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

