from __future__ import absolute_import, division
from __future__ import print_function

from dxtbx.format.FormatSMVADSC import FormatSMVADSC

class FormatSMVHamamatsu(FormatSMVADSC):

  @staticmethod
  def understand(image_file):

    size, header = FormatSMVHamamatsu.get_smv_header(image_file)

    wanted_header_items = ['DETECTOR_NAME']

    for header_item in wanted_header_items:
      if not header_item in header:
        return 0

    return header["DETECTOR_NAME"].lower().find("hamamatsu")>=0

  def _start(self):

    FormatSMVADSC._start(self)

  def detectorbase_start(self):
    from iotbx.detectors.hamamatsu import HamamatsuImage
    self.detectorbase = HamamatsuImage(self._image_file)
    self.detectorbase.open_file = self.open_file
    self.detectorbase.readHeader()

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print(FormatSMVADSC.understand(arg))
