from __future__ import absolute_import, division

from dxtbx.format.FormatSMVADSC import FormatSMVADSC

class FormatSMVADSCSN(FormatSMVADSC):
  '''A class for reading SMV format ADSC images, with detector serial number'''

  @staticmethod
  def understand(image_file):
    # The header must include a DATE and an integer-valued DETECTOR_SN
    # for this format to apply.

    size, header = FormatSMVADSC.get_smv_header(image_file)

    if 'DATE' not in header.keys():
      return False
    try:
      int(header['DETECTOR_SN'])
    except (KeyError, ValueError):
      return False

    return True

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    FormatSMVADSC.__init__(self, image_file, **kwargs)

    return

  def _start(self):

    FormatSMVADSC._start(self)

  def detectorbase_start(self):

    from iotbx.detectors.adsc import ADSCImage
    self.detectorbase = ADSCImage(self._image_file)
    self.detectorbase.open_file = self.open_file
    self.detectorbase.readHeader()

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatSMVADSCSN.understand(arg)
