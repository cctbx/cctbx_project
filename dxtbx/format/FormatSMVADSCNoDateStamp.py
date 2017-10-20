from __future__ import absolute_import, division

from dxtbx.format.FormatSMVADSC import FormatSMVADSC

class FormatSMVADSCNoDateStamp(FormatSMVADSC):
  '''A class for reading SMV format ADSC images, with detector serial number'''

  @staticmethod
  def understand(image_file):

    # assert for this that the image file has to be a file not a URL

    import os
    if not os.path.exists(image_file):
      return False

    size, header = FormatSMVADSC.get_smv_header(image_file)

    wanted_header_items = ['TIME']

    for header_item in wanted_header_items:
      if not header_item in header:
        return False

    unwanted_header_items = ['DATE']

    for header_item in unwanted_header_items:
      if header_item in header:
        return False

    return True

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    FormatSMVADSC.__init__(self, image_file, **kwargs)

  def _start(self):

    FormatSMVADSC._start(self)

  def _scan(self):
    '''Return the scan information for this image, using the timestamp
    from the file rather than the header.'''

    format = self._scan_factory.format('SMV')
    exposure_time = float(self._header_dictionary['TIME'])

    import os
    epoch = float(os.stat(self._image_file)[8])

    osc_start = float(self._header_dictionary['OSC_START'])
    osc_range = float(self._header_dictionary['OSC_RANGE'])

    return self._scan_factory.single(
        self._image_file, format, exposure_time,
        osc_start, osc_range, epoch)

  def detectorbase_start(self):
    from iotbx.detectors.adsc import ADSCImage
    self.detectorbase = ADSCImage(self._image_file)
    self.detectorbase.open_file = self.open_file
    self.detectorbase.readHeader()

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatSMVADSCNoDateStamp.understand(arg)
