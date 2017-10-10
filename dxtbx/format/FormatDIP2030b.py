from __future__ import absolute_import, division
from __future__ import print_function

from dxtbx.format.Format import Format

class FormatDIP2030b(Format):

  @staticmethod
  def understand(image_file):
    # for MacScience DIP2030b only, file size is exactly 18001024 bytes
    headerstart = 3000 * 3000 * 2
    try:
      F = FormatDIP2030b.open_file(image_file, 'rb')
      F.seek(headerstart)
      rawheader = F.read(1024)
      eof = F.read(1)# end of file
      F.close()
    except IOError as e:
      return False

    return eof == "" and rawheader[0:3] == "DIP"

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file.'''

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    Format.__init__(self, image_file, **kwargs)

  def detectorbase_start(self): pass
  def _start(self):
    from iotbx.detectors.macscience import DIPImage
    self.detectorbase = DIPImage(self._image_file)
    self.detectorbase.readHeader()

  def _goniometer(self):

    return self._goniometer_factory.single_axis()

  def _detector(self):
    '''Return a model for a simple detector'''

    twotheta = self.detectorbase.parameters["TWOTHETA"]
    # At present, ignore non-zero two theta for the dxtbx model
    # XXX Return to this issue later.

    return self._detector_factory.simple(
        sensor = 'IMAGE_PLATE',
        distance = self.detectorbase.parameters["DISTANCE"],
        beam_centre = (self.detectorbase.parameters["BEAM_CENTER_X"],
                       self.detectorbase.parameters["BEAM_CENTER_Y"]),
        fast_direction = '+x',
        slow_direction = '-y',
        pixel_size = (self.detectorbase.parameters["PIXEL_SIZE"],
                      self.detectorbase.parameters["PIXEL_SIZE"]),
        image_size = (self.detectorbase.parameters["SIZE1"],
                      self.detectorbase.parameters["SIZE2"]),
        trusted_range = (0, self.detectorbase.parameters["CCD_IMAGE_SATURATION"]),
        mask = [])  # a list of dead rectangles

  def _beam(self):
    '''Return a simple model for the beam.'''

    return self._beam_factory.simple(self.detectorbase.parameters["WAVELENGTH"])

  def _scan(self):
    '''Return the scan information for this image.'''

    return self._scan_factory.single(
      filename = self._image_file,
      format = "DIP",
      exposure_times = self.detectorbase.parameters["TIME"],
      osc_start = self.detectorbase.parameters["OSC_START"],
      osc_width = self.detectorbase.parameters["OSC_RANGE"],
      epoch = None)

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print(FormatDIP2030b.understand(arg))
