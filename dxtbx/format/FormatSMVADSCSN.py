from __future__ import absolute_import, division, print_function

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

    # Mapping of serial numbers to models for known detectors
    self._sn_to_model = {401: 'Q4U',
                         402: 'Q4',
                         414: 'Q4',
                         423: 'Q4R',
                         428: 'Q4R',
                         429: 'Q4', # or Q4R?
                         441: 'Q210',
                         442: 'Q210',
                         443: 'Q210',
                         444: 'Q210', # or Q210R?
                         445: 'Q210',
                         446: 'Q210',
                         447: 'Q210',
                         448: 'Q210',
                         457: 'Q210R',
                         471: 'Q270',
                         472: 'Q270',
                         474: 'Q270',
                         901: 'Q210',
                         905: 'Q315',
                         907: 'Q315R', # or Q315?
                         913: 'Q315',
                         917: 'Q315R',
                         923: 'Q315R',
                         925: 'Q315',
                         926: 'Q315R',
                         928: 'Q315R',
                         931: 'Q315R',
                         933: 'Q315R',
                         }

    FormatSMVADSC.__init__(self, image_file, **kwargs)

  def _adsc_module_gain(self, model=None):
    '''Overload to look the model number up from the serial number table'''

    if model is None:
      sn = int(self._header_dictionary['DETECTOR_SN'])
      model = self._sn_to_model.get(sn)
    return super(FormatSMVADSCSN, self)._adsc_module_gain(model=model)

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
    print(FormatSMVADSCSN.understand(arg))
