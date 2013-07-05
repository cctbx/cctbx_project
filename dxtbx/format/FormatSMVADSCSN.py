from __future__ import division

from dxtbx.format.FormatSMVADSC import FormatSMVADSC

class FormatSMVADSCSN(FormatSMVADSC):
    '''A class for reading SMV format ADSC images, with detector serial number'''

    @staticmethod
    def understand(image_file):

        size, header = FormatSMVADSC.get_smv_header(image_file)

        wanted_header_items = ['DETECTOR_SN', 'DATE']

        for header_item in wanted_header_items:
            if not header_item in header:
                return 0

        return True

    def __init__(self, image_file):
        '''Initialise the image structure from the given file, including a
        proper model of the experiment.'''

        assert(self.understand(image_file))

        FormatSMVADSC.__init__(self, image_file)

        return

    def _start(self):

        from iotbx.detectors.adsc import ADSCImage
        self.detectorbase = ADSCImage(self._image_file)
        self.detectorbase.readHeader()

        FormatSMVADSC._start(self)


if __name__ == '__main__':

    import sys

    for arg in sys.argv[1:]:
        print FormatSMVADSCSN.understand(arg)
