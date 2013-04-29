from __future__ import division
from dxtbx.format.Format import Format
from dxtbx.format.FormatMultiImage import FormatMultiImage

class FormatHDF5(Format, FormatMultiImage):

    def __init__(self, image_file):
        assert(self.understand(image_file))
        Format.__init__(self, image_file)

    @staticmethod
    def understand(image_file):
        try:
            tag = FormatHDF5.open_file(image_file, 'rb').read(8)
        except IOError,e:
            return False

        return tag == "\211HDF\r\n\032\n"

if __name__ == '__main__':
    import sys
    for arg in sys.argv[1:]:
        print FormatHDF5.understand(arg)
