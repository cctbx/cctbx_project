from __future__ import division
from dxtbx.format.Format import Format
from dxtbx.format.FormatMultiImage import FormatMultiImage
from dxtbx.model import Beam # import dependency
from dxtbx.model import Detector # import dependency
from dxtbx.model import Goniometer # import dependency
from dxtbx.model import Scan # import dependency

class FormatHDF5(Format, FormatMultiImage):

  def __init__(self, image_file):
    assert(self.understand(image_file))
    FormatMultiImage.__init__(self)
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
