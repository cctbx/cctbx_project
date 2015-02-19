from __future__ import division

from dxtbx.format.FormatStill import FormatStill
from dxtbx.format.FormatPYunspecified import FormatPYunspecified

class FormatPYunspecifiedStill(FormatStill, FormatPYunspecified):

  @staticmethod
  def understand(image_file):
    """Seems like the static method wastes a lot of effort here; it's not possible to
    just read the first few bytes; instead understand() reads the entire first data
    item in the file; an entire binary image.  This data is then read again in the
    _start() method and again in the detectorbase constructor.
    """

    try:
      stream = FormatPYunspecified.open_file(image_file, 'rb')
      import cPickle as pickle
      data = pickle.load(stream)
    except IOError,e:
      return False

    wanted_header_items = ['OSC_START','OSC_RANGE']

    for header_item in wanted_header_items:
      if not header_item in data:
        return True

    return data['OSC_RANGE'] <= 0


  def __init__(self, image_file):
    '''Initialise the image structure from the given file.'''

    assert(self.understand(image_file))

    FormatPYunspecified.__init__(self, image_file)


if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatPYunspecifiedStill.understand(arg)
