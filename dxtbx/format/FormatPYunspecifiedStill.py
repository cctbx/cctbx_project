from __future__ import division

from dxtbx.format.FormatStill import FormatStill
from dxtbx.format.FormatPYunspecified import FormatPYunspecified
from dxtbx.format.FormatPYunspecified import FormatPYunspecifiedInMemory

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

    if 'OSC_START' not in data or 'OSC_RANGE' not in data:
      return True

    return data['OSC_RANGE'] <= 0


  def __init__(self, image_file):
    '''Initialise the image structure from the given file.'''

    assert(self.understand(image_file))

    FormatPYunspecified.__init__(self, image_file)

class FormatPYunspecifiedStillInMemory(FormatStill, FormatPYunspecifiedInMemory):
  """ Overrides the Format object's init method to accept an image dictionary
      instead of a file name. Used with XFELs when it is desirable to never write
      a file to disk, but to process it only in memory.
  """

  @staticmethod
  def understand(image_file):
    """ If it's an image dictionary, we understand this """
    data = image_file

    try:
      if 'OSC_START' not in data.keys() or 'OSC_RANGE' not in data.keys():
        return True
    except AttributeError:
      return False

    return data['OSC_RANGE'] <= 0

  def __init__(self, data):
    """ @param data In memory image dictionary, alredy initialized """
    FormatPYunspecifiedInMemory.__init__(self, data)

    import copy
    self._image_file = copy.deepcopy(data)

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatPYunspecifiedStill.understand(arg)
