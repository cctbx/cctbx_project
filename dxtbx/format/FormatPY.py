# Implementation of a base class to read a pickled Python dictionary.

from __future__ import absolute_import, division
from __future__ import print_function

from builtins import chr
from dxtbx.format.Format import Format

class FormatPY(Format):
  '''Let's take an educated guess as to how to recognize a Python
  pickle file containing a dictionary.  Not easy because there are
  three pickle protocols in Python 2.7.  Dangerous because the lowest
  pickle format only gives us two unique bytes by which to recognize
  a dictionary.  Could possibly conflict with other image formats.'''

  @staticmethod
  def understand(image_file):
    try:
      tag = FormatPY.open_file(image_file, 'rb').read(4)
    except IOError as e:
      return False

    return tag[0:2]=="(d" or tag[0:2]=="}q" or \
           tag[0:4]==chr(0o200)+chr(0o002)+"}q"

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file.'''

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)
    Format.__init__(self, image_file, **kwargs)
    return

if __name__ == '__main__':
  import sys
  for arg in sys.argv[1:]:
    print(FormatPY.understand(arg))
