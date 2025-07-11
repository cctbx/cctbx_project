"""Edit a reflection file"""
from __future__ import absolute_import, division, print_function

import sys
from iotbx import reflection_file_editor

if __name__ == "__main__" :
  reflection_file_editor.run(sys.argv[1:])
