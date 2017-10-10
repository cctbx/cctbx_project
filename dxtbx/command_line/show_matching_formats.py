from __future__ import absolute_import, division
from __future__ import print_function

def recurse(fmt, arg):
  for _fmt in fmt._children:
    if _fmt.understand(arg):
      print(_fmt.__name__, _fmt.understand(arg))
      recurse(_fmt, arg)
    else:
      print(fmt.__name__, fmt.understand(arg))




def show_matching_formats():
  import sys
  from dxtbx.format.Registry import Registry

  for arg in sys.argv[1:]:
    print('=== %s ===' % arg)
    for fmt in Registry.get():
      if fmt.understand(arg):
        print(fmt.__name__, fmt.understand(arg))
        recurse(fmt, arg)
      else:
        print(fmt.__name__, fmt.understand(arg))

if __name__ == '__main__':
  show_matching_formats()
