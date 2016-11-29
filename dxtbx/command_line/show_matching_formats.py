from __future__ import division
def show_matching_formats():
  import sys
  from dxtbx.format.Registry import Registry

  for arg in sys.argv[1:]:
    print '=== %s ===' % arg
    for fmt in Registry.get():
      print fmt.__name__, fmt.understand(arg)

if __name__ == '__main__':
  show_matching_formats()
