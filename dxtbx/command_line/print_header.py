from __future__ import division
def print_header():
  import sys
  from dxtbx.format.Registry import Registry
  from scitbx.array_family import flex

  # this will do the lookup for every frame - this is strictly not needed
  # if all frames are from the same instrument

  for arg in sys.argv[1:]:
    print '=== %s ===' % arg
    format_class = Registry.find(arg)
    print 'Using header reader: %s' % format_class.__name__
    i = format_class(arg)
    print i.get_beam()
    print i.get_goniometer()
    print i.get_detector()
    print i.get_scan()
    from dxtbx.format.FormatMultiImage import FormatMultiImage
    if not issubclass(format_class, FormatMultiImage):
      try:
        d = i.get_raw_data()
        d = d.set_selected((d < 0), 0)
        print 'Total Counts: %d' % flex.sum(d)
      except AttributeError, e:
        print e

if __name__ == '__main__':
  print_header()
