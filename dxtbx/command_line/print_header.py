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
    beam = i.get_beam()
    goniometer = i.get_goniometer()
    detector = i.get_detector()
    scan = i.get_scan()
    if beam is None:
      print 'No beam model found'
    else:
      print beam
    if detector is None:
      print 'No detector model found'
    else:
      print detector
    if goniometer is None:
      print 'No goniometer model found'
    else:
      print goniometer
    if scan is None:
      print 'No scan model found'
    else:
      print scan
    from dxtbx.format.FormatMultiImage import FormatMultiImage
    if not issubclass(format_class, FormatMultiImage):
      try:
        d = i.get_raw_data()
        d.set_selected((d < 0), 0)
        print 'Total Counts: %d' % flex.sum(d)
      except AttributeError, e:
        print "Could not read image data"

if __name__ == '__main__':
  print_header()
