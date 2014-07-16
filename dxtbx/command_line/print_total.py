from __future__ import division
def print_total():
  import sys
  from dxtbx.format.Registry import Registry

  # this will do the lookup for every frame - this is strictly not needed
  # if all frames are from the same instrument

  for arg in sys.argv[1:]:
    print '=== %s ===' % arg
    format_class = Registry.find(arg)
    print 'Using header reader: %s' % format_class.__name__
    i = format_class(arg)
    image_size = i.get_detector()[0].get_image_size()
    data = i.get_raw_data()
    selection = (data < 0)
    data = data.set_selected(selection, 0)
    total = sum(data)
    print 'Total Counts: %d' % total
    print 'Average Counts: %.2f' % (total / (image_size[0] * image_size[1]))

if __name__ == '__main__':
  print_total()
