from __future__ import absolute_import, division
from __future__ import print_function

from builtins import range
def saturation(image_file):
  from dxtbx import load
  i = load(image_file)
  d = i.get_detector()
  raw_data = i.get_raw_data()
  if not isinstance(raw_data, tuple):
    raw_data = (raw_data,)
  if i.get_scan() is None:
    return 0, \
         max([max(raw_data[pid]) / d[pid].get_trusted_range()[1] for pid in range(len(d))])
  else:
    return i.get_scan().get_image_range()[0], \
         max([max(raw_data[pid]) / d[pid].get_trusted_range()[1] for pid in range(len(d))])

if __name__ == '__main__':
  import sys

  for image_file in sys.argv[1:]:
    i, s = saturation(image_file)
    print('%6d %.6f' % (i, s))
