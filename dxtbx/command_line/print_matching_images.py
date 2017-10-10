from __future__ import absolute_import, division
from __future__ import print_function

def print_matching_images(image):
  from dxtbx.sweep_filenames import find_matching_images
  matching_images = find_matching_images(image)
  for mi in matching_images:
    print(mi)

if __name__ == '__main__':

  import sys
  print_matching_images(sys.argv[1])
