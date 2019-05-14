from __future__ import absolute_import, division, print_function
from six.moves import range
import os

def match_dials_files(list1, list2, suffix1, suffix2):
  ids1 = [os.path.join(os.path.dirname(item), os.path.basename(item).split(suffix1)[0])
          for item in list1]
  ids2 = [os.path.join(os.path.dirname(item), os.path.basename(item).split(suffix2)[0])
          for item in list2]
  matches = [(ids1[i] + suffix1,ids2[ids2.index(ids1[i])] + suffix2)
             for i in range(len(ids1)) if ids1[i] in ids2]
  return (zip(*matches))
