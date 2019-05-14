from __future__ import absolute_import, division, print_function
from libtbx import easy_pickle
import logging

class SingleFrame:
  """ Class that creates single-image agregate metrics/scoring that can then be
  used in downstream clustering or filtering procedures.
  """
  def __init__(self, path, filename, crystal_num=0):
    try:
      # Warn on error, but continue directory traversal.
      d = easy_pickle.load(path)
      self.miller_array = d['observations'][crystal_num]
      self.path = path
      self.name = filename
      self.pg = d['pointgroup']
      self.uc = d['current_orientation'][crystal_num].unit_cell() \
                                                     .niggli_cell() \
                                                     .parameters()
      self.total_i = d['observations'][crystal_num].sum()
      logging.debug("Extracted image {}".format(filename))
    except KeyError:
      logging.warning("Could not extract point group and unit cell from %s\n" % path)
    except IOError:
      logging.warning("Could not read %s. It may not be a pickle file.\n" % path)
