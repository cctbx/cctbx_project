""" Transfer of scalepack reflection files to flex arrays.
"""

# Sample scalepack OUTPUT FILE
#    1
# -987
#    34.698    71.491    54.740    90.000   106.549    90.000 p21
#   0   0   4  3617.6   287.2
#   0   1   6 12951.7  1583.6 12039.2  1665.8
#
# Format: (3I4, 4F8.1)

import exceptions
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import crystal
from cctbx import miller
from cctbx.array_family import flex

class ScalepackFormatError(exceptions.Exception): pass

class scalepack_reader:

  def __init__(self, file_handle):
    self.file_handle = file_handle
    line = self.file_handle.readline()
    if (line.rstrip() != "    1"):
      raise ScalepackFormatError, "line 1: expecting '    1'"
    line = self.file_handle.readline()
    if (line.rstrip() != " -987"):
      raise ScalepackFormatError, "line 2: expecting ' -987'"
    line_error = "line 3: expecting unit cell parameters and space group label"
    line = self.file_handle.readline()
    if (len(line) < 63 or line[60] != ' '):
      raise ScalepackFormatError, line_error
    try:
      uc_params = [float(line[i * 10 : (i + 1) * 10]) for i in xrange(6)]
    except:
      raise ScalepackFormatError, line_error
    self.unit_cell = uctbx.unit_cell(uc_params)
    self.space_group_symbol = line[61:].strip()
    if (len(self.space_group_symbol) == 0):
      raise ScalepackFormatError, line_error
    try:
      self.space_group_info = sgtbx.space_group_info(self.space_group_symbol)
    except:
      self.space_group_info = None
    self.miller_indices = flex.miller_index()
    self.fobs = flex.double()
    self.sigmas = flex.double()
    self.anomalous = 0
    line_count = 3
    while 1:
      line = self.file_handle.readline()
      line_count += 1
      line_error = "line %d: expecting (3I4, 4F8.1)" % line_count
      if (line == ""): break
      line = line.rstrip() + (" " * 44)
      flds = []
      used = 0
      for width in (4,4,4,8,8,8,8):
        next_used = used + width
        flds.append(line[used:next_used].strip())
        used = next_used
      try:
        h = [int(flds[i]) for i in xrange(3)]
      except:
        raise ScalepackFormatError, line_error
      for i in (0,1):
        j = 3+2*i
        if (len(flds[j])):
          try:
            fobs, sigma = (float(flds[j]), float(flds[j+1]))
          except:
            raise ScalepackFormatError, line_error
          if (i):
            h = [-e for e in h]
            self.anomalous = 1
          self.miller_indices.append(h)
          self.fobs.append(fobs)
          self.sigmas.append(sigma)

  def redefine_unit_cell(self, unit_cell):
    self.unit_cell = unit_cell

  def redefine_space_group(self, space_group_info):
    self.space_group_info = space_group_info

  def as_miller_array(self, info=0):
    assert self.space_group_info != None
    crystal_symmetry = crystal.symmetry(
      self.unit_cell, space_group_info=self.space_group_info)
    miller_set = miller.set(
      crystal_symmetry, self.miller_indices, self.anomalous)
    return miller.array(miller_set, self.fobs, self.sigmas, info=info)

def run():
  import sys, os
  import cPickle
  import pickle
  binary_pickle = 1
  to_pickle = "--pickle" in sys.argv[1:]
  for file_name in sys.argv[1:]:
    if (file_name.startswith("--")): continue
    f = open(file_name, "r")
    s = scalepack_reader(f)
    f.close()
    miller_array = s.as_miller_array(info="From file: "+file_name)
    miller_array.show_summary()
    if (to_pickle):
      pickle_file_name = file_name + ".pickle"
      f = open(pickle_file_name, "wb")
      pickle.dump(miller_array, f, binary_pickle)
      f.close()
      f = open(pickle_file_name, "rb")
      cPickle.load(f)
      f.close()
    print
  t = os.times()
  print "u+s,u,s: %.2f %.2f %.2f" % (t[0] + t[1], t[0], t[1])

if (__name__ == "__main__"):
  run()
