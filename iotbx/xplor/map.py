"""Example format for Xplor Maps:

       2 !NTITLE
 REMARKS FILENAME=""
 REMARKS scitbx.flex.double to Xplor map format
      24       0      24     120       0     120      54       0      54
 3.20420E+01 1.75362E+02 7.96630E+01 9.00000E+01 9.00000E+01 9.00000E+01
ZYX
       0
-2.84546E-01-1.67775E-01-5.66095E-01-1.18305E+00-1.49559E+00-1.31942E+00
-1.01611E+00-1.00873E+00-1.18992E+00-1.02460E+00-2.72099E-01 5.94242E-01
<deleted>
   -9999
  0.0000E+00  1.0000E+00
That is:
...a blank line
...an integer giving the number of title lines, with mandatory !NTITLE
...title lines in %-264s format
...X, Y, and Z sections giving:
       sections per unit cell, in the given direction
       ordinal of first section in file
       ordinal of last section in file
...unit cell dimensions
...slow, medium, fast section order, always ZYX
...for each slow section, the section number
...sectional data in special fortran format shown
...-9999
...map average and standard deviation
"""
from __future__ import absolute_import, division, print_function

import iotbx.xplor.ext as ext
from cctbx import miller
from cctbx import maptbx
from cctbx import uctbx
from cctbx.array_family import flex
import sys
from six.moves import range
from six.moves import zip

class gridding(object):

  def __init__(self, n, first, last):
    self.n = tuple(n)
    self.first = tuple(first)
    self.last = tuple(last)
    for x,ni,fi,li in zip("XYZ",n,first,last):
      if (ni < 1 or fi > li): raise RuntimeError(
        "Illegal xplor map gridding for dimension %s: "
        "gridding=%d, first=%d, last=%d" % (x,ni,fi,li))

  def format_9i8(self):
    result = ""
    for triple in zip(self.n,self.first,self.last):
      result += ("%8d"*3) % triple
    return result

  def as_flex_grid(self):
    return flex.grid(self.first, self.last, False)

  def is_compatible_flex_grid(self, flex_grid, is_p1_cell=False):
    if (flex_grid.nd() != 3): return False
    if (is_p1_cell is None): # XXX temporary flag allowing any cell
      self_as_flex_grid = self.as_flex_grid()
      if (flex_grid.origin() != self_as_flex_grid.origin()): return False
      if (flex_grid.last() != self_as_flex_grid.last()): return False
    elif (not is_p1_cell):
      self_as_flex_grid = self.as_flex_grid()
      if (flex_grid.origin() != self_as_flex_grid.origin()): return False
      if (flex_grid.last() != self_as_flex_grid.last()): return False
    else:
      if (not flex_grid.is_0_based()): return False
      if (flex_grid.focus() != self.n): return False
    return True

class reader(object):

  def __init__(self, file_name, header_only=False):
    with open(file_name, "r") as f:
      f.readline()
      self.title_lines = []
      ntitle = int(f.readline().strip().split("!")[0])
      self.title_lines=[]
      for x in range(ntitle):
        line = f.readline().rstrip()
        self.title_lines.append(line)
      line = f.readline()
      values = [int(line[i:i+8]) for i in range(0,72,8)]
      self.gridding = gridding(
        n     = [values[i] for i in range(0,9,3)],
        first = [values[i] for i in range(1,9,3)],
        last  = [values[i] for i in range(2,9,3)])
      line = f.readline()
      params = [float(line[i:i+12]) for i in range(0,72,12)]
      self.unit_cell = uctbx.unit_cell(params)
      order = f.readline().strip()
      assert order == "ZYX"

    if (header_only):
      self.data = None
      self.average = None
      self.standard_deviation = None
    else:
      ext_reader = ext.map_reader(
        file_name=file_name,
        n_header_lines=len(self.title_lines)+5,
        grid=self.gridding.as_flex_grid())
      self.data = ext_reader.data
      self.average = ext_reader.average
      self.standard_deviation = ext_reader.standard_deviation

  def show_summary(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    print(prefix+"Title lines:", len(self.title_lines), file=out)
    for line in self.title_lines:
      print(prefix+"  "+line.rstrip(), file=out)
    g = self.gridding
    print(prefix+"Gridding:", file=out)
    print(prefix+"  n:    ", g.n, file=out)
    print(prefix+"  first:", g.first, file=out)
    print(prefix+"  last: ", g.last, file=out)
    print(prefix+"Total number of data points:", self.data.size(), file=out)
    stats = maptbx.statistics(self.data)
    print(prefix+"  min:   %.6g" % stats.min(), file=out)
    print(prefix+"  max:   %.6g" % stats.max(), file=out)
    print(prefix+"  mean:  %.6g" % stats.mean(), file=out)
    print(prefix+"  sigma: %.6g" % stats.sigma(), file=out)

def writer(file_name, title_lines, unit_cell, gridding,
           data, is_p1_cell=False,
           average=-1,
           standard_deviation=-1):
  assert gridding.is_compatible_flex_grid(
    flex_grid=data.accessor(),
    is_p1_cell=is_p1_cell)
  with open(file_name, "w") as f:
    f.write("\n")
    f.write("%8d !NTITLE\n" % len(title_lines))
    for line in title_lines:
      f.write("%-264s\n" % line)
    f.write("%s\n" % gridding.format_9i8())

  if (is_p1_cell is None): # XXX temporary flag allowing any cell
    ext.map_writer(
      file_name=file_name,
      unit_cell=unit_cell,
      data=data,
      average=average,
      standard_deviation=standard_deviation)
  elif (not is_p1_cell):
    ext.map_writer(
      file_name=file_name,
      unit_cell=unit_cell,
      data=data,
      average=average,
      standard_deviation=standard_deviation)
  else:
    ext.map_writer(
      file_name=file_name,
      unit_cell=unit_cell,
      gridding_first=gridding.first,
      gridding_last=gridding.last,
      data=data,
      average=average,
      standard_deviation=standard_deviation)

def cctbx_miller_fft_map_as_xplor_map(
      self,
      file_name,
      title_lines=["cctbx.miller.fft_map"],
      gridding_first=None,
      gridding_last=None,
      average=None,
      standard_deviation=None):
  if (gridding_first is None): gridding_first = (0,0,0)
  if (gridding_last is None): gridding_last = self.n_real()
  gridding_ = gridding(
    n=self.n_real(),
    first=gridding_first,
    last=gridding_last)
  data = self.real_map()
  if (average is None or standard_deviation is None):
    statistics = maptbx.statistics(data)
    if (average is None): average = statistics.mean()
    if (standard_deviation is None): standard_deviation = statistics.sigma()
  writer(
    file_name=file_name,
    title_lines=title_lines,
    unit_cell=self.unit_cell(),
    gridding=gridding_,
    data=data,
    is_p1_cell=True,
    average=average,
    standard_deviation=standard_deviation)

# injecting
miller.fft_map.as_xplor_map = cctbx_miller_fft_map_as_xplor_map
