import cctbx.array_family.flex

import boost.python
ext = boost.python.import_ext("iotbx_xplor_ext")
from iotbx_xplor_ext import *

from cctbx import uctbx

__doc__='''Example format for Xplor Maps:

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
'''

class XplorMap(ext.XplorMap):
  def __init__(self):
    self.title = []
    self.sections = []
    self.unitcell = None
    self.order = None
    self.data = None
    self.average = 0.0
    self.stddev = 1.0
    ext.XplorMap.__init__(self)

  def read(self,arg):
    f = open(arg,"r")
    f.readline()

    # title
    ntitle = int(f.readline().strip().split("!")[0])
    self.title=[]
    for x in xrange(ntitle):
      line = f.readline().rstrip()
      self.title.append(line)

    # sections
    line = f.readline()
    self.sections=[]
    for x in xrange(0,72,8):
      literal = line[x:x+8]
      self.sections.append(int(literal))
    sec = self.sections

    #These assertions imply that the map covers the whole unit
    #cell.  So: get rid of the assertions in short order.
    assert (sec[2]-sec[1] == sec[0] or sec[2]-sec[1] == sec[0]-1)
    assert (sec[5]-sec[4] == sec[3] or sec[5]-sec[4] == sec[3]-1)
    assert (sec[8]-sec[7] == sec[6] or sec[8]-sec[7] == sec[6]-1)

    # unit cell
    line = f.readline()
    ucparams=[]
    for x in xrange(0,72,12):
      literal = line[x:x+12]
      ucparams.append(float(literal))
    self.unitcell = uctbx.unit_cell(tuple(ucparams))

    # ordering information ZYX
    self.order = f.readline().strip()
    f.close()

    self.data = self.ReadXplorMap(arg,len(self.title)+5,self.sections)

    # average and standard deviation
    f = open(arg,"r")
    while (f.readline().strip()!='-9999'):
      pass
    (self.average,self.stddev)=tuple(
      [float(z) for z in f.readline().strip().split()])
    return self

  def write(self,arg):
    assert self.order=='ZYX'
    #but the sections and the data are given as XYZ:
    for x in xrange(3):
      assert self.data.focus()[x]==self.sections[3*x+2]-self.sections[3*x+1]+1

    f = open(arg,'wb')
    f.write("\n")
    f.write("%8d !NTITLE\n"%len(self.title))
    for line in self.title:
      f.write("%-264s\n"%line)
    f.write("%8d%8d%8d%8d%8d%8d%8d%8d%8d\n"%tuple(self.sections))
    f.close()
    self.WriteXplorMap(self.unitcell,self.data,self.average,self.stddev,arg)
