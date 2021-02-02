from __future__ import absolute_import, division, print_function
from six.moves import range
import struct

header_struct = [
  ('Device',10,'s'),
  ('Version',10,'s'),
  ('Crystal',20,'s'),
  ('CrystalSystem',12,'s'),
  (None,24),
  ('SpaceGroup',12,'s'),
  ('mosaic1',4,'!f'),
  ('memo',80,'s'),
  ('reserve1',84,'s'),
  ('date',12,'s'),
  ('operatorname',20,'s'),
  ('target',4,'s'),
  ('wavelength',4,'!f'),
  ('monotype',20,'s'),
  ('mono2theta',4,'!f'),
  ('collimator',20,'s'),
  ('filter',4,'s'),
  ('distance',4,'!f'),
  ('Kv',4,'!f'),
  ('mA',4,'!f'),
  ('focus',12,'s'),
  ('Xmemo',80,'s'),
  ('cyl',4,'!i'),
  (None,60),
  ('Spindle',4,'s'),          # Crystal mount axis closest to spindle axis
  ('Xray_axis',4,'s'),        # Crystal mount axis closest to beam axis
  ('phidatum',4,'!f'),
  ('phistart',4,'!f'),
  ('phiend',4,'!f'),
  ('noscillations',4,'!i'),
  ('minutes',4,'!f'),         # Exposure time in minutes?
  ('beampixels_x',4,'!f'),
  ('beampixels_y',4,'!f'),    # Direct beam position in pixels
  ('omega',4,'!f'),
  ('chi',4,'!f'),
  ('twotheta',4,'!f'),
  ('Mu',4,'!f'),              # Spindle inclination angle?
  ('ScanTemplate',204,'s'),   # This space is now used for storing the scan
                              # templates information
  ('nFast',4,'!i'),
  ('nSlow',4,'!i'),           # Number of fast, slow pixels
  ('sizeFast',4,'!f'),
  ('sizeSlow',4,'!f'),        # Size of fast, slow direction in mm
  ('record_length',4,'!i'),   # Record length in bytes
  ('number_records',4,'!i'),  # number of records
  ('Read_start',4,'!i'),      # For partial reads, 1st read line
  ('IP_num',4,'!i'),          # Which imaging plate 1, 2 ?
  ('Ratio',4,'!f'),           # Output ratio for high value pixels
  ('Fading_start',4,'!f'),    # Fading time to start of read
  ('Fading_end',4,'!f'),      # Fading time to end of read
  ('computer',10,'s'),        # Type of computer "IRIS", "VAX", "SUN", etc
  ('plate_type',10,'s'),      # Type of IP
  ('Dr',4,'!i'),
  ('Dx',4,'!i'),
  ('Dz',4,'!i'),              # IP scanning codes??
  ('PixShiftOdd',4,'!f'),     # Pixel shift to odd lines
  ('IntRatioOdd',4,'!f'),     # Intensity ratio to odd lines
  ('MagicNum',4,'!i'),        # Magic number to indicate next values are legit
  ('NumGonAxes',4,'!i'),      # Number of goniometer axes
  ('a5x3fGonVecs',60,'!fffffffffffffff'),# Goniometer axis vectors
  ('a5fGonStart',20,'!fffff'),# Start angles for each of 5 axes
  ('a5fGonEnd',20,'!fffff'),  # End angles for each of 5 axes
  ('a5fGonOffset',20,'!fffff'),# Offset values for each of 5 axes
  ('ScanAxisNum',4,'!i'),     # Which axis is the scan axis?
  ('AxesNames',40,'s'),       # Names of the axes (space or comma separated?)'''
]
class Raxis(object):
  def __init__(self,file):
    self.file = file

  def readHeader(self,verbose=0):
    with open(self.file,'rb') as F:
      self.head={}
      seek = 0
      for item in header_struct:
        if item[0]==None:
          F.read(item[1])
        elif item[2]=='s':
          self.head[item[0]]=F.read(item[1])[0:item[1]]
          if verbose:print(item[0],self.head[item[0]])
        elif len(item[2])>2:
          rawdata = F.read(item[1])
          assert len(rawdata)==struct.calcsize(item[2])
          self.head[item[0]] = struct.unpack(item[2],rawdata)
          if verbose:print(item[0],self.head[item[0]])
        else:
          rawdata = F.read(item[1])
          assert len(rawdata)==struct.calcsize(item[2])
          self.head[item[0]] = struct.unpack(item[2],rawdata)[0]
          if verbose:print(item[0],self.head[item[0]])
        seek+=item[1]

  def data(self):
    Dim0 = self.head['nFast'] #number of fast pixels
    ToRead = self.head['record_length']
    ReadLines = self.head['number_records']

    with open(self.file,'rb') as F:
      F.seek(ToRead)
      raw_data = F.read(ToRead * ReadLines)

    # For a normal image, there should be no padding per line
    # Each line might be padded, so figure this out
    BytesPerLine = Dim0 * 2;
    if BytesPerLine < ToRead:
      # Remove all padding bytes
      raw_data = b"".join(
        raw_data[record * ToRead : record * ToRead + BytesPerLine]
        for record in range(ReadLines)
      )

    self.CharTemp = raw_data

  def dump(self):
    ptr = 0
    for x in range(0,len(CharTemp),2):
      unsigned_int = struct.unpack( "!H",self.CharTemp[x:x+2] )[0]
      if unsigned_int <= 32767:
        print(float(unsigned_int))
      else:
        print(( float(unsigned_int)+32768.0 ) * self.head['Ratio'])

if __name__=='__main__':
  R = Raxis('H-x071_0001.osc')
  R.readHeader()
  R.data()
  R.dump()
