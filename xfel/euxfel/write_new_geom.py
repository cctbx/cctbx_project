from __future__ import absolute_import, division, print_function
from dxtbx.datablock import DataBlockFactory, DataBlockDumper
import sys
from six.moves import range

fin = sys.argv[1]
px_size = 0.05 # in mm
dpx = [0,0,1,1,0,1,1,0] # in px
dpy = [0,0,1,1,0,1,1,0] # in px

db=DataBlockFactory.from_json_file(fin)
db0=db[0]
dd=db0.to_dict()
for i in range(8):
  x = dd['detector'][0]['panels'][i]['origin'][0]+dpx[i]*px_size
  y = dd['detector'][0]['panels'][i]['origin'][1]+dpy[i]*px_size
  dd['detector'][0]['panels'][i]['origin']=(x,y,0)
xx=DataBlockFactory.from_dict(dd)
yy=DataBlockDumper(xx)
yy.as_file('new_detector_geom.json')
