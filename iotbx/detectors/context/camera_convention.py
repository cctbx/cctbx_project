from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex

"Transformation specifically meant for use by spotfinder & labelit"
def Cameras(index):
  '''All ADSC detectors: index=1
         MarCCD:         index=1
  '''
  from cctbx.array_family import flex
  if index==0:  return flex.double((1,0,0,1))
  if index==1:  return flex.double((-1,0,0,-1))
  if index==2:  return flex.double((-1,0,0,1))
  if index==3:  return flex.double((1,0,0,-1))
  if index==4:  return flex.double((0,1,1,0))
  if index==5:  return flex.double((0,-1,-1,0))
  if index==6:  return flex.double((0,1,-1,0))
  if index==7:  return flex.double((0,-1,1,0))
