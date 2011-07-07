'''
Determined 8/25/2004 for Q210 images collected at beamline 8.2.1.

The "film coordinates" needed for autoindexing are defined as follows:

MOSFLM, HKL2000 and LABELIT.  The origin is the top corner of the image
plate farthest from the storage ring, that is the top left corner when
viewed from the source.  X points down, Y points towards the storage
ring.

ADXV.  The origin is the bottom corner of the image plate farthest from
the storage ring, that is the bottom left corner when viewed from the
source.  X points toward the storage ring, Y points up.


The GUI view for each program is as follows:

HKL2000, LABELIT and ADXV.  The image emulates the view of an observer
standing at the source, i.e., with the storage ring closest to the right
edge of the image.

MOSFLM.  The image emulates the view of an observer looking toward the
source, and in addition, the image is rotated 90 clockwise.  Therefore,
the storage ring is nearest the top edge of the image, and the top edge
of the detector is at the right edge of the image.
'''

class spot_xy_convention:
  def __init__(self,W1,W2):
    self.W1 = W1; self.W2=W2

  def select(self,spot,index):
    if index%2==1: assert self.W1==self.W2
    if index==0:  return (spot[0],spot[1],spot[2])
    if index==1:  return (spot[1],spot[0],spot[2])
    if index==2:  return (self.W1-spot[0],spot[1],spot[2])
    if index==3:  return (spot[1],self.W1-spot[0],spot[2])
    if index==4:  return (spot[0],self.W2-spot[1],spot[2])
    if index==5:  return (self.W1-spot[1],spot[0],spot[2])
    if index==6:  return (self.W1-spot[0],self.W2-spot[1],spot[2])
    if index==7:  return (self.W1-spot[1],self.W1-spot[0],spot[2])
    if index==8:  return (spot[0],spot[1],-spot[2])
    if index==9:  return (spot[1],spot[0],-spot[2])
    if index==10:  return (self.W1-spot[0],spot[1],-spot[2])
    if index==11:  return (spot[1],self.W1-spot[0],-spot[2])
    if index==12:  return (spot[0],self.W2-spot[1],-spot[2])
    if index==13:  return (self.W1-spot[1],spot[0],-spot[2])
    if index==14:  return (self.W1-spot[0],self.W2-spot[1],-spot[2])
    if index==15:  return (self.W1-spot[1],self.W1-spot[0],-spot[2])
    raise

  def inverse(self,tspot):
    if self.index%2==1: assert self.W1==self.W2
    if self.index==0:  return (tspot[0],tspot[1],tspot[2])
    if self.index==1:  return (tspot[1],tspot[0],tspot[2])
    if self.index==2:  return (self.W1-tspot[0],tspot[1],tspot[2])
    if self.index==3:  return (tspot[1],self.W1-tspot[0],tspot[2])
    if self.index==4:  return (tspot[0],self.W2-tspot[1],tspot[2])
    if self.index==5:  return (self.W1-tspot[1],tspot[0],tspot[2])
    if self.index==6:  return (self.W1-tspot[0],self.W2-tspot[1],tspot[2])
    if self.index==7:  return (self.W1-tspot[1],self.W1-tspot[0],tspot[2])
    if self.index==8:  return (tspot[0],tspot[1],-tspot[2])
    if self.index==9:  return (tspot[1],tspot[0],-tspot[2])
    if self.index==10:  return (self.W1-tspot[0],tspot[1],-tspot[2])
    if self.index==11:  return (tspot[1],self.W1-tspot[0],-tspot[2])
    if self.index==12:  return (tspot[0],self.W2-tspot[1],-tspot[2])
    if self.index==13:  return (self.W1-tspot[1],tspot[0],-tspot[2])
    if self.index==14:  return (self.W1-tspot[0],self.W2-tspot[1],-tspot[2])
    if self.index==15:  return (self.W1-tspot[1],self.W1-tspot[0],-tspot[2])
    raise

