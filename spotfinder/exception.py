from __future__ import division
from exceptions import Exception

class SpotfinderError(Exception):
  def __init__(self,message,processdict=None):
    Exception.__init__(self,message)
    self.classname="Spotfinder Problem"
    self.parameters = processdict
