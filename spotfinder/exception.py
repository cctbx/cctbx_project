from __future__ import absolute_import, division, print_function

class SpotfinderError(Exception):
  def __init__(self,message,processdict=None):
    Exception.__init__(self,message)
    self.classname="Spotfinder Problem"
    self.parameters = processdict
