from __future__ import division
from dxtbx_model_ext import *
from beam import *
from goniometer import *
from detector import *
from scan import *

# definition of __getinitargs__ for Detector class
def detector_getinitargs(self):
    panel_list = [panel for panel in self]
    return (panel_list,)

# now inject __getinitargs__
Detector.__getinitargs__ = detector_getinitargs
