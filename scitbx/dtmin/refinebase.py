from __future__ import division
from scitbx.dtmin.compulsory import Compulsory
from scitbx.dtmin.optional import Optional
from scitbx.dtmin.min_logging import Logging
from scitbx.dtmin.auxiliary import Auxiliary

class RefineBase(Compulsory, Optional, Logging, Auxiliary):

  def __init__(self):
    self.nmp = 0 # number of macrocycle parameters
