#!/usr/bin/env python
# FormatStill.py
#
# Root class for still shots.  A still shot has no goniomter and no
# scan in their model, as these constructs are not meaningful.
#
from __future__ import absolute_import, division
from dxtbx.format.Format import Format
from dxtbx.model.detector import Detector
from dxtbx.model.beam import Beam
import exceptions

class FormatStill(Format):
  def setup(self):
    '''Read the image file, construct the information which we will be
    wanting about the experiment from this. N.B. in your implementation
    of this you will probably want to make use of the static methods
    below and probably add some format parsing code too. Please also keep
    in mind that your implementation may be further subclassed by
    someone else.

    Do not create scan or or goniometer objects'''
    self._start()
    try:
      detector_instance = self._detector()
      assert(isinstance(detector_instance, Detector))
      self._detector_instance = detector_instance

      beam_instance = self._beam()
      assert(isinstance(beam_instance, Beam))
      self._beam_instance = beam_instance

    except exceptions.Exception:
      # FIXME ideally should not squash the errors here...
      pass
    finally:
      self._end()

  def _goniometer(self):
    '''Not sensible for still shot data'''

    return None


  def _scan(self):
    '''Not sensible for still shot data'''

    return None
