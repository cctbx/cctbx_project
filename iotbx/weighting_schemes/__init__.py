""" Weighting scheme mixing to be inherited from to provide I/O behaviours,
and actual weighting schemes for builders to store the parsed information.
They all feature:

  o __str__: a string representation of the weighting scheme in a format
             that is appropriate for the CIF item _refine_ls_weighting_details;

  o type: the kind of weighting scheme,
          that is appropriate for the CIF item _refine_ls_weighting_scheme
"""

from __future__ import division

import libtbx

class mainstream_shelx_weighting_mixin(object):

  def __str__(self):
    """ A string representation of the weighting scheme in a format that is
        appropriate for the CIF item _refine_ls_weighting_details.
    """
    if round(self.a, 4) in (0.1, 0.2):
      a = "%.1f" %self.a
    else:
      a = "%.4f" %self.a
    if round(self.b, 4) == 0: b_part=""
    else: b_part = "+%.4fP" %self.b
    return ("w=1/[\s^2^(Fo^2^)+(%sP)^2^%s]"
            " where P=(Fo^2^+2Fc^2^)/3" %(a, b_part))

  def type(self):
    return "calc"


class mainstream_shelx_weighting(mainstream_shelx_weighting_mixin):

  def __init__(self, a, b, c=0, d=0, e=0, f=1/3):
    libtbx.adopt_init_args(self, locals())


class unit_weighting_mixin(object):

  def __str__(self):
    return "w=1"

  def type(self):
    return "unit"


class unit_weighting(unit_weighting_mixin): pass


class sigma_weighting_mixin(object):

  def __str__(self): return "w=1/sigma^2"

  def type(self): return "sigma"


class sigma_weighting(sigma_weighting_mixin): pass
