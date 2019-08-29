from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry

"""
Jiffy code to patch phil files from v3 to v4 of XFEL gui
"""

def sync_phil(params, unused):
  for scope in unused:
    value = scope.object.extract()[0]
    if scope.path == 'experiment':
      params.facility.lcls.experiment = value
    elif scope.path == 'web.user':
      params.facility.lcls.web.user = value
    elif scope.path == 'web.password':
      params.facility.lcls.web.password = value
    elif scope.path == 'web.enforce80':
      params.facility.lcls.web.enforce80 = value
    elif scope.path == 'web.enforce81':
      params.facility.lcls.web.enforce81 = value
    elif scope.path == 'use_ffb':
      params.facility.lcls.use_ffb = value
    elif scope.path == 'dump_shots':
      params.facility.lcls.dump_shots = value
    else:
      raise Sorry("Cannot understand argument %s=%s"%(scope.path, value))
