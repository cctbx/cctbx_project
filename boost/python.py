from __future__ import absolute_import, division, print_function
from boost_adaptbx.boost.python import *
import warnings

warnings.warn(
  "importing from boost.python is deprecated; this module will be removed shortly. "
  "import from boost_adaptbx.boost.python instead. "
  "Please see https://github.com/cctbx/cctbx_project/issues/458 for more information.",
  FutureWarning,
  stacklevel=2
)

_skip_warning = True
class injector(object):
  '''Deprecated function. Instead of

  class CrystalExt(boost.python.injector, Crystal):

  please use

  @boost_adaptbx.boost.python.inject_into(Crystal)
  class _(object):
  '''
  class __metaclass__(meta_class):

    def __init__(self, name, bases, dict):
      if not _skip_warning:
        warnings.warn(
          "boost.python.injector is deprecated and does not work on Python 3. "
          "Please see https://github.com/cctbx/cctbx_project/pull/386 for more information.",
          FutureWarning,
          stacklevel=2
        )
      if (len(bases) > 1):
        # bases[0] is this injector
        target = bases[1] # usually a Boost.Python class
        def setattr_from_dict(d):
          for k,v in d.items():
            if (k not in (
                  "__init__",
                  "__del__",
                  "__module__",
                  "__file__",
                  "__dict__")):
              if (k != "__doc__" or v is not None):
                setattr(target,k,v)
        setattr_from_dict(dict)
        for b in bases[2:]: # usually mix-in classes, if any
          setattr_from_dict(b.__dict__)
      return type.__init__(self, name, (), {})
_skip_warning = False
