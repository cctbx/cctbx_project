from __future__ import absolute_import, division, print_function
from boost_adaptbx.boost.std_pair import *
import warnings

warnings.warn(
  "importing from boost.std_pair is deprecated; this module will be removed shortly. "
  "import from boost_adaptbx.boost.std_pair instead. "
  "Please see https://github.com/cctbx/cctbx_project/issues/458 for more information.",
  FutureWarning,
  stacklevel=2
)
