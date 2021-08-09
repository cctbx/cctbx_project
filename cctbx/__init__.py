from __future__ import absolute_import, division, print_function
from libtbx.version import get_version

__version__ = get_version()

factor_kev_angstrom = 6.62607015 * 2.99792458 / 1.602176634
factor_ev_angstrom  = factor_kev_angstrom * 1000

# (ab)use miller extension to provide hendrickson_lattman constructors
_as_hendrickson_lattman = None
def hendrickson_lattman(*args, **kw):
  global _as_hendrickson_lattman
  if (_as_hendrickson_lattman is None):
    import cctbx.miller
    _as_hendrickson_lattman = cctbx.miller.as_hendrickson_lattman
  return _as_hendrickson_lattman(*args, **kw)
