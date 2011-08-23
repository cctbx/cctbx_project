factor_kev_angstrom = 6.6260755 * 2.99792458 / 1.60217733
factor_ev_angstrom  = 6626.0755 * 2.99792458 / 1.60217733

# (ab)use miller extension to provide hendrickson_lattman constructors
_as_hendrickson_lattman = None
def hendrickson_lattman(*args, **kw):
  global _as_hendrickson_lattman
  if (_as_hendrickson_lattman is None):
    import cctbx.miller
    _as_hendrickson_lattman = cctbx.miller.as_hendrickson_lattman
  return _as_hendrickson_lattman(*args, **kw)
