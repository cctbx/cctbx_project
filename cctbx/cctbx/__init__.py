# (ab)use miller extension to provide hendrickson_lattman constructors
_as_hendrickson_lattman = None
def hendrickson_lattman(*args, **kw):
  global _as_hendrickson_lattman
  if (_as_hendrickson_lattman is None):
    import cctbx.miller
    _as_hendrickson_lattman = cctbx.miller.as_hendrickson_lattman
  return _as_hendrickson_lattman(*args, **kw)
