from cctbx.xray import ext
import boost.python
import sys

class gradient_flags(ext.gradient_flags):

  def __init__(self, site=None,
                     u_iso=None,
                     u_aniso=None,
                     occupancy=None,
                     fp=None,
                     fdp=None,
                     u=None,
                     sqrt_u_iso=False,
                     tan_b_iso_max=False,
                     default=False):
    if (u is not None): assert u_iso is None and u_aniso is None
    if (u is None): u = default
    if (site is None): site = default
    if (u_iso is None): u_iso = u
    if (u_aniso is None): u_aniso = u
    if (occupancy is None): occupancy = default
    if (fp is None): fp = default
    if (fdp is None): fdp = default
    ext.gradient_flags.__init__(self,
      site, u_iso, u_aniso, occupancy, fp, fdp, sqrt_u_iso, tan_b_iso_max)

class _(boost.python.injector, ext.gradient_flags):

  def copy(self):
    return ext.gradient_flags(self)

  def show_summary(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "gradient_flags:"
    print >> f, "  site:", self.site
    print >> f, "  u_iso:", self.u_iso
    print >> f, "  u_aniso:", self.u_aniso
    print >> f, "  occupancy:", self.occupancy
    print >> f, "  fp:", self.fp
    print >> f, "  fdp:", self.fdp
    print >> f, "  sqrt_u_iso:", self.sqrt_u_iso
    print >> f, "  tan_b_iso_max:", self.tan_b_iso_max

  def string_of_true(self):
    result = []
    if (self.site): result.append("site")
    if (self.u_iso): result.append("u_iso")
    if (self.u_aniso): result.append("u_aniso")
    if (self.occupancy): result.append("occupancy")
    if (self.fp): result.append("fp")
    if (self.fdp): result.append("fdp")
    if (self.sqrt_u_iso): result.append("sqrt_u_iso")
    if (self.tan_b_iso_max): result.append("tan_b_iso_max")
    return ",".join(result)
