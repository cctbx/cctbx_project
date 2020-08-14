from __future__ import absolute_import, division, print_function
from cctbx.xray import ext
import boost_adaptbx.boost.python as bp
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

@bp.inject_into(ext.gradient_flags)
class _():

  def copy(self):
    return ext.gradient_flags(self)

  def show_summary(self, f=None):
    if (f is None): f = sys.stdout
    print("gradient_flags:", file=f)
    print("  site:", self.site, file=f)
    print("  u_iso:", self.u_iso, file=f)
    print("  u_aniso:", self.u_aniso, file=f)
    print("  occupancy:", self.occupancy, file=f)
    print("  fp:", self.fp, file=f)
    print("  fdp:", self.fdp, file=f)
    print("  sqrt_u_iso:", self.sqrt_u_iso, file=f)
    print("  tan_b_iso_max:", self.tan_b_iso_max, file=f)

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
