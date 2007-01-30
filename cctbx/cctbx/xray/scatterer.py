from cctbx.xray import ext
import cctbx.eltbx.xray_scattering
import cctbx.eltbx.tiny_pse
from cctbx import eltbx
from cctbx import adptbx
from cctbx.array_family import flex
from libtbx.str_utils import show_string
import boost.python
from cStringIO import StringIO
import sys

class scatterer(ext.scatterer):

  def __init__(self, label="",
                     site=(0,0,0),
                     u=None,
                     occupancy=1,
                     scattering_type=None,
                     fp=0,
                     fdp=0,
                     b=None):
    assert u is None or b is None
    if   (b is not None): u = adptbx.b_as_u(b)
    elif (u is None): u = 0
    if (scattering_type is None):
      scattering_type = eltbx.xray_scattering.get_standard_label(
        label=label, exact=False)
    ext.scatterer.__init__(
      self, label, site, u, occupancy, scattering_type, fp, fdp)

class _scatterer(boost.python.injector, ext.scatterer):

  def customized_copy(self,
        label=None,
        site=None,
        u=None,
        b=None,
        occupancy=None,
        scattering_type=None,
        fp=None,
        fdp=None):
    assert u is None or b is None
    if (b is not None): u = adptbx.b_as_u(b)
    if (label is None): label = self.label
    if (site is None): site = self.site
    if (u is None):
      if (self.anisotropic_flag): u = self.u_star
      else: u = self.u_iso
    if (occupancy is None): occupancy = self.occupancy
    if (scattering_type is None): scattering_type = self.scattering_type
    if (fp is None): fp = self.fp
    if (fdp is None): fdp = self.fdp
    return scatterer(
      label=label,
      site=site,
      u=u,
      occupancy=occupancy,
      scattering_type=scattering_type,
      fp=fp,
      fdp=fdp)

  def element_symbol(self):
    try:
      return eltbx.tiny_pse.table(self.scattering_type).symbol()
    except RuntimeError:
      return None

  def show(self, f=None, unit_cell=None):
    if (f is None): f = sys.stdout
    print >> f, "%-4s" % self.label,
    print >> f, "%-4s" % self.scattering_type,
    print >> f, "%3d" % self.multiplicity(),
    print >> f, "(%7.4f %7.4f %7.4f)" % self.site,
    print >> f, "%4.2f" % self.occupancy,
    if (not self.anisotropic_flag):
      print >> f, "%6.4f" % self.u_iso
    else:
      assert unit_cell is not None
      u_cart = adptbx.u_star_as_u_cart(unit_cell, self.u_star)
      print >> f, "%6.4f" % adptbx.u_cart_as_u_iso(u_cart)
      print >> f, "     u_cart =", ("%6.3f " * 5 + "%6.3f") % u_cart
    if (self.fp != 0 or self.fdp != 0):
      print >> f, "     fp,fdp = %6.4f,%6.4f" % (
        self.fp,
        self.fdp)

class anomalous_scatterer_group:

  def __init__(self,
        iselection,
        f_prime,
        f_double_prime,
        fix,
        selection_string=None):
    self.iselection = iselection
    self.f_prime = f_prime
    self.f_double_prime = f_double_prime
    self.fix_f_prime = False
    self.fix_f_double_prime = False
    for fix_item in fix:
      assert fix_item in ["f_prime", "f_double_prime"]
      if (fix_item == "f_prime"):
        self.fix_f_prime = True
      else:
        self.fix_f_double_prime = True
    self.selection_string = selection_string

  def labels_fixed(self):
    result = []
    if (self.fix_f_prime): result.append("f_prime")
    if (self.fix_f_double_prime): result.append("f_double_prime")
    return result

  def show_summary(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    print >> out, prefix+"Anomalous scatterer group:"
    if (self.selection_string is not None):
      print >> out, prefix+"  Selection:", show_string(self.selection_string)
    print >> out, prefix+"  Number of selected scatterers:", \
      self.iselection.size()
    print >> out, prefix+"  f_prime:        %.6g" % self.f_prime
    print >> out, prefix+"  f_double_prime: %.6g" % self.f_double_prime
    labels_fixed = self.labels_fixed()
    print >> out, prefix+"  fix:",
    if (len(labels_fixed) == 0): print >> out, "None"
    else: print >> out, " ".join(labels_fixed)

  def copy_to_scatterers_in_place(self, scatterers):
    for i_seq in self.iselection:
      scatterers[i_seq].fp = self.f_prime
      scatterers[i_seq].fdp = self.f_double_prime

  def extract_from_scatterers_in_place(self, scatterers, tolerance=1.e-4):
    fps = flex.double()
    fdps = flex.double()
    fps.reserve(self.iselection.size())
    fdps.reserve(fps.size())
    for i_seq in self.iselection:
      fps.append(scatterers[i_seq].fp)
      fdps.append(scatterers[i_seq].fdp)
      for values,label in [(fps, "f_prime"), (fdps, "f_double_prime")]:
        stats = flex.min_max_mean_double(values)
        if (stats.max - stats.min <= tolerance):
          setattr(self, label, stats.mean)
        else:
          msg = ["Anomalous scatterer group with significantly different %s:"
            % label]
          if (self.selection_string is not None):
            msg.append("  Selection: %s" % show_string(self.selection_string))
          msg.append("  Number of selected scatterers: %d" % stats.n)
          s = StringIO()
          stats.show(out=s, prefix="  %s " % label, show_n=False)
          msg.extend(s.getvalue().splitlines())
          msg.append("  tolerance: %.6g" % tolerance)
          raise RuntimeError("\n".join(msg))
