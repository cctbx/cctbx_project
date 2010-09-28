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
        fdp=None,
        flags=None,
        u_iso=None,
        u_star=None):
    assert u is None or b is None
    if (b is not None): u = adptbx.b_as_u(b)
    del b
    if (u is not None):
      assert u_iso is None and u_star is None
      if (flags is None):
        assert self.flags.use_u_iso() ^ self.flags.use_u_aniso()
      else:
        assert flags.use_u_iso() ^ flags.use_u_aniso()
    result = ext.scatterer(other=self)
    if (label is not None): result.label = label
    if (site is not None): result.site = site
    if (flags is not None): result.flags = flags
    if (u is not None):
      if (result.flags.use_u_iso()):
        result.u_iso = u
      else:
        result.u_star = u
    if (u_iso is not None): result.u_iso = u_iso
    if (u_star is not None): result.u_star = u_star
    if (occupancy is not None): result.occupancy = occupancy
    if (scattering_type is not None): result.scattering_type = scattering_type
    if (fp is not None): result.fp = fp
    if (fdp is not None): result.fdp = fdp
    return result

  def __repr__(self):
    """ The returned string does usually eval to self, except if the both of
    self.flags.use_u_iso() and self.flags.use_u_aniso() are True.
    """
    r = []
    r.append("  label='%s'" % self.label)
    if (eltbx.xray_scattering.get_standard_label(label=self.label,
                                                 exact=False)
        != self.scattering_type):
      r.append("  scattering_type='%s'" % self.scattering_type)
    r.append("  site=(%.6f, %.6f, %.6f)"  % self.site)
    if self.flags.use_u_iso():
      r.append("  u=%.6f" % self.u_iso)
    if self.flags.use_u_aniso():
      r.append("  u=(%.6f, %.6f, %.6f,\n"
               "     %.6f, %.6f, %.6f)" % self.u_star)
    if self.occupancy != 1:
      r.append("  occupancy=%.6f" % self.occupancy)
    if self.fp != 0 or self.fdp != 0:
      r.append("  fp=%.6f" % self.fp)
      r.append("  fdp=%.6f" % self.fdp)
    return "xray.scatterer(\n%s)" % ",\n".join(r)

  def element_and_charge_symbols(self, exact=False):
    return eltbx.xray_scattering.get_element_and_charge_symbols(
      scattering_type=self.scattering_type, exact=exact)

  def element_symbol(self, exact=False):
    e, c = self.element_and_charge_symbols(exact=exact)
    if (len(e) == 0): return None
    return e

  def show(self, f=None, unit_cell=None):
    if (f is None): f = sys.stdout
    print >> f, "%-4s" % self.label,
    print >> f, "%-4s" % self.scattering_type,
    print >> f, "%3d" % self.multiplicity(),
    print >> f, "(%7.4f %7.4f %7.4f)" % self.site,
    print >> f, "%4.2f" % self.occupancy,
    if self.flags.use_u_iso():
      print >> f, "%6.4f" % self.u_iso,
    else:
      print >> f, '[ - ]',
    if self.flags.use_u_aniso():
      assert unit_cell is not None
      u_cart = adptbx.u_star_as_u_cart(unit_cell, self.u_star)
      print >> f, "%6.4f" % adptbx.u_cart_as_u_iso(u_cart)
      print >> f, "     u_cart =", ("%6.3f " * 5 + "%6.3f") % u_cart,
    else:
      print >> f, '[ - ]',
    if (self.fp != 0 or self.fdp != 0):
      print >> f, "\n     fp,fdp = %6.4f,%6.4f" % (
        self.fp,
        self.fdp),
    print >> f

class anomalous_scatterer_group:

  def __init__(self,
        iselection,
        f_prime,
        f_double_prime,
        refine,
        selection_string=None):
    self.iselection = iselection
    self.f_prime = f_prime
    self.f_double_prime = f_double_prime
    self.refine_f_prime = False
    self.refine_f_double_prime = False
    for refine_item in refine:
      assert refine_item in ["f_prime", "f_double_prime"]
      if (refine_item == "f_prime"):
        self.refine_f_prime = True
      else:
        self.refine_f_double_prime = True
    self.selection_string = selection_string

  def labels_refine(self):
    result = []
    if (self.refine_f_prime): result.append("f_prime")
    if (self.refine_f_double_prime): result.append("f_double_prime")
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
    labels_refine = self.labels_refine()
    print >> out, prefix+"  refine:",
    if (len(labels_refine) == 0): print >> out, "None"
    else: print >> out, " ".join(labels_refine)

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
