from __future__ import absolute_import, division, print_function
import boost_adaptbx.boost.python as bp
from six.moves import zip
ext = bp.import_ext("cctbx_xray_ext")
from cctbx_xray_ext import *

from cctbx.array_family import flex
import sys

@bp.inject_into(scattering_type_registry)
class _():

  def type_count_dict(self):
    result = {}
    unique_counts = list(self.unique_counts)
    for t,i in self.type_index_pairs_as_dict().items():
      result[t] = unique_counts[i]
    return result

  def sorted_type_index_pairs(self, heaviest_first=True):
    ugs = self.unique_gaussians_as_list()
    pairs = []
    sf0s = flex.double()
    for t,i in self.type_index_pairs_as_dict().items():
      pairs.append((t,i))
      gaussian = ugs[i]
      if (gaussian is None):
        sf0s.append(0)
      else:
        sf0s.append(gaussian.at_stol(0))
    perm = flex.sort_permutation(sf0s, reverse=heaviest_first)
    return flex.select(pairs, permutation=perm)

  def show_summary(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    unique_gaussians = self.unique_gaussians_as_list()
    unique_counts = list(self.unique_counts)
    parts = []
    for t,i in self.sorted_type_index_pairs():
      gaussian = unique_gaussians[i]
      if (gaussian is None):
        gn = "None"
      else:
        gn = str(gaussian.n_terms())
        if (gaussian.c() != 0):
          gn += "+c"
      parts.append("%s%s:%s*%d" % (prefix, t, gn, unique_counts[i]))
      prefix = ""
    print(" ".join(parts), file=out)

  def warning_if_any(self,std_lbl):
    if( self.last_table()=="neutron" ):
      from cctbx.eltbx.neutron import neutron_news_1992_table
      scattering_info = neutron_news_1992_table(std_lbl, True)
      b = scattering_info.bound_coh_scatt_length()
      if(b.imag != 0.0):
        return " WARNING! anomalous neutron scattering ("+str(b.imag) \
            +") is ignored!"
    return ""

  def show(self,
        header="Number of scattering types:",
        show_sf0=True,
        show_gaussians=True,
        out=None,
        prefix=""):
    if (out is None): out = sys.stdout
    unique_gaussians = self.unique_gaussians_as_list()
    unique_counts = list(self.unique_counts)
    tips = self.sorted_type_index_pairs()
    if (header is not None):
      print(prefix + header, len(tips), file=out)
    if (len(tips) == 0):
      print(prefix + "  Empty scattering-type registry.", file=out)
    else:
      nt = max(3,max([len(t) for t,i in tips]))
      nt_fmt = "%%-%ds " % nt
      nc = max(5,len(str(max(unique_counts))))
      nc_fmt = "%%%dd" % nc
      line = prefix + "  Type%s %sNumber" % (" "*(nt-3), " "*(nc-5))
      if (show_sf0): line += "    sf(0)"
      if (show_gaussians): line += "   Gaussians"
      print(line, file=out)
      for t,i in tips:
        line = prefix + "   " \
             + nt_fmt%t \
             + nc_fmt%unique_counts[i] + " "
        gaussian = unique_gaussians[i]
        if (show_sf0):
          if (gaussian is None):
            line += "     None"
          else:
            line += " %8.2f" % gaussian.at_stol(0)
        if (show_gaussians):
          if (gaussian is None):
            line += "      None"
          else:
            line += " %7s" % str(gaussian.n_terms())
            if (gaussian.c() != 0): line += "+c"
        print("%s%s" % (line.rstrip(), self.warning_if_any(t)), file=out)
      if (show_sf0):
        print(prefix \
          + "  sf(0) = scattering factor at diffraction angle 0.", file=out)

  def sum_of_scattering_factors_at_diffraction_angle_0(self):
    result = 0
    for g,c in zip(self.unique_gaussians_as_list(), self.unique_counts):
      result += g.at_stol(0) * c
    return result

  def wilson_dict(self):
    result = {}
    unique_counts = list(self.unique_counts)
    for t,i in self.sorted_type_index_pairs():
      result[t] = unique_counts[i]
    return result

  def as_type_gaussian_dict(self):
    result = {}
    ugs = self.unique_gaussians_as_list()
    for t,i in self.type_index_pairs_as_dict().items():
      result[t] = ugs[i]
    return result

@bp.inject_into(sampled_model_density)
class _():

  def real_map_unpadded(self):
    from cctbx import maptbx
    result = self.real_map()
    if (not result.is_padded()): return result
    return maptbx.copy(result, flex.grid(result.focus()))
