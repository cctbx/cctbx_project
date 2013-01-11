
"""
WARNING: IN DEVELOPMENT AND INCOMPLETE.

Tools for formatting the standard "Table 1" in MX papers.  Relies on other
methods (phenix.refine, phenix.model_vs_data, phenix.merging_statistics, etc.)
to extract statistics for display.
"""

from __future__ import division
from libtbx import slots_getstate_setstate
from libtbx import str_utils
import re

angstrom = u"\u00C5".encode("utf-8", "strict").strip()

# XXX are these complete, and if not, what is missing?  this appears to cover
# the template suggested by NSMB plus more, but not clear what other journals
# expect.  the VTF paper does not list any explicit Table 1 requirements, but
# the forthcoming PDB validation reports may change expectations.
# TODO fill in CIF tags
keywords = [
  # (attr, label, format,  cif_tag)
  # data-only statistics
  ("wavelength", "Wavelength (%s)" % angstrom, "%.4g", None),
  ("d_max_min", "Resolution range (%s)" % angstrom, "%.3g - %.3g", None),
  ("space_group", "Space group", "%s", None),
  ("unit_cell", "Unit cell", "%g %g %g %g %g %g", None),
  ("n_refl_all", "Total reflections", "%d", None),
  ("n_refl", "Unique reflections", "%d", None),
  ("multiplicity", "Multiplicity", "%.1f", None),
  ("completeness", "Completeness (%)", "%.2f", None),
  ("i_over_sigma", "Mean I/sigma(I)", "%.2f", None),
  ("wilson_b", "Wilson B-factor", "%.2f", None),
  ("r_sym", "R-merge", "%.4g", None),
  ("r_meas", "R-meas", "%.4g", None),
  ("cc_one_half", "CC1/2", "%.3g", None),
  ("cc_star", "CC*", "%.3g", None),
  # model-based statistics
  ("r_work", "R-work", "%.4f", None),
  ("r_free", "R-free", "%.4f", None),
  ("cc_work", "CC(work)", "%.3f", None),
  ("cc_free", "CC(free)", "%.3f", None),
  ("n_atoms", "Number of atoms", "%d", None),
  ("n_macro_atoms", "  macromolecules", "%d", None),
  ("n_ligand_atoms", "  ligands", "%d", None),
  ("n_waters", "  water", "%d", None),
  ("n_residues", "Protein residues", "%d", None),
  ("bond_rmsd", "RMS(bonds)", "%.3f", None),
  ("angle_rmsd", "RMS(angles)", "%.2f", None),
  ("rama_favored", "Ramachandran favored (%)", "%.2g", None),
  ("rama_outliers", "Ramachandran outliers (%)", "%.2g", None),
  ("clashscore", "Clashscore", "%.2f", None),
  ("adp_mean", "Average B-factor", "%.2f", None),
  ("adp_mean_mm", "  macromolecules", "%.2f", None),
  ("adp_mean_lig", "  ligands", "%.2f", None),
  ("adp_mean_wat", "  solvent", "%.2f", None),
  ("solvent_content", "Solvent content", "%.1f%%", None),
]

keyword_formats = dict([ (kw, fs) for (kw, label, fs, cif_tag) in keywords ])

# XXX is this superfluous?
def format_value (fs, value) :
  try :
    val_str = str_utils.format_value(fs, value, replace_none_with="").strip()
  except Exception, e :
    raise RuntimeError("Formatting error: %s, %s" % (fs, value))
  else :
    return val_str

def format_d_max_min (d_max_min) :
  """Format a resolution range (e.g. '30 - 2.56')"""
  if (d_max_min is None) :
    return ""
  else :
    (d_max, d_min) = d_max_min
    d_max_str = "%.4g " % d_max
    d_min_str = re.sub("\.$", ".0", re.sub("0*$", "", "%.3f" % d_min))
    return "%s - %s" % (d_max_str, d_min_str)

class column (slots_getstate_setstate) :
  """
  Statistics for a single structure, including optional high-resolution
  shell.  Any combination of standard attributes is permitted as keyword
  arguments to the constructor.  (Note that the high-resolution shell is
  itself another instance of this class, with fewer keywords specified.)
  """

  __slots__ = [ "outer_shell", "label" ] + [ kw[0] for kw in keywords ]

  def __init__ (self, **kwds) :
    kwds = dict(kwds)
    for name in self.__slots__ :
      if (name in kwds.keys()) :
        setattr(self, name, kwds[name])
      else :
        setattr(self, name, None)
    self.outer_shell = None

  def add_outer_shell (self, **kwds) :
    self.outer_shell = column(**kwds)

  def format_stat (self, name) :
    value = getattr(self, name, None)
    shell_value = getattr(self.outer_shell, name, None)
    fs = keyword_formats[name]
    if (name == "d_max_min") :
      cell_value = format_d_max_min(value)
    else :
      cell_value = format_value(fs, value)
    if (shell_value is not None) :
      if (name == "d_max_min") :
        subvalue = format_d_max_min(shell_value)
      else :
        subvalue = format_value(fs, shell_value)
      cell_value += " (%s)" % subvalue
    return cell_value

  def __repr__ (self) :
    lines = []
    for (stat_name, label, format_string) in keywords :
      value = getattr(self, stat_name, None)
      if (value is not None) :
        lines.append("%s %s" % (stat_name, value))
    return "\n".join(lines)
