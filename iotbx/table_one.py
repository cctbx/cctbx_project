"""
Tools for formatting the standard "Table 1" in MX papers.  Relies on other
methods (phenix.refine, phenix.model_vs_data, phenix.merging_statistics, etc.)
to extract statistics for display.
"""

from __future__ import absolute_import, division, print_function
from libtbx import slots_getstate_setstate
from libtbx.utils import Sorry
from libtbx import str_utils
import re
from six.moves import range

angstrom = u"\u00C5".encode("utf-8", "strict").strip()

# XXX are these complete, and if not, what is missing?  this appears to cover
# the template suggested by NSMB plus more, but not clear what other journals
# expect.  the VTF paper does not list any explicit Table 1 requirements, but
# the forthcoming PDB validation reports may change expectations.
# TODO fill in missing CIF tags - some of these may not have standard
# tag names yet
keywords = [
  # (attr, label, format,  cif_tag)
  # data-only statistics
  ("wavelength", "Wavelength", "%.4g", None),
  ("d_max_min", "Resolution range", "%.3g - %.3g", None),
  ("space_group", "Space group", "%s", None),
  ("unit_cell", "Unit cell", "%g %g %g %g %g %g", None),
  ("n_refl_all", "Total reflections", "%d", "_reflns.pdbx_number_measured_all"),
  ("n_refl", "Unique reflections", "%d", "_reflns.number_obs"),
  ("multiplicity", "Multiplicity", "%.1f", "_reflns.pdbx_redundancy"),
  ("completeness", "Completeness (%)", "%.2f", "_reflns.percent_possible_obs"),
  ("i_over_sigma", "Mean I/sigma(I)", "%.2f", "_reflns.pdbx_netI_over_sigmaI"),
  ("wilson_b", "Wilson B-factor", "%.2f", "_reflns.B_iso_Wilson_estimate"),
  ("r_sym", "R-merge", "%.4g", "_reflns.pdbx_Rmerge_I_obs"),
  ("r_meas", "R-meas", "%.4g", "_reflns.pdbx_Rrim_I_obs"),
  ("r_pim", "R-pim", "%.4g", "_reflns.pdbx_Rpim_I_obs"),
  ("cc_one_half", "CC1/2", "%.3g", "_reflns.phenix_cc_1/2"),
  ("cc_star", "CC*", "%.3g", "_reflns.phenix_cc_star"),
  # refinement statistics
  # TODO figure out how to extract this...
  ("n_refl_refine", "Reflections used in refinement", "%d", None),
    # XXX this is problematic - I am not sure if our usage is consistent with
    # the intended purpose of this tag
  #  "_refine.ls_number_reflns_obs"),
  ("n_free", "Reflections used for R-free", "%d",
    "_refine.ls_number_reflns_R_free"),
  ("r_work", "R-work", "%.4f", "_refine.ls_R_factor_R_work"),
  ("r_free", "R-free", "%.4f", "_refine.ls_R_factor_R_free"),
  ("cc_work", "CC(work)", "%.3f", None),
  ("cc_free", "CC(free)", "%.3f", None),
  ("n_atoms", "Number of non-hydrogen atoms", "%d", None),
  ("n_macro_atoms", "  macromolecules", "%d", None),
  ("n_ligand_atoms", "  ligands", "%d", None),
  ("n_waters", "  solvent", "%d", None),
  ("n_residues", "Protein residues", "%d", None),
  ("n_nuc", "Nucleic acid bases", "%d", None),
  ("bond_rmsd", "RMS(bonds)", "%.3f", None),
  ("angle_rmsd", "RMS(angles)", "%.2f", None),
  ("rama_favored", "Ramachandran favored (%)", "%.2f", None),
  ("rama_allowed", "Ramachandran allowed (%)", "%.2f", None),
  ("rama_outliers", "Ramachandran outliers (%)", "%.2f", None),
  ("rota_outliers", "Rotamer outliers (%)", "%.2f", None),
  ("clashscore", "Clashscore", "%.2f", None),
  ("adp_mean", "Average B-factor", "%.2f", None),
  ("adp_mean_mm", "  macromolecules", "%.2f", None),
  ("adp_mean_lig", "  ligands", "%.2f", None),
  ("adp_mean_wat", "  solvent", "%.2f", None),
  ("n_tls_groups", "Number of TLS groups", "%d", None),
#  ("solvent_content", "Solvent content", "%.1f%%", None),
]

# statistics that we don't want to show if not applicable (mainly for different
# molecule types)
optional_if_none = set([
  "n_residues", "n_ligand_atoms", "n_waters", "n_nuc",
  "rama_favored", "rama_allowed", "rama_outliers", "rota_outliers",
  "adp_mean_lig", "adp_mean_wat", "cc_work", "cc_free", "cc_star",
  "n_tls_groups",
])

keyword_formats = dict([ (kw, fs) for (kw, label, fs, cif_tag) in keywords ])

class column(slots_getstate_setstate):
  """
  Statistics for a single structure, including optional high-resolution
  shell.  Any combination of standard attributes is permitted as keyword
  arguments to the constructor.  (Note that the high-resolution shell is
  itself another instance of this class, with fewer keywords specified.)
  """

  __slots__ = [ "outer_shell", "label", "anomalous_flag" ] + \
              [ kw[0] for kw in keywords ]

  def __init__(self, **kwds):
    kwds = dict(kwds)
    for name in self.__slots__ :
      if (name in kwds):
        setattr(self, name, kwds[name])
      else :
        setattr(self, name, None)
    self.outer_shell = None

  def add_outer_shell(self, **kwds):
    self.outer_shell = column(**kwds)
    return self

  def format_stat(self, name):
    value = getattr(self, name, None)
    if (value is None) and (name in optional_if_none):
      return None
    shell_value = getattr(self.outer_shell, name, None)
    fs = keyword_formats[name]
    if (name == "d_max_min"):
      cell_value = format_d_max_min(value)
    else :
      cell_value = format_value(fs, value)
    if (shell_value is not None):
      if (name == "d_max_min"):
        subvalue = format_d_max_min(shell_value)
      else :
        subvalue = format_value(fs, shell_value)
      cell_value += " (%s)" % subvalue
    return cell_value

  def __repr__(self):
    lines = []
    for (stat_name, label, format_string, cif_tag) in keywords :
      value = getattr(self, stat_name, None)
      if (value is not None):
        lines.append("%s %s" % (stat_name, value))
    return "\n".join(lines)

class table(slots_getstate_setstate):
  """
  Combined table of statistics for one or more structures published together.
  """
  __slots__ = ["text_field_separation", "columns"]

  def __init__(self, text_field_separation=4):
    self.text_field_separation = text_field_separation
    self.columns = []

  def add_column(self, col):
    self.columns.append(col)

  def format_as_txt(self):
    rows = []
    rows.append([""] + [ column.label for column in self.columns ])
    for (stat_name, label, fstring, cif_tag) in keywords :
      row = []
      for column in self.columns :
        row.append(column.format_stat(stat_name))
      if (row == [ None for x in range(len(row)) ]):
        continue
      row.insert(0, label)
      rows.append(row)
    n_rows = len(rows)
    n_cols = len(self.columns) + 1
    columns = [ [ row[i] for row in rows ] for i in range(n_cols) ]
    columns = [ resize_column(col, "right") for col in columns ]
    table = [ [ col[j] for col in columns ] for j in range(n_rows) ]
    sep = " " * self.text_field_separation
    out = "\n".join([ sep.join(row) for row in table ])
    return out

  def format_as_csv(self):
    rows = []
    rows.append([""] + [ column.label if column.label is not None else "" for column in self.columns ])
    for (stat_name, label, fstring, cif_tag) in keywords :
      row = []
      for column in self.columns :
        stat = column.format_stat(stat_name)
        if stat is None:
          stat = ""
        row.append(stat)
      if ( (row == [ None for x in range(len(row)) ]) or
           (row is None) ):
        continue
      if label is None:
        label = ""
      row.insert(0, label)
      rows.append(row)
    return "\n".join([ ",".join(row) for row in rows ])

  def format_as_rtf(self):
    try :
      import PyRTF
    except ImportError :
      raise Sorry("The PyRTF module is not available.")
    doc = PyRTF.Document()
    ss = doc.StyleSheet
    section = PyRTF.Section()
    doc.Sections.append(section)
    p = PyRTF.Paragraph( ss.ParagraphStyles.Heading2 )
    p.append("Table 1.  Data collection and refinement statistics.")
    section.append(p)
    n_cols = len(self.columns) + 1
    col_widths = [ PyRTF.TabPS.DEFAULT_WIDTH * 3 ] * n_cols
    table = PyRTF.Table(*col_widths)
    header = [ PyRTF.Cell(PyRTF.Paragraph("")) ]
    anomalous_flag = False
    for column in self.columns :
      label = column.label
      if (column.anomalous_flag):
        label += "*"
        anomalous_flag = True
      column_label = PyRTF.Paragraph(ss.ParagraphStyles.Heading2, label)
      header.append(PyRTF.Cell(column_label))
    table.AddRow(*header)
    for (stat_name, label, fstring, cif_tag) in keywords :
      row = [ PyRTF.Cell(PyRTF.Paragraph(ss.ParagraphStyles.Heading2, label)) ]
      n_none = 0
      for column in self.columns :
        txt = column.format_stat(stat_name)
        if (txt is None):
          n_none += 1
        p = PyRTF.Paragraph(ss.ParagraphStyles.Normal,
          PyRTF.ParagraphPS(alignment=2))
        p.append(txt)
        row.append(PyRTF.Cell(p))
      if (n_none == len(row) - 1):
        continue
      table.AddRow(*row)
    section.append(table)
    p = PyRTF.Paragraph(ss.ParagraphStyles.Normal)
    p.append("Statistics for the highest-resolution shell are shown in "+
      "parentheses.")
    section.append(p)
    return doc

  def save_txt(self, file_name):
    f = open(file_name, "w")
    out = self.format_as_txt()
    f.write(out)
    f.close()

  def save_csv(self, file_name):
    f = open(file_name, "w")
    out = self.format_as_csv()
    f.write(out)
    f.close()

  def save_rtf(self, file_name):
    import PyRTF
    DR = PyRTF.Renderer()
    doc = self.format_as_rtf()
    f = open(file_name, "w")
    DR.Write(doc, f)
    f.close()

#-----------------------------------------------------------------------
# UTILITY FUNCTIONS

# XXX is this superfluous?
def format_value(fs, value):
  try :
    val_str = str_utils.format_value(fs, value, replace_none_with="").strip()
  except Exception as e :
    raise RuntimeError("Formatting error: %s, %s" % (fs, value))
  else :
    return val_str

def format_d_max_min(d_max_min):
  """Format a resolution range (e.g. '30 - 2.56')"""
  if (d_max_min is None) or (d_max_min == (None, None)):
    return ""
  else :
    (d_max, d_min) = d_max_min
    d_max_str = "%.4g " % d_max
    d_min_str = re.sub(r"\.$", ".0", re.sub("0*$", "", "%.3f" % d_min))
    return "%s - %s" % (d_max_str, d_min_str)

def resize_column(cell_values, alignment="right"):
  max_width = max([ len(str(cell)) for cell in cell_values ])
  if (alignment == "right"):
    fs = "%%%ds" % max_width
    return [ fs % cell for cell in cell_values ]
  elif (alignment == "left"):
    fs = "%%%-ds" % max_width
    return [ fs % cell for cell in cell_values ]
  else :
    raise RuntimeError("Alignemnt '%s' not supported." % alignment)
