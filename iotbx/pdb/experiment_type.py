"""Class for holding information about experiment method (type)."""
from __future__ import absolute_import, division, print_function

class experiment_type(object):
  """ Class for holding information about experiment method (type).

  It is recorded in EXPDTA field of PDB format and _exptl.method
  of mmCIF format. Additional information:
  https://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#EXPDTA
  https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_exptl.method.html
  """

  def __init__(self, lines):
    """
      Initialization

    Args:
        lines (list of lines): list of lines with methods extracted
          from appropriate places
    """
    assert isinstance(lines, list)
    self.lines = []
    for l in lines:
      self.lines.append(l.strip().upper())

  def __repr__(self):
    return "; ".join(self.lines)

  def is_xray(self):
    return "X-RAY DIFFRACTION" in self.lines

  def is_electron_microscopy(self):
    return "ELECTRON MICROSCOPY" in self.lines

  def is_neutron(self):
    return "NEUTRON DIFFRACTION" in self.lines

  def is_join_xn(self):
    return self.is_xray() and self.is_neutron()

  def is_empty(self):
    return len(self.lines) == 0
