from cctbx import xray
from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
import sys

class atom(object):

  def __init__(self, label, site, connectivity):
    self.label = label
    self.site = tuple(site)
    self.connectivity = connectivity

  def as_xray_scatterer(self):
    return xray.scatterer(label=self.label, site=self.site)

class read_entry(object):

  def __init__(self, f):
    line = f.readline()
    if (len(line) == 0):
      self.tag = False
      return
    if (line[0] != "*"):
      self.tag = None
      return
    self.tag = line[1:].strip()
    assert len(self.tag) != 0
    self.title = f.readline().strip()
    self.reference = f.readline().strip()
    self.space_group_symbol = f.readline().strip()
    assert len(self.space_group_symbol) != 0
    self.unit_cell_parameters = [float(p) for p in f.readline().split()]
    assert len(self.unit_cell_parameters) != 0
    self.atoms = []
    while 1:
      line = f.readline()
      if (len(line) == 0): break
      if (line.startswith("--")): break
      atom_line_fields = line.split("#", 1)[0].split()
      if (len(atom_line_fields) == 0): continue
      assert len(atom_line_fields) in (4,5), line.rstrip()
      site = [float(x) for x in atom_line_fields[1:4]]
      if (len(atom_line_fields) > 4):
        connectivity = int(atom_line_fields[4])
      else:
        connectivity = None
      self.atoms.append(atom(atom_line_fields[0], site, connectivity))
    self._space_group_info = None
    self._unit_cell = None

  def space_group_info(self):
    if (self._space_group_info is None):
      self._space_group_info = sgtbx.space_group_info(
        self.space_group_symbol, table_id="A1983")
    return self._space_group_info

  def space_group(self):
    return self.space_group_info().group()

  def _derived_laue_group_symbol(self):
    return str(sgtbx.space_group_info(
      group=self.space_group().build_derived_laue_group())).replace(" ", "")

  def unit_cell(self):
    self._unit_cell = uctbx.infer_unit_cell_from_symmetry(
      self.unit_cell_parameters, self.space_group())
    if (self._unit_cell is None):
      raise RuntimeError, "Cannot interpret unit cell parameters."
    return self._unit_cell

  def crystal_symmetry(self):
    return crystal.symmetry(
      unit_cell=self.unit_cell(),
      space_group_info=self.space_group_info())

  def show(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "Tag:", self.tag
    print >> f, "Title:", self.title
    print >> f, "Reference:", self.reference
    self.crystal_symmetry().show_summary(f=f)
    for atm in self.atoms:
      print >> f, "%-8s" % atm.label,
      print "%8.5f %8.5f %8.5f" % atm.site,
      if (atm.connectivity != None):
        print "%2d" % atm.connectivity,
      print

  def connectivities(self, all_or_nothing=False):
    result = [atom.connectivity for atom in self.atoms]
    if (all_or_nothing):
      n_none = result.count(None)
      if (n_none == len(result)):
        result = None
      elif (n_none != 0):
        raise AssertionError("Tag %s: %d atom%s missing the bond count."
          % (self.tag, n_none, [" is","s are"][int(n_none!=1)]))
    return result

  def as_xray_structure(self, min_distance_sym_equiv=0.5):
    result = xray.structure(
      special_position_settings=crystal.special_position_settings(
        crystal_symmetry=self.crystal_symmetry(),
        min_distance_sym_equiv=min_distance_sym_equiv))
    for atm in self.atoms:
      result.add_scatterer(atm.as_xray_scatterer())
    return result

class read_all_entries(object):

  def __init__(self, f):
    self.entries = []
    while 1:
      entry = read_entry(f)
      if (entry.tag is None):
        continue
      if (entry.tag == False):
        break
      self.entries.append(entry)

  def show(self, f=None):
    if (f is None): f = sys.stdout
    for entry in self.entries:
      entry.show(f=f)
      print >> f

  def __call__(self):
    return self.entries

  def get(self, tag):
    for entry in self.entries:
      if (entry.tag == tag):
        return entry
    return None
