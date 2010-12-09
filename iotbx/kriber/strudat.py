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
    if (self._unit_cell is None):
      if (len(self.unit_cell_parameters) == 6):
        self._unit_cell = uctbx.unit_cell(self.unit_cell_parameters)
      else:
        crystal_system = self.space_group().crystal_system()
        unit_cell_parameters = self.unit_cell_parameters
        if (crystal_system == "Cubic"):
          assert len(unit_cell_parameters) == 1
          a = unit_cell_parameters[0]
          self._unit_cell = uctbx.unit_cell((a,a,a,90,90,90))
        elif (crystal_system in ("Hexagonal", "Trigonal")):
          assert len(unit_cell_parameters) == 2
          is_rhombohedral = False
          if (crystal_system == "Trigonal"):
            laue_group = self._derived_laue_group_symbol()
            if (laue_group in ("R-3m:R", "R-3:R")):
              is_rhombohedral = True
          if (is_rhombohedral):
            a = unit_cell_parameters[0]
            angle = unit_cell_parameters[1]
            self._unit_cell = uctbx.unit_cell((a,a,a,angle,angle,angle))
          else:
            a = unit_cell_parameters[0]
            c = unit_cell_parameters[1]
            self._unit_cell = uctbx.unit_cell((a,a,c,90,90,120))
        elif (crystal_system == "Tetragonal"):
          assert len(unit_cell_parameters) == 2
          a = unit_cell_parameters[0]
          c = unit_cell_parameters[1]
          self._unit_cell = uctbx.unit_cell((a,a,c,90,90,90))
        elif (crystal_system == "Orthorhombic"):
          assert len(unit_cell_parameters) == 3
          a = unit_cell_parameters[0]
          b = unit_cell_parameters[1]
          c = unit_cell_parameters[2]
          self._unit_cell = uctbx.unit_cell((a,b,c,90,90,90))
        elif (crystal_system == "Monoclinic"):
          assert len(unit_cell_parameters) == 4
          a = unit_cell_parameters[0]
          b = unit_cell_parameters[1]
          c = unit_cell_parameters[2]
          angle = unit_cell_parameters[3]
          laue_group = self._derived_laue_group_symbol()
          if (laue_group == "P12/m1"):
            self._unit_cell = uctbx.unit_cell((a,b,c,90,angle,90))
          elif (laue_group == "P112/m"):
            self._unit_cell = uctbx.unit_cell((a,b,c,90,90,angle))
          elif (laue_group == "P2/m11"):
            self._unit_cell = uctbx.unit_cell((a,b,c,angle,90,90))
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
