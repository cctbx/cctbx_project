from cctbx import uctbx
from cctbx import sgtbx
from scitbx.python_utils import easy_pickle
from scitbx.python_utils.misc import adopt_init_args
import os

class sdb_site:

  def __init__(self, action, segid, type, x, y, z, b, q, g):
    adopt_init_args(self, locals())

  def as_xray_scatterer(self, unit_cell=None):
    from cctbx import adptbx
    from cctbx import xray
    from cctbx.eltbx.caasf import wk1995
    caasf = None
    try: caasf = wk1995(self.type)
    except:
      try: caasf = wk1995(self.segid)
      except: pass
    if (caasf is None): caasf = wk1995("const")
    site = (self.x, self.y, self.z)
    if (unit_cell is not None): site = unit_cell.fractionalize(site)
    return xray.scatterer(
      label="_".join((self.segid, self.type)),
      site=site,
      u=adptbx.b_as_u(self.b),
      occupancy=self.q,
      caasf=caasf)

class sdb_file:

  def __init__(self, file_name, unit_cell, space_group_info, sites):
    adopt_init_args(self, locals())

  def as_xray_structure(self, min_distance_sym_equiv=0.5):
    assert self.unit_cell is not None
    assert self.space_group_info is not None
    from cctbx import crystal
    from cctbx import xray
    symmetry = crystal.symmetry(
      unit_cell=self.unit_cell,
      space_group_info=self.space_group_info)
    structure = xray.structure(crystal.special_position_settings(
      crystal_symmetry=symmetry,
      min_distance_sym_equiv=min_distance_sym_equiv))
    for site in self.sites:
      structure.add_scatterer(site.as_xray_scatterer(self.unit_cell))
    return structure

def generic_add_str(m, buffer):
  if (not m): return
  seq_no = int(m.group(1))
  assert seq_no == len(buffer) + 1
  buffer.append(m.group(2))

def generic_add_float(m, buffer):
  if (not m): return
  seq_no = int(m.group(1))
  assert seq_no == len(buffer) + 1
  buffer.append(float(m.group(2)))

class raw_parameters:

  def __init__(self, file_name):
    self.file_name = file_name
    self.unit_cell = None
    self.space_group_info = None
    self.action = []
    self.segid = []
    self.type = []
    self.x = []
    self.y = []
    self.z = []
    self.b = []
    self.q = []
    self.g = []

  def add_action(self, m): generic_add_str(m, self.action)
  def add_segid(self, m): generic_add_str(m, self.segid)
  def add_type(self, m): generic_add_str(m, self.type)
  def add_x(self, m): generic_add_float(m, self.x)
  def add_y(self, m): generic_add_float(m, self.y)
  def add_z(self, m): generic_add_float(m, self.z)
  def add_b(self, m): generic_add_float(m, self.b)
  def add_q(self, m): generic_add_float(m, self.q)
  def add_g(self, m): generic_add_str(m, self.g)

  def as_sdb_sites(self):
    assert len(self.segid) == len(self.action)
    assert len(self.type) == len(self.action)
    assert len(self.x) == len(self.action)
    assert len(self.y) == len(self.action)
    assert len(self.z) == len(self.action)
    assert len(self.b) == len(self.action)
    assert len(self.q) == len(self.action)
    assert len(self.g) == len(self.action)
    sites = []
    for i in xrange(len(self.x)):
      sites.append(sdb_site(
        self.action[i], self.segid[i], self.type[i],
        self.x[i], self.y[i], self.z[i],
        self.b[i], self.q[i],
        self.g[i]))
    return sdb_file(
      self.file_name, self.unit_cell, self.space_group_info, sites)

def multi_sdb_parser(lines):
  # Parser for one or more cns sdb files.
  # Lines interpreted:
  #   {+ file: heavy_search_1.sdb +}
  #   sg= P6 a= 116.097 b= 116.097 c= 44.175 alpha= 90 beta= 90 gamma= 120
  #   {===>} site.action_1="refine";
  #   {===>} site.segid_1="SITE"; site.type_1="SE";
  #   {===>} site.x_1=18.7869; site.y_1=12.1257; site.z_1=0.163635;
  #   {===>} site.b_1=65.6408; site.q_1=1; site.g_1="";
  # Sites must be sorted.
  import re
  sdb_files = []
  p = 0
  for line in lines:
    m = re.search(r'\{\+\s+file:\s*(\S*)', line)
    if (m):
      if (p): sdb_files.append(p.as_sdb_sites())
      p = raw_parameters(m.group(1))
    if (not p): continue
    m = re.match(  r'sg=\s*(\S+)\s*a=\s*(\S+)\s*b=\s*(\S+)\s*c=\s*(\S+)'
                 + r'\s*alpha=\s*(\S+)\s*beta=\s*(\S+)\s*gamma=\s*(\S+)',
                 line)
    if (m):
      p.unit_cell = uctbx.unit_cell(
        [float(m.group(i+2)) for i in xrange(6)])
      p.space_group_info = sgtbx.space_group_info(m.group(1))
    else:
      m = re.match(r'\{===>\}\s*sg=\s*"(\S+)"\s*;', line)
      if (m):
        p.space_group_info = sgtbx.space_group_info(m.group(1))
    p.add_action(re.search(r'site\.action_(\d+)\s*=\s*"([^"]*)"', line))
    p.add_segid(re.search(r'site\.segid_(\d+)\s*=\s*"([^"]*)"', line))
    p.add_type(re.search(r'site\.type_(\d+)\s*=\s*"([^"]*)"', line))
    p.add_x(re.search(r'site\.x_(\d+)\s*=\s*([^\s;]*)', line))
    p.add_y(re.search(r'site\.y_(\d+)\s*=\s*([^\s;]*)', line))
    p.add_z(re.search(r'site\.z_(\d+)\s*=\s*([^\s;]*)', line))
    p.add_b(re.search(r'site\.b_(\d+)\s*=\s*([^\s;]*)', line))
    p.add_q(re.search(r'site\.q_(\d+)\s*=\s*([^\s;]*)', line))
    p.add_g(re.search(r'site\.g_(\d+)\s*=\s*"([^"]*)"', line))
  if (p): sdb_files.append(p.as_sdb_sites())
  return sdb_files

def run(args):
  show_raw = "--show_raw" in args
  write_pickle = "--pickle" in args
  unit_cell = None
  for arg in args:
    if (arg.startswith("--unit_cell=")):
      params = arg.split("=", 1)[1]
      unit_cell = uctbx.unit_cell(params)
  for file_name in args:
    if (file_name.startswith("--")): continue
    f = open(file_name, "r")
    lines = f.readlines()
    f.close()
    sdb_files = multi_sdb_parser(lines)
    for sdb in sdb_files:
      if (unit_cell is not None): sdb.unit_cell = unit_cell
      print "file:", sdb.file_name
      if (sdb.unit_cell is not None):
        print "unit cell:", sdb.unit_cell.parameters()
      if (sdb.space_group_info is not None):
        print "space group:", sdb.space_group_info
      if (show_raw):
        for site in sdb.sites:
          print site.action, site.segid, site.type, site.g
          print " ", site.x, site.y, site.z, site.b, site.q
      else:
        xray_structure = sdb.as_xray_structure()
        xray_structure.show_summary().show_scatterers()
        if (write_pickle):
          file_name_pickle = os.path.split(sdb.file_name)[1] + ".pickle"
          print "Writing:", file_name_pickle
          easy_pickle.dump(file_name_pickle, xray_structure)
      print
