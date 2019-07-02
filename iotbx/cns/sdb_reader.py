from __future__ import absolute_import, division, print_function
from iotbx.cns.crystal_symmetry_utils import \
  re_sg_uc, crystal_symmetry_from_re_match
from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import adptbx
from cctbx import xray
import cctbx.eltbx.xray_scattering
from cctbx import eltbx
from libtbx import adopt_init_args
from libtbx import easy_pickle
import re
import os
from six.moves import range

class sdb_site(object):

  def __init__(self, action, segid, type, x, y, z, b, q, g):
    adopt_init_args(self, locals())

  def as_xray_scatterer(self, unit_cell=None):
    scattering_type = eltbx.xray_scattering.get_standard_label(
      label=self.type, exact=False, optional=True)
    if (scattering_type is None):
      scattering_type = eltbx.xray_scattering.get_standard_label(
        label=self.segid, exact=False, optional=True)
    if (scattering_type is None): scattering_type = "unknown"
    site = (self.x, self.y, self.z)
    if (unit_cell is not None): site = unit_cell.fractionalize(site)
    return xray.scatterer(
      label="_".join((self.segid, self.type)),
      site=site,
      u=adptbx.b_as_u(self.b),
      occupancy=self.q,
      scattering_type=scattering_type)

class sdb_file(object):

  def __init__(self, file_name, unit_cell, space_group_info, sites):
    adopt_init_args(self, locals())

  def crystal_symmetry(self):
    return crystal.symmetry(
      unit_cell=self.unit_cell,
      space_group_info=self.space_group_info)

  def as_xray_structure(self, crystal_symmetry=None, force_symmetry=False,
                              min_distance_sym_equiv=0.5):
    crystal_symmetry = self.crystal_symmetry().join_symmetry(
      other_symmetry=crystal_symmetry,
      force=force_symmetry)
    assert crystal_symmetry.unit_cell() is not None
    assert crystal_symmetry.space_group_info() is not None
    structure = xray.structure(crystal.special_position_settings(
      crystal_symmetry=crystal_symmetry,
      min_distance_sym_equiv=min_distance_sym_equiv))
    for site in self.sites:
      structure.add_scatterer(site.as_xray_scatterer(structure.unit_cell()))
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

class raw_parameters(object):

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
    for i in range(len(self.x)):
      sites.append(sdb_site(
        self.action[i], self.segid[i], self.type[i],
        self.x[i], self.y[i], self.z[i],
        self.b[i], self.q[i],
        self.g[i]))
    return sdb_file(
      self.file_name, self.unit_cell, self.space_group_info, sites)

def multi_sdb_parser(lines, file_name=None, max_characters=1000000):
  # Parser for one or more cns sdb files.
  # Lines interpreted:
  #   {+ file: heavy_search_1.sdb +}
  #   sg= P6 a= 116.097 b= 116.097 c= 44.175 alpha= 90 beta= 90 gamma= 120
  #   {===>} site.action_1="refine";
  #   {===>} site.segid_1="SITE"; site.type_1="SE";
  #   {===>} site.x_1=18.7869; site.y_1=12.1257; site.z_1=0.163635;
  #   {===>} site.b_1=65.6408; site.q_1=1; site.g_1="";
  # Sites must be sorted.
  sdb_files = []
  block_name = None
  current_symmetry = None
  n_characters = 0
  p = 0
  for line in lines:
    if (max_characters != 0):
      n_characters += len(line)
      if (n_characters > max_characters): break
    m = re.search(r'\{\+\s+file:\s*(\S*)', line)
    if (m):
      block_name = m.group(1)
    m = re.match(re_sg_uc, line)
    if (m):
      current_symmetry = crystal_symmetry_from_re_match(m=m)
    m = re.search(
      r'\{\-\s+begin\s+block\s+parameter\s+definition\s+\-\}', line)
    if (m):
      if (block_name is None):
        i = len(sdb_files) + 1
        if (p): i += 1
        if (file_name is None):
          block_name = "block_%d" % i
        else:
          block_name = file_name + "_%d" % i
      if (p): sdb_files.append(p.as_sdb_sites())
      p = raw_parameters(block_name)
      if (current_symmetry is None):
        p.unit_cell = None
        p.space_group_info = None
      else:
        p.unit_cell = current_symmetry.unit_cell()
        p.space_group_info = current_symmetry.space_group_info()
        current_symmetry = None
      block_name = None
    if (not p): continue
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
  space_group_info = None
  for arg in args:
    if (arg.startswith("--unit_cell=")):
      params = arg.split("=", 1)[1]
      unit_cell = uctbx.unit_cell(params)
    elif (arg.startswith("--space_group=")):
      symbol = arg.split("=", 1)[1]
      space_group_info = sgtbx.space_group_info(symbol)
  for file_name in args:
    if (file_name.startswith("--")): continue
    f = open(file_name, "r")
    lines = f.readlines()
    f.close()
    sdb_files = multi_sdb_parser(lines, file_name)
    for sdb in sdb_files:
      if (unit_cell is not None): sdb.unit_cell = unit_cell
      if (space_group_info is not None): sdb.space_group_info=space_group_info
      print("file:", sdb.file_name)
      if (sdb.unit_cell is not None):
        print("unit cell:", sdb.unit_cell.parameters())
      if (sdb.space_group_info is not None):
        print("space group:", sdb.space_group_info)
      if (show_raw):
        for site in sdb.sites:
          print(site.action, site.segid, site.type, site.g)
          print(" ", site.x, site.y, site.z, site.b, site.q)
      else:
        xray_structure = sdb.as_xray_structure()
        xray_structure.show_summary().show_scatterers()
        if (write_pickle):
          file_name_pickle = os.path.split(sdb.file_name)[1] + ".pickle"
          print("Writing:", file_name_pickle)
          easy_pickle.dump(file_name_pickle, xray_structure)
      print()
