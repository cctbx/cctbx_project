from cctbx.misc import python_utils

class sdb_site:
  def __init__(self, action, segid, type, x, y, z, b, q, g):
    python_utils.adopt_init_args(self, locals())

class sdb_file:
  def __init__(self, file_name, sites):
    python_utils.adopt_init_args(self, locals())

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
    return sdb_file(self.file_name, sites)

def multi_sdb_parser(lines):
  # Parser for one or more cns sdb files.
  # Lines interpreted:
  #   {+ file: heavy_search_1.sdb +}
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

if (__name__ == "__main__"):
  import sys
  for file_name in sys.argv[1:]:
    f = open(file_name, "r")
    lines = f.readlines()
    f.close()
    sdb_files = multi_sdb_parser(lines)
    for sdb in sdb_files:
      print "file:", sdb.file_name
      for site in sdb.sites:
        print site.action, site.segid, site.type, site.g
        print " ", site.x, site.y, site.z, site.b, site.q
