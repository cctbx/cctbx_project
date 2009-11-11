import sys, os
op = os.path

nodes_32bit = set("""\
ribbon
longnose
sharptail
""".splitlines())

nodes_64bit = set("""\
krait
anaconda
rosie
chevy
firtree
""".splitlines())

def node_os_bits(node):
  if (node in nodes_32bit): return "32b"
  if (node in nodes_64bit): return "64b"
  raise RuntimeError

class table_entry(object):

  __slots__ = [
    "current_platform",
    "current_node",
    "build_platform",
    "build_node",
    "gcc_static",
    "compiler_versions",
    "n_scatt",
    "n_refl",
    "all_utimes"]

  def __init__(O, lines):
    assert len(lines) == 16
    for i_line in xrange(5):
      line = lines[i_line]
      slot = O.__slots__[i_line]
      assert line.startswith(slot+": ")
      setattr(O, slot, line[len(slot)+2:])
    O.compiler_versions = []
    for i_line in xrange(5, 10):
      line = lines[i_line]
      assert line.startswith("compiler: ")
      O.compiler_versions.append(line[10:])
    O.current_node = O.current_node.replace(".lbl.gov", "")
    O.build_node = O.build_node.replace(".lbl.gov", "")
    line = lines[10]
    assert line.startswith("n_scatt * n_refl: ")
    flds = line.split()
    assert len(flds) == 6
    O.n_scatt = int(flds[3])
    O.n_refl = int(flds[5])
    O.all_utimes = []
    for i_line in xrange(11, 16):
      line = lines[i_line]
      flds = line.split()
      assert len(flds) == 8
      O.all_utimes.append([float(fld) for fld in flds])

  def format_utimes(O, i_compiler):
    return " ".join(["%6.2f" % u for u in O.all_utimes[i_compiler]])

  def ifort_label(O):
    ic = 0
    if (O.all_utimes[0][0] < 0):
      return None, None
    v = O.compiler_versions[ic].split()[2]
    assert v in ["9.1", "11.1"]
    if (v == "9.1"): v = "_"+v
    return ic, "if%s_%s" % (v, node_os_bits(O.build_node))

  def gfort_label(O):
    lbls = {
      "GNU Fortran (GCC 3.2 20020903 (Red Hat Linux 8.0 3.2-7))"
        " 3.2 20020903 (Red Hat Linux 8.0 3.2-7)": "gf32",
      "GNU Fortran (GCC) 3.4.2 20041017 (Red Hat 3.4.2-6.fc3)": "gf34",
      "GNU Fortran 95 (GCC 4.0.0 20041019 (Red Hat 4.0.0-0.8))": "gf40",
      "GNU Fortran 95 (GCC 4.0.0 20050519 (Red Hat 4.0.0-8))": "gf40",
      "GNU Fortran (GCC) 4.1.2 20070925 (Red Hat 4.1.2-33)": "gf41",
      "GNU Fortran (GCC) 4.2.2": "gf42",
      "GNU Fortran (GCC) 4.3.4": "gf43",
      "GNU Fortran (GCC) 4.4.2": "gf44"}
    if (O.all_utimes[1][0] >= 0):
      ic = 1
    elif (O.all_utimes[2][0] >= 0):
      ic = 2
    else:
      return None, None
    v = lbls[O.compiler_versions[ic]]
    return ic, "%s_%s" % (v, node_os_bits(O.build_node))

  def gpp_label(O):
    lbls = {
      "g++ (GCC) 3.2 20020903 (Red Hat Linux 8.0 3.2-7)": "gc32",
      "g++ (GCC) 3.3.4 (pre 3.3.5 20040809)": "gc33",
      "g++ (GCC) 3.4.2 20041017 (Red Hat 3.4.2-6.fc3)": "gc34",
      "g++ (GCC) 4.0.0 20050519 (Red Hat 4.0.0-8)": "gc40",
      "g++ (GCC) 4.1.2 20070925 (Red Hat 4.1.2-33)": "gc41",
      "g++ (GCC) 4.2.2": "gc42",
      "g++ (GCC) 4.3.4": "gc43",
      "g++ (GCC) 4.4.2": "gc44"}
    ic = 4
    v = lbls[O.compiler_versions[ic]]
    return ic, "%s_%s" % (v, node_os_bits(O.build_node))

def process_time_tables():
  import libtbx.load_env
  time_tables_path = libtbx.env.find_in_repositories(
    relative_path="compcomm/newsletter09/time_tables",
    test=op.isfile,
    optional=False)
  result = []
  lines = open(time_tables_path).read().splitlines()
  assert len(lines) == 22 * 17
  for i_line in xrange(0, len(lines), 17):
    assert lines[i_line] == ""
    result.append(table_entry(lines=lines[i_line+1:i_line+17]))
  return result

def run(args):
  assert len(args) == 0
  table_entries = process_time_tables()
  #
  print "Intel Fortran times on chevy:"
  tab = {}
  for entry in table_entries:
    if (entry.current_node == "chevy"):
      ic, lbl = entry.ifort_label()
      if (ic is None): continue
      if (    lbl == "if_9.1_64b"
          and entry.build_node != "anaconda"):
        continue
      if (    lbl == "if_9.1_32b"
          and entry.build_node != "longnose"):
        continue
      tab[lbl] = (entry.format_utimes(i_compiler=ic), entry.build_node)
  for key in ["if_9.1_32b", "if_9.1_64b", "if11.1_64b"]:
    print key, " ".join(tab[key])
  print
  #
  print "GNU Fortran times on chevy:"
  tab = {}
  for entry in table_entries:
    if (entry.current_node == "chevy"):
      ic, lbl = entry.gfort_label()
      if (ic is None): continue
      tab[lbl] = (entry.format_utimes(i_compiler=ic), entry.build_node)
  for key in sorted(tab.keys()):
    print key, " ".join(tab[key])
  print
  #
  print "Intel C++ times on chevy:"
  for entry in table_entries:
    if (entry.current_node == "chevy" and entry.build_node == "chevy"):
      if (entry.all_utimes[3][0] >= 0):
        print entry.format_utimes(i_compiler=3), \
          entry.compiler_versions[3].replace(" (ICC)", "")
  print
  #
  print "GNU C++ times on chevy:"
  tab = {}
  for entry in table_entries:
    if (entry.current_node == "chevy"):
      ic, lbl = entry.gpp_label()
      if (ic is None): continue
      tab[lbl] = (entry.format_utimes(i_compiler=ic), entry.build_node)
  for key in sorted(tab.keys()):
    print key, " ".join(tab[key])
  print
  #
  print "GNU C++ 3.2 32-bit (Red Hat 8.0) executable on all platforms:"
  for entry in table_entries:
    if (entry.build_node == "ribbon"):
      print entry.format_utimes(i_compiler=4), entry.current_node
  print

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
