import iotbx.phil
from libtbx.test_utils import show_diff
from libtbx import Auto
import os
import sys

def exercise():
  master = iotbx.phil.parse(input_string="""\
u=10,12 13 80,90 100
  .type=unit_cell
s=19
  .type=space_group
U=None
  .type=unit_cell
S=None
  .type=space_group
ua=Auto
  .type=unit_cell
sa=Auto
  .type=space_group
""")
  assert not show_diff(master.format(master.extract()).as_str(), """\
u = 10 12 13 80 90 100
s = "P 21 21 21"
U = None
S = None
ua = Auto
sa = Auto
""")
  custom = iotbx.phil.parse(input_string="""\
ua = None
s = Auto
U = 1,2,3
u = Auto
S = 33
sa = None
""")
  assert not show_diff(
    master.fetch(source=custom).extract_format().as_str(), """\
u = Auto
s = Auto
U = 1 2 3 90 90 90
S = "P n a 21"
ua = None
sa = None
""")
  #
  master = iotbx.phil.parse(input_string="""\
s=None
  .type=space_group
  .multiple=True
""")
  custom = iotbx.phil.parse(input_string="""\
s=p212121
s=19
""")
  assert not show_diff(master.fetch(source=custom).as_str(), """\
s = 19
""")
  custom = iotbx.phil.parse(input_string="""\
s=19
s=p212121
""")
  assert not show_diff(master.fetch(source=custom).as_str(), """\
s = p212121
""")
  custom = iotbx.phil.parse(input_string="""\
s=18
s=p212121
""")
  assert not show_diff(master.fetch(source=custom).as_str(), """\
s = 18
s = p212121
""")
  #
  master = iotbx.phil.parse(input_string="""\
sel = None
  .type = atom_selection
  .multiple = True
""")
  clai = master.command_line_argument_interpreter()
  assert not show_diff(clai.process(
    "sel=altloc  ' ' or  name C\\* ").as_str(), """\
sel = altloc ' ' or name C\*
""")
  user = iotbx.phil.parse(input_string="""\
sel = "altloc ' ' or name C\*"
sel = altloc ' ' or name D\*
""")
  work = master.fetch(user)
  assert not show_diff(work.as_str(), """\
sel = "altloc ' ' or name C\\\\*"
sel = altloc ' ' or name D\\*
""")
  ex = work.extract()
  assert len(ex.sel) == 2
  assert ex.sel[0] == "altloc ' ' or name C\*"
  assert ex.sel[1] == "altloc ' ' or name D\*"
  #
  pcl = iotbx.phil.process_command_line(args=["s = 230"], master_string="""\
s=None
  .type=space_group
""")
  assert not show_diff(str(pcl.work.extract().s), "I a -3 d")
  #
  master_phil_str = """
model = None
  .type = path
use_geometry_restraints = False
  .type = bool
"""
  open("model.pdb", "w").write("""\
ATOM      1  O   HOH     1      53.448  18.599 -10.134  1.00 20.00
""")
  pcl = iotbx.phil.process_command_line_with_files(
    args=["model.pdb", "--use_geometry_restraints"],
    master_phil_string=master_phil_str,
    pdb_file_def="model")
  params = pcl.work.extract()
  assert (params.model == os.path.join(os.getcwd(), "model.pdb"))
  assert (params.use_geometry_restraints)
  master_phil_str = """
data_dir = None
  .type = path
d_min = None
  .type = float
nproc = None
  .type = int
"""
  pcl = iotbx.phil.process_command_line_with_files(
    args=["/var/tmp", "3.0", "5"],
    master_phil_string=master_phil_str,
    directory_def="data_dir",
    integer_def="nproc",
    float_def="d_min")
  params = pcl.work.extract()
  if (sys.platform in ["linux2", "darwin"]) :
    assert (params.data_dir == "/var/tmp")
  assert (params.d_min == 3.0)
  assert (params.nproc == 5)
  print "OK"

if (__name__ == "__main__"):
  exercise()
