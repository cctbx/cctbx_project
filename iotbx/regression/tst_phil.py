from __future__ import absolute_import, division, print_function
import iotbx.phil
from libtbx.test_utils import show_diff, Exception_expected
from libtbx.utils import Sorry
from libtbx import Auto
import os

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
  # make sure Unicode strings work!
  custom = iotbx.phil.parse(input_string=u"""\
U = 1 2 3 90 105 90
S = P21
""")
  assert not show_diff(
    master.fetch(source=custom).extract_format().as_str(), """\
u = 10 12 13 80 90 100
s = "P 21 21 21"
U = 1 2 3 90 105 90
S = "P 1 21 1"
ua = Auto
sa = Auto
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
sel = altloc ' ' or name C\\*
""")
  user = iotbx.phil.parse(input_string="""\
sel = "altloc ' ' or name C\\*"
sel = altloc ' ' or name D\\*
""")
  work = master.fetch(user)
  assert not show_diff(work.as_str(), """\
sel = "altloc ' ' or name C\\\\*"
sel = altloc ' ' or name D\\*
""")
  ex = work.extract()
  assert len(ex.sel) == 2
  assert ex.sel[0] == "altloc ' ' or name C\\*"
  assert ex.sel[1] == "altloc ' ' or name D\\*"
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
  with open("model.pdb", "w") as f:
    f.write("""\
ATOM      1  O   HOH     1      53.448  18.599 -10.134  1.00 20.00
""")
  pcl = iotbx.phil.process_command_line_with_files(
    args=["model.pdb", "--use_geometry_restraints"],
    master_phil_string=master_phil_str,
    pdb_file_def="model")
  assert (pcl.get_file_type_count("pdb") == 1)
  assert (pcl.get_file_type_count("hkl") == 0)
  params = pcl.work.extract()
  pdb_in = pcl.get_cached_file(params.model)
  assert (pdb_in is not None)
  assert (params.model == os.path.join(os.getcwd(), "model.pdb"))
  assert (params.use_geometry_restraints)
  model_in = pcl.get_file(params.model, force_type="pdb")
  master_phil_str = """
data_dir = None
  .type = path
d_min = None
  .type = float
nproc = None
  .type = int
"""
  pcl = iotbx.phil.process_command_line_with_files(
    args=["/", "3.0", "5"],
    master_phil_string=master_phil_str,
    directory_def="data_dir",
    integer_def="nproc",
    float_def="d_min")
  params = pcl.work.extract()
  assert (params.data_dir == "/")
  assert (params.d_min == 3.0)
  assert (params.nproc == 5)
  # shelx file hack
  with open("tst_iotbx_phil.hkl", "w") as f:
    f.write("""\
   1   2  -1  -23.34    4.56   1
   2  -3   9   12.45    6.12   2
99999999999999999.9999999.999999
-999-999-999-9999.99-9999.99-999
   0   0   0    0.00    0.00   0""")
  master_phil_str = """
data = None
  .type = path
space_group = None
  .type = space_group
unit_cell = None
  .type = unit_cell
"""
  pcl = iotbx.phil.process_command_line_with_files(
    args=[
      "tst_iotbx_phil.hkl",
      "P6122",
      "50,50,40,90,90,120",
    ],
    master_phil_string=master_phil_str,
    reflection_file_def="data",
    space_group_def="space_group",
    unit_cell_def="unit_cell")
  params = pcl.work.extract()
  assert (str(params.space_group) == "P 61 2 2")
  assert (str(params.unit_cell) == "(50, 50, 40, 90, 90, 120)")
  from iotbx.file_reader import any_file
  try :
    hkl_in = any_file(params.data, force_type="hkl")
    print(hkl_in.file_server.miller_arrays[0].is_xray_intensity_array())
  except Sorry as s :
    assert ("Unresolved amplitude/intensity ambiguity" in str(s))
  else :
    raise Exception_expected
  pcl = iotbx.phil.process_command_line_with_files(
    args=["tst_iotbx_phil.hkl=hklf3"],
    master_phil_string=master_phil_str,
    reflection_file_def="data")
  params = pcl.work.extract()
  assert (params.data == "tst_iotbx_phil.hkl=hklf3")
  hkl_in = any_file(params.data, force_type="hkl")
  ma = hkl_in.file_server.miller_arrays[0]
  assert ma.is_xray_amplitude_array()
  assert ma.anomalous_flag() is False

  with open("tst_iotbx_phil.ncs", "w") as f:
    f.write("""
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2  0.623490  0.781831  0.000000        0.00000
REMARK 350   BIOMT2   2 -0.781831  0.623490  0.000000        0.00000
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000

""")
  master_phil_str = """
ncs_file = None
  .type = path
"""
  pcl = iotbx.phil.process_command_line_with_files(
    args=[
      "tst_iotbx_phil.ncs",
    ],
    master_phil_string=master_phil_str,
    ncs_file_def="ncs_file")
  params = pcl.work.extract()
  assert (os.path.split(str(params.ncs_file))[-1] == "tst_iotbx_phil.ncs")

  with open("tst_iotbx_phil.ncs_spec", "w") as f:
    f.write("""
Summary of NCS information
Wed Mar 16 15:02:22 2016
/Users/terwill/Desktop/working/Jan_2015/phenix/cryo-pdb/3j9c/build

source_info ncs_biomtr.ncs




new_ncs_group
new_operator

rota_matrix    1.0000    0.0000    0.0000
rota_matrix    0.0000    1.0000    0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    0.0000    0.0000    0.0000
new_operator

rota_matrix    0.6235   -0.7818    0.0000
rota_matrix    0.7818    0.6235   -0.0000
rota_matrix   -0.0000    0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

""")
  master_phil_str = """
ncs_file = None
  .type = path
"""
  pcl = iotbx.phil.process_command_line_with_files(
    args=[
      "tst_iotbx_phil.ncs_spec",
    ],
    master_phil_string=master_phil_str,
    ncs_file_def="ncs_file")
  params = pcl.work.extract()
  assert (os.path.split(str(params.ncs_file))[-1] == "tst_iotbx_phil.ncs_spec")



  print("OK")

if (__name__ == "__main__"):
  exercise()
