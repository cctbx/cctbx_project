
from __future__ import division
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from io import StringIO

def exercise () :
  from mmtbx.command_line import b_factor_statistics
  from mmtbx.regression import make_fake_anomalous_data
  pdb_file = make_fake_anomalous_data.write_pdb_input_cd_cl(
    file_base="tst_b_factor_stats")
  out = StringIO()
  stats = b_factor_statistics.run(args=[pdb_file], out=out)
  assert ("""\
| all     : 79     0      7.73    42.64   13.18    None  None   None   |""" in
  out.getvalue())
  out = StringIO()
  stats = b_factor_statistics.run(args=[pdb_file, "selection=element CD"],
    out=out)
  assert ("""\
| all(noH): 1      0      15.00   15.00   15.00    None  None   None   |""" in
    out.getvalue())
  print("OK")

if (__name__ == "__main__") :
  exercise()
