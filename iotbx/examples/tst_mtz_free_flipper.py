import libtbx.load_env
if (libtbx.env.has_module("ccp4io")):
  from iotbx.examples import mtz_free_flipper
  from iotbx.examples import mtz_convert_free_to_work
  import iotbx.mtz as iotbx_mtz
else:
  iotbx_mtz = None
from libtbx.test_utils import show_diff
from libtbx.utils import format_cpu_times
from cStringIO import StringIO
import os

def exercise():
  if (iotbx_mtz is None):
    print \
      "Skipping iotbx/examples/tst_mtz_free_flipper.py: ccp4io not available"
    return
  input_file_name = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/l.mtz",
    test=os.path.isfile)
  if (input_file_name is None):
    print "Skipping exercise(): input file not available"
    return
  label = "R-free-flags(-)"
  mtz_free_flipper.run(args=[input_file_name], label=label)
  mtz_free_flipper.run(args=["free_flipped_l.mtz"], label=label)
  mtz_convert_free_to_work.run(args=[input_file_name], label=label)
  spreadsheets = []
  for file_name, expected in [
        (input_file_name, (13469,1065,2323)),
        ("free_flipped_l.mtz", (1065,13469,2323)),
        ("free_flipped_free_flipped_l.mtz", (13469,1065,2323)),
        ("less_free_l.mtz", (14002,532,2323))]:
    s = StringIO()
    mtz_obj = iotbx_mtz.object(file_name=file_name)
    mtz_obj.show_column_data(out=s, format="spreadsheet")
    s = s.getvalue()
    spreadsheets.append(s)
    n_0, n_1, n_rest = 0, 0, 0
    for line in s.splitlines()[1:]:
      if (line.endswith(",0")): n_0 += 1
      elif (line.endswith(",1")): n_1 += 1
      else:
        assert line.endswith(",")
        n_rest += 1
    assert mtz_obj.n_reflections() == n_0 + n_1 + n_rest
    assert (n_0, n_1, n_rest) == expected
  assert not show_diff(spreadsheets[0], spreadsheets[2])
  assert spreadsheets[0] != spreadsheets[1]

def run():
  exercise()
  print format_cpu_times()

if (__name__ == "__main__"):
  run()
