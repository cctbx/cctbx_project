"""reduce2 prints the peak memory usage at the end of its run, just before
the "Job complete" block printed by run_program. Other programs do not.
"""
from __future__ import absolute_import, division, print_function
from six.moves import cStringIO as StringIO
from libtbx.utils import multi_out
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2, pdbtools

pdb_str = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      11.204   2.315  10.312  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.555   2.835  10.528  1.00 20.00           C
ATOM      3  C   ALA A   1      13.588   1.726  10.381  1.00 20.00           C
ATOM      4  O   ALA A   1      13.325   0.567  10.721  1.00 20.00           O
ATOM      5  CB  ALA A   1      12.866   3.965   9.556  1.00 20.00           C
END
"""

def run(program_class, prefix, extra_args=()):
  fn_in = prefix + ".pdb"
  with open(fn_in, "w") as f:
    f.write(pdb_str)
  log = StringIO()
  logger = multi_out()
  logger.register("log", log)
  run_program(program_class=program_class, logger=logger,
    args=[fn_in, "overwrite=True"] + list(extra_args))
  return log.getvalue().splitlines()

def test_reduce2_peak_memory():
  lines = run(reduce2.Program, "tst_reduce2_peak_memory_1")
  i = lines.index("Job complete")
  # peak memory line, blank line, "=====" line, "Job complete"
  line = lines[i-3]
  fields = line.split()
  assert fields[:2] == ["peak", "memory:"], lines[i-4:i+1]
  assert float(fields[2]) > 0, line
  assert fields[3] in ["MB", "GB"], line

def test_other_program_no_peak_memory():
  lines = run(pdbtools.Program, "tst_reduce2_peak_memory_2",
    extra_args=["output.prefix=tst_reduce2_peak_memory_2_out"])
  assert lines[-1].startswith("wall clock time:"), lines[-3:]
  assert not [l for l in lines if l.startswith("peak memory:")]

if __name__ == "__main__":
  test_reduce2_peak_memory()
  test_other_program_no_peak_memory()
  print("OK")
