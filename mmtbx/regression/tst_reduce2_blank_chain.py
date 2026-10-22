"""reduce2 Optimizer on a model with a blank chain ID and a clique of two
interacting Movers. The clique report summed initial scores parsed by position
from each Mover's description; the blank chain ID drops a token and the run
crashed with an IndexError. Model: 2ERL Ser 34 and EOH 200 (blank chain ID);
in altloc B, Ser HG and EOH HO form a clique.
"""
from __future__ import absolute_import, division, print_function
import os, re
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2

pdb_str = """\
CRYST1   53.910   23.100   23.100  90.00 110.40  90.00 C 2           4
ATOM    493  N   SER    34      12.844   0.772  16.083  1.00  9.93           N
ATOM    494  CA  SER    34      13.557   1.963  16.575  1.00 10.18           C
ATOM    495  C   SER    34      12.950   2.566  17.788  1.00  9.15           C
ATOM    496  O   SER    34      13.608   3.309  18.555  1.00 12.19           O
ATOM    497  CB ASER    34      13.553   3.076  15.438  0.48 13.27           C
ATOM    498  CB BSER    34      13.638   2.928  15.354  0.52 14.77           C
ATOM    499  OG ASER    34      12.413   3.945  15.615  0.48 13.17           O
ATOM    500  OG BSER    34      14.493   2.380  14.346  0.52 17.06           O
HETATM  600  C1 AEOH   200      12.974   5.558  13.017  0.54 29.75           C
HETATM  601  C1 BEOH   200      14.446   4.322  12.408  0.46 32.58           C
HETATM  602  C2 AEOH   200      12.259   6.535  13.707  0.54 26.57           C
HETATM  603  C2 BEOH   200      13.905   4.879  13.576  0.46 32.29           C
HETATM  604  O  AEOH   200      12.271   5.065  11.838  0.54 27.18           O
HETATM  605  O  BEOH   200      15.921   4.304  12.511  0.46 29.83           O
END
"""

def run(model_text, prefix, extra_args=()):
  fn_in = prefix + ".pdb"
  fn_out = prefix + "_out.pdb"
  fn_desc = prefix + "_out.txt"
  with open(fn_in, "w") as f:
    f.write(model_text)
  for fn in (fn_out, fn_desc):
    if os.path.exists(fn): os.remove(fn)
  args = [fn_in, "output.filename=" + fn_out,
    "output.description_file_name=" + fn_desc, "overwrite=True"]
  run_program(program_class=reduce2.Program, logger=null_out(),
    args=args + list(extra_args))
  with open(fn_out) as f:
    model_out = f.read()
  with open(fn_desc) as f:
    return model_out, f.read()

def check_clique(desc):
  # Totals line, then one line per Mover with its initial score
  m = re.search(r"Set of 2 Movers:\s+Totals: initial score (\S+), final score \S+\n"
    r"(.*)\n(.*)\n", desc)
  assert m is not None, desc
  scores = [float(re.search(r"Initial score: (\S+)", l).group(1))
    for l in m.group(2, 3)]
  assert "BSER 34 HG" in m.group(2) + m.group(3), m.group(0)
  assert "BEOH 200 HO" in m.group(2) + m.group(3), m.group(0)
  assert abs(float(m.group(1)) - sum(scores)) < 0.006, m.group(0)

def test_add():
  check_clique(run(pdb_str, "tst_reduce2_blank_chain_1")[1])

def test_optimize():
  model_h = run(pdb_str, "tst_reduce2_blank_chain_2")[0]
  assert " HO BEOH   200" in model_h
  check_clique(run(model_h, "tst_reduce2_blank_chain_3",
    ["approach=optimize"])[1])

if __name__ == "__main__":
  test_add()
  test_optimize()
  print("OK")
