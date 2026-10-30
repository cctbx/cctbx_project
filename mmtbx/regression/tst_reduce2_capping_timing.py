"""Manual timing script (not part of run_tests.py) for reduce2 capping.

Performance target: on a ribosome-sized model, detection + capping together
(capping.find_dangling_atoms + capping.cap_dangling_atoms) must take less
than one restraint build (model.process(make_restraints=True)).

Usage, from an empty scratch directory (the model is downloaded there).
Step (d) dominates the runtime (hours on a ribosome); pass skip_d to omit it:
  libtbx.python tst_reduce2_capping_timing.py [PDB_CODE] [skip_d]   (default 4ybb)

Measures, in seconds: (a) model.manager + process(make_restraints=True),
(b) find_dangling_atoms, (c) cap_dangling_atoms on a second, identically
built model (its find_dangling_atoms is not timed), (d) a full
reduce2 qmready=True capping=True output.write_files=False run. Prints PASS
if (b)+(c) < (a), else SLOW; FAIL if (c) raised.

Result (2026-10-23, Intel Core i9-9880H 2.3 GHz, macOS, 4ybb, 295060 atoms;
two instances of the script overlapped on this machine, so times are
somewhat pessimistic):
  (a) manager + process(make_restraints=True)  164.12 s
  (b) find_dangling_atoms                       13.30 s  (30 dangling atoms)
  (c) cap_dangling_atoms                         9.41 s  (36 H added)
  (d) reduce2 qmready=True capping=True: not completed, Sorry after
      6585.85 s: Bad clash(es) prevents qmready=True: pdb=" HE  ARGDP  25 "
      - pdb=" H2  HOHDP 316 " 0.459 A (an H-H clash between Arg HE and
      water H2, neither a cap H by name; not investigated)
  (b)+(c) = 22.71 s < (a): PASS
"""
from __future__ import absolute_import, division, print_function
import os
import sys
import time
import iotbx.pdb
import mmtbx.model
from iotbx.pdb.fetch import fetch_and_write
from iotbx.cli_parser import run_program
from libtbx.utils import null_out, Sorry
from mmtbx.hydrogens import capping
from mmtbx.programs import reduce2

def get_file(code):
  name = "%s.cif" % code.lower()
  if not os.path.isfile(name):
    print("Fetching %s" % code)
    name = fetch_and_write(code.lower(), entity="model_cif", log=sys.stdout)
    if name is None:
      raise Sorry("Could not fetch %s" % code)
  return os.path.abspath(name)

def build_model(path):
  inp = iotbx.pdb.input(file_name=path)
  model = mmtbx.model.manager(model_input=inp, log=null_out(),
    stop_for_unknowns=False)
  model.process(make_restraints=True)
  return model

def run(code, skip_d=False):
  path = get_file(code)
  print("Model: %s" % path)
  t = time.time()
  model = build_model(path)
  t_a = time.time() - t
  print("(a) model.manager + process(make_restraints=True): %.2f s" % t_a)
  print("    atoms: %d" % model.get_number_of_atoms())

  t = time.time()
  result = capping.find_dangling_atoms(model)
  t_b = time.time() - t
  print("(b) find_dangling_atoms: %.2f s" % t_b)

  model2 = build_model(path)
  result2 = capping.find_dangling_atoms(model2)
  print("    dangling atoms found: %d" % len(result2.records))
  t = time.time()
  c_failed = False
  try:
    capped = capping.cap_dangling_atoms(model2, result2)
    t_c = time.time() - t
    print("(c) cap_dangling_atoms: %.2f s" % t_c)
    print("    H added: %d" % capped.n_added)
  except Sorry as e:
    t_c = time.time() - t
    c_failed = True
    print("(c) cap_dangling_atoms raised Sorry after %.2f s: %s" % (t_c, e))

  t = time.time()
  if skip_d:
    print("(d) reduce2 run skipped")
  else:
    try:
      run_program(program_class=reduce2.Program, args=[path, "qmready=True",
        "capping=True", "output.write_files=False", "overwrite=True"],
        logger=null_out())
      print("(d) reduce2 qmready=True capping=True: %.2f s" % (
        time.time() - t))
    except Sorry as e:
      print("(d) reduce2 not completed (Sorry after %.2f s): %s" % (
        time.time() - t, e))

  if c_failed:
    verdict = "FAIL"
  else:
    verdict = "PASS" if (t_b + t_c) < t_a else "SLOW"
  print("(b)+(c) = %.2f s vs (a) = %.2f s: %s" % (t_b + t_c, t_a, verdict))

if __name__ == "__main__":
  args = sys.argv[1:]
  skip_d = "skip_d" in args
  codes = [a for a in args if a != "skip_d"]
  run(codes[0] if codes else "4ybb", skip_d=skip_d)
