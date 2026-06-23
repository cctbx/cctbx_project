"""Behavior-equivalence regression test for cctbx.crystal.find_best_cell.find_best_cell.

find_best_cell.__init__ is the dominant runtime cost of dials.stills_process indexing (it is
reached via crystal.symmetry.change_of_basis_op_to_best_cell, called millions of times). The
find-best-cell-perf branch rewrote it with a per-space-group "plan" cache and lazy materialization
of symmetry()/all_cells(), claiming bit-identical results, but the production class had no direct
correctness coverage (the existing tst_find_best_cell.py only exercises alternative_find_best_cell).

This test pins cb_op(), symmetry(), and all_cells() to a golden baseline generated from master
(stored in find_best_cell_golden.pickle) and asserts the current code reproduces it. The baseline
is a modest, curated subset of inputs that still spans all four find_best_cell branches.

Inputs: a subset of the niggli-reduction cell population (harvested verbatim from
cctbx/regression/tst_niggli_reduction_cpp.py, including its edge cases). The niggli cells carry no
space group, so each is classified by its angle pattern into a lattice metric system and paired
only with space groups that metric actually satisfies: all-right-angle cells get orthorhombic
settings; a single non-right angle gets the monoclinic settings of the matching unique axis
(b-unique reference plus a-/c-unique non-reference, which exercise the cb_op_std_inp derivation);
cells with no usable pair of right angles get only triclinic. This keeps every (cell, space group)
pairing physically sensible -- compatibility is additionally enforced by the default
crystal.symmetry construction. Inputs are subsampled per lattice system to keep the baseline small
while covering the niggli (P 1), identity, monoclinic, and orthorhombic branches.

Usage (run through the cctbx interpreter):

  libtbx.python cctbx/regression/tst_find_best_cell_equivalence.py
      match the current code against the committed golden baseline (prints OK).

  libtbx.python cctbx/regression/tst_find_best_cell_equivalence.py --generate [--out DIR]
      regenerate the golden baseline from the CURRENTLY IMPORTED find_best_cell. To refresh it
      from master, run this on a master checkout (the committed baseline was generated that way).

With no arguments it runs the match, so it is safe as a registered regression test.
"""
from __future__ import absolute_import, division, print_function

import os
import sys
import importlib

from cctbx import sgtbx, crystal
from libtbx import easy_pickle

# Import the submodule explicitly: cctbx.crystal re-exports the find_best_cell *class* under the
# same name, which shadows the submodule for `import ... as` / attribute access.
_fbc_module = importlib.import_module("cctbx.crystal.find_best_cell")

GOLDEN_BASENAME = "find_best_cell_golden.pickle"

# An angle within this many degrees of 90 is treated as a right angle when classifying a cell's
# lattice metric system. Much tighter than crystal.symmetry's own compatibility tolerance, so a
# cell classified into a system always passes that system's construction gate.
ANGLE_TOL = 1.0e-3

# Space groups grouped by the lattice metric they require, so each cell is only ever paired with
# space groups its angles actually satisfy (see classify_lattice_system). Monoclinic groups are
# split by unique axis -- the axis whose angle differs from 90. b-unique is the reference
# monoclinic setting; a- and c-unique are non-reference settings that exercise find_best_cell's
# cb_op_std_inp derivation. Numbers: ==1 niggli; 2 identity (<3); 3-15 monoclinic; 16-74
# orthorhombic. (High-symmetry >=75 groups are omitted: they hit the SAME identity early-return as
# P -1, so add no code coverage while demanding tetragonal/hexagonal/cubic metrics.) The set is
# kept lean -- one primitive b-unique monoclinic (P 1 2 1, the common case, exercising the
# best-monoclinic-beta path) plus centered a-/c-unique settings (far fewer candidate ops, so a much
# smaller all_cells) -- to keep the golden baseline compact.
LATTICE_SPACE_GROUPS = {
  "triclinic":    ("P 1", "P -1"),
  "monoclinic_a": ("C 2 1 1",),
  "monoclinic_b": ("P 1 2 1",),
  "monoclinic_c": ("I 1 1 2",),
  "orthorhombic": ("P 2 2 2", "P 21 21 21", "C 2 2 21", "F 2 2 2"),
}

# best_monoclinic_beta only changes the result for b-unique monoclinic (do_best_beta needs
# unique_axis == 1); that group gets both flags (exercising the do_best_beta True and False paths),
# everything else the production default only.
BETA_FLAGS_BY_SYSTEM = {"monoclinic_b": (True, False)}
DEFAULT_BETA_FLAGS = (True,)

ANGULAR_TOLERANCE = 3

# Cap on cells per lattice system. Edge-case sources are always kept; the large sweeps fill any
# remaining room (strided for determinism). Capping per system (rather than striding the raw grid
# uniformly) guarantees the rare all-right-angle orthorhombic cells and each monoclinic unique-axis
# are represented. Kept small so the committed baseline stays compact.
PER_SYSTEM_CAP = {
  "triclinic": 12, "monoclinic_a": 2, "monoclinic_b": 2,
  "monoclinic_c": 2, "orthorhombic": 5,
}

# Genuine edge-case sources (always kept in full); the rest are large sweeps that fill cap room.
EDGE_SOURCES = frozenset((
  "exercise_textbook_examples", "exercise_problem_parameters", "exercise_extreme",
  "exercise_real_world_examples", "exercise_a3_a4_zero_component", "exercise_zero_epsilon"))

# The cell-producing exercise_* functions of the niggli suite (exercise_iteration_limit is
# omitted: it does not call compare_reductions and produces no cell).
CELL_SOURCES = (
  "exercise_textbook_examples",
  "exercise_problem_parameters",
  "exercise_extreme",
  "exercise_real_world_examples",
  "exercise_grid",
  "exercise_bravais_plus",
  "exercise_gruber_types",
  "exercise_a3_a4_zero_component",
  "exercise_zero_epsilon",
)


def harvest_niggli_cells():
  """Collect the canonical niggli cell population, grouped by source.

  Reuses the exact cells the niggli suite exercises (including the seeded random Gruber-type set
  and every hand-collected edge case) by monkeypatching the suite's compare_reductions hooks with
  a collector -- no duplication of the cell definitions here.
  """
  from cctbx.regression import tst_niggli_reduction_cpp as niggli
  collected = {}
  current = []
  def collector(unit_cell, *args, **kwargs):
    current.append(unit_cell)
  orig_cr = niggli.compare_reductions
  orig_crl = niggli.compare_reductions_with_limit
  niggli.compare_reductions = collector
  niggli.compare_reductions_with_limit = collector
  try:
    for source in CELL_SOURCES:
      del current[:]
      getattr(niggli, source)()
      collected[source] = list(current)
  finally:
    niggli.compare_reductions = orig_cr
    niggli.compare_reductions_with_limit = orig_crl
  return collected


def classify_lattice_system(unit_cell):
  """Classify a cell into the lattice metric system its angles support (see LATTICE_SPACE_GROUPS).

  all right angles -> orthorhombic; exactly one non-right angle -> monoclinic of the matching
  unique axis (a/b/c); otherwise -> triclinic.
  """
  alpha, beta, gamma = unit_cell.parameters()[3:6]
  off = [abs(angle - 90.0) > ANGLE_TOL for angle in (alpha, beta, gamma)]
  n_off = sum(off)
  if n_off == 0:
    return "orthorhombic"
  if n_off == 1:
    return ("monoclinic_a", "monoclinic_b", "monoclinic_c")[off.index(True)]
  return "triclinic"


def select_cells():
  """Return [(system, source, unit_cell), ...].

  Every cell is classified by its lattice metric system. Edge-case sources are always kept; the
  large sweeps fill each system's remaining cap (strided for determinism).
  """
  collected = harvest_niggli_cells()
  by_system = dict((system, {"edge": [], "large": []})
                   for system in LATTICE_SPACE_GROUPS)
  for source in CELL_SOURCES:
    bucket = "edge" if source in EDGE_SOURCES else "large"
    for uc in collected[source]:
      by_system[classify_lattice_system(uc)][bucket].append((source, uc))
  selected = []
  for system in LATTICE_SPACE_GROUPS:
    edge = by_system[system]["edge"]
    large = by_system[system]["large"]
    room = max(0, PER_SYSTEM_CAP[system] - len(edge))
    if room == 0:
      sampled = []
    elif len(large) <= room:
      sampled = large
    else:
      sampled = large[::(len(large) // room)][:room]
    for source, uc in edge + sampled:
      selected.append((system, source, uc))
  return selected


def build_inputs():
  """Return the deterministic, ordered list of inputs:
  [(label, crystal.symmetry, angular_tolerance, best_monoclinic_beta), ...].

  Each niggli cell is paired with the space groups its lattice metric supports and the applicable
  beta flag(s). crystal.symmetry is built with its default flags, so the construction itself
  asserts the cell is compatible with the space group (metric-matched by construction; this is a
  belt-and-suspenders gate). Pairs on which crystal.symmetry raises are skipped. This function does
  NOT depend on find_best_cell, so the same inputs are reproduced at generate time (on master) and
  at match time (on the branch).
  """
  inputs = []
  for system, source, uc in select_cells():
    betas = BETA_FLAGS_BY_SYSTEM.get(system, DEFAULT_BETA_FLAGS)
    for symbol in LATTICE_SPACE_GROUPS[system]:
      sgi = sgtbx.space_group_info(symbol)
      try:
        xs = crystal.symmetry(unit_cell=uc, space_group_info=sgi)
      except Exception:
        continue
      for beta in betas:
        label = "%s|%s|%s|beta=%s" % (system, source, symbol, beta)
        inputs.append((label, xs, ANGULAR_TOLERANCE, beta))
  return inputs


def canonical_outputs(fbc, sg_cache):
  """Exact, order-preserving serialization of the three public accessors.

  cb_op as its canonical xyz string; symmetry and every all_cells entry as
  (unit-cell parameters, space-group lookup_symbol). Space-group symbols are memoized in sg_cache
  (keyed on the hashable space_group object) to avoid repeated space_group_type builds.
  """
  def sg_symbol(symmetry):
    sg = symmetry.space_group()
    key = sg_cache.get(sg)
    if key is None:
      key = sg.type().lookup_symbol()
      sg_cache[sg] = key
    return key
  sym = fbc.symmetry()
  return (
    fbc.cb_op().as_xyz(),
    tuple(sym.unit_cell().parameters()),
    sg_symbol(sym),
    tuple((tuple(c.unit_cell().parameters()), sg_symbol(c))
          for c in fbc.all_cells()),
  )


def result_for(find_best_cell_cls, xs, ang_tol, beta, sg_cache):
  """Return ("ok", canonical_outputs) or ("raised", exception_type_name).

  Some inputs legitimately raise (e.g. niggli reduction of the extreme high-iteration cell exceeds
  its limit); the raising behavior is itself pinned, so a difference in it is a regression."""
  try:
    fbc = find_best_cell_cls(
      xs, angular_tolerance=ang_tol, best_monoclinic_beta=beta)
  except Exception as e:
    return ("raised", type(e).__name__)
  return ("ok", canonical_outputs(fbc, sg_cache))


def golden_path(out_dir=None):
  if out_dir is None:
    out_dir = os.path.dirname(os.path.abspath(__file__))
  return os.path.join(out_dir, GOLDEN_BASENAME)


def generate(out_dir=None):
  cls = _fbc_module.find_best_cell
  sg_cache = {}
  rows = [(label, result_for(cls, xs, ang_tol, beta, sg_cache))
          for label, xs, ang_tol, beta in build_inputs()]
  path = golden_path(out_dir)
  if out_dir is not None and not os.path.isdir(out_dir):
    os.makedirs(out_dir)
  easy_pickle.dump(path, rows)
  size_kb = os.path.getsize(path) / 1024.0
  print("wrote %d inputs -> %s (%.1f KB)" % (len(rows), path, size_kb))


def _describe_diff(label, expected, got):
  """Short description of the first divergence between two result_for() values, or None."""
  if expected[0] != got[0]:
    return "%s: kind got %r, expected %r" % (label, got[0], expected[0])
  if expected[0] == "raised":
    return None  # both raised the same exception type
  exp_cb, exp_sym_uc, exp_sym_sg, exp_all = expected[1]
  got_cb, got_sym_uc, got_sym_sg, got_all = got[1]
  if got_cb != exp_cb:
    return "%s: cb_op got %s, expected %s" % (label, got_cb, exp_cb)
  if got_sym_uc != exp_sym_uc:
    return "%s: symmetry unit cell got %s, expected %s" % (label, got_sym_uc, exp_sym_uc)
  if got_sym_sg != exp_sym_sg:
    return "%s: symmetry space group got %s, expected %s" % (label, got_sym_sg, exp_sym_sg)
  if len(got_all) != len(exp_all):
    return "%s: all_cells length got %d, expected %d" % (label, len(got_all), len(exp_all))
  for i, (g, e) in enumerate(zip(got_all, exp_all)):
    if g != e:
      return "%s: all_cells[%d] got %s, expected %s" % (label, i, g, e)
  return None


def match(out_dir=None):
  path = golden_path(out_dir)
  if not os.path.isfile(path):
    raise SystemExit(
      "golden baseline not found: %s\n"
      "  generate it (on a master checkout) with:\n"
      "    libtbx.python %s --generate" % (path, os.path.basename(__file__)))
  expected_rows = easy_pickle.load(path)
  inputs = build_inputs()
  if len(inputs) != len(expected_rows):
    raise AssertionError(
      "input set changed (%d now vs %d in the golden baseline); regenerate the baseline"
      % (len(inputs), len(expected_rows)))
  sg_cache = {}
  mismatches = []
  for (label, xs, ang_tol, beta), (exp_label, expected) in zip(inputs, expected_rows):
    if label != exp_label:
      raise AssertionError(
        "input order changed (%r vs golden %r); regenerate the baseline" % (label, exp_label))
    got = result_for(_fbc_module.find_best_cell, xs, ang_tol, beta, sg_cache)
    diff = _describe_diff(label, expected, got)
    if diff is not None:
      mismatches.append(diff)
  if mismatches:
    for d in mismatches[:20]:
      print("MISMATCH " + d)
    if len(mismatches) > 20:
      print("... and %d more" % (len(mismatches) - 20))
    raise AssertionError(
      "find_best_cell output diverged from the golden baseline in %d of %d inputs"
      % (len(mismatches), len(inputs)))
  print("OK (%d inputs reproduced the golden baseline)" % len(inputs))


def run(argv):
  if "--generate" in argv:
    out_dir = argv[argv.index("--out") + 1] if "--out" in argv else None
    generate(out_dir=out_dir)
  else:
    match()


if __name__ == "__main__":
  run(sys.argv[1:])
