"""Search a PDB unit-cell database for the entries nearest a query cell,
correctly accounting for reduction-boundary settings.

This is a demo CLI that ports the logic prototyped in search_full.py (see
/Users/dwpaley/daily/20251120/search_full.py) onto the new
crystal.symmetry.change_of_basis_op_to_nearest_setting API, which fixes a
change-of-basis convention bug present in the original prototype (see the
TODO in search_full.py for details).
"""
# LIBTBX_SET_DISPATCHER_NAME cctbx.search_nearest_cell

from __future__ import absolute_import, division, print_function

import sys
import time
from six.moves import cStringIO as StringIO

import numpy as np

import iotbx.phil
from cctbx import crystal
from cctbx import sgtbx
from libtbx.utils import Sorry

master_phil_str = """
unit_cell = None
  .type = unit_cell
  .help = "Query unit cell, e.g. unit_cell=39.741,183.767,140.649,90,90,90"
space_group = None
  .type = space_group
  .help = "Query space group, e.g. space_group=C222"
data_file = None
  .type = path
  .help = "Augmented PDB cell CSV (15 columns: pdb_id,a,b,c,alpha,beta,gamma,Z,space_group,min_cell_a,min_cell_b,min_cell_c,min_cell_alpha,min_cell_beta,min_cell_gamma), written by cctbx.fetch_pdb_cells, e.g. data/pdb_crystallography_data.csv"
n_results = 10
  .type = int
  .help = "Number of top matches to display"
length_tolerance = 0.03
  .type = float
  .help = "Fractional length tolerance passed to nearest_setting"
angle_tolerance = 3.0
  .type = float
  .help = "Angle tolerance (degrees) passed to nearest_setting"
max_search_index = 10
  .type = int
  .help = "Currently unused: change_of_basis_op_to_nearest_setting does not" \
          "expose this knob; the internal find_near_minimum_settings default applies"
"""

master_phil = iotbx.phil.parse(master_phil_str)

usage_string = """\
cctbx.search_nearest_cell unit_cell=39.741,183.767,140.649,90,90,90 \\
  space_group=C222 data_file=pdb_crystallography_data.csv n_results=10

Search a PDB unit-cell database (CSV) for the entries nearest a query cell.
The nearly-reduced settings of the query cell are computed once; each database
cell is then reduced to its minimum cell and compared against all of those
settings (via crystal.symmetry.change_of_basis_op_to_nearest_setting), so that
matches which are only visible across a cell-reduction boundary are not missed.
Each match is reported transformed into the query's setting."""


def load_pdb_cells(path):
  """Parse an augmented PDB cell CSV (15 columns, written by
  cctbx.fetch_pdb_cells) into (rows, n_skipped, min_cell_params).

  rows: list of (pdb_code, crystal.symmetry) built from columns 0-8, same
    semantics as before (skip-and-count on unparsable pdb_id/cell/space
    group).
  n_skipped: int, same semantics as before.
  min_cell_params: (N, 6) float ndarray, N == len(rows), row-aligned with
    `rows`. A row whose min-cell columns (9-14) are blank (precompute-time
    reduction failure) has all 6 entries as np.nan; the caller derives
    (min_cell_params_valid, valid_indices) via
    `~np.isnan(min_cell_params).any(axis=1)`.

  Assumes the augmented 15-column format unconditionally -- there is no
  legacy plain-CSV fallback. A file lacking the 6 trailing min-cell columns
  (e.g. an old 9-column data/pdb_crystallography_data.csv snapshot) will
  raise IndexError on the first data row, since the min-cell fields are read
  by direct positional indexing (entry[9]..entry[14]); regenerate such a
  file with cctbx.fetch_pdb_cells.
  """
  rows = []
  n_skipped = 0
  min_cell_rows = []
  sg_info_cache = {}
  with open(path) as f:
    for i, line in enumerate(f):
      if i == 0:
        continue  # header row
      entry = line.strip().split(',')
      try:
        pdb_code = entry[0]
        a, b, c, al, be, ga = [float(x) for x in entry[1:7]]
        sg_str = entry[8]
        if sg_str not in sg_info_cache:
          sg_info_cache[sg_str] = sgtbx.space_group_info(symbol=sg_str)
        # Mirrors fetch_pdb_cells.py's _crystal_symmetry_from_row (own inline
        # copy, not shared -- see tst_move_nearest.py's equivalence test that
        # pins the two together).
        try:
          cs = crystal.symmetry(
            unit_cell=(a, b, c, al, be, ga),
            space_group_info=sg_info_cache[sg_str])
        except Exception:
          sg_str_r = sg_str + ' :R'
          if sg_str_r not in sg_info_cache:
            sg_info_cache[sg_str_r] = sgtbx.space_group_info(symbol=sg_str_r)
          cs = crystal.symmetry(
            unit_cell=(a, b, c, al, be, ga),
            space_group_info=sg_info_cache[sg_str_r])
      except Exception:
        n_skipped += 1
        continue
      rows.append((pdb_code, cs))
      # Direct indexing (not a entry[9:15] slice) so a short/truncated data
      # row raises IndexError immediately rather than silently producing a
      # shorter-than-6 (or all-NaN) min-cell tuple.
      min_cell_entry = []
      for j in (9, 10, 11, 12, 13, 14):
        x = entry[j]
        try:
          min_cell_entry.append(float(x))
        except Exception:
          min_cell_entry.append(np.nan)
      min_cell_rows.append(tuple(min_cell_entry))
  min_cell_params = np.array(min_cell_rows, dtype=float) if min_cell_rows \
    else np.empty((0, 6))
  return rows, n_skipped, min_cell_params


def run(args, out=sys.stdout):
  if (len(args) == 0):
    params_out = StringIO()
    master_phil.show(out=params_out, prefix="  ")
    print(usage_string, file=out)
    print("", file=out)
    print("Full parameters:", file=out)
    print(params_out.getvalue(), file=out)
    return

  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_phil_str,
    space_group_def="space_group",
    unit_cell_def="unit_cell",
    usage_string=usage_string)
  params = cmdline.work.extract()

  if (params.unit_cell is None):
    raise Sorry("Query unit_cell not specified.")
  if (params.space_group is None):
    raise Sorry("Query space_group not specified.")
  if (params.data_file is None):
    raise Sorry("data_file not specified.")

  cs_ref = crystal.symmetry(
    unit_cell=params.unit_cell,
    space_group_info=params.space_group)

  print("Query: unit_cell=%s  space_group=%s" % (
    cs_ref.unit_cell(), cs_ref.space_group_info()), file=out)

  rows, n_skipped, precomputed = load_pdb_cells(params.data_file)
  print("Loaded %d cells from %s (%d skipped: unparsable space group or cell)" % (
    len(rows), params.data_file, n_skipped), file=out)

  print("Comparing nearly-reduced settings of the query against all database cells...",
        file=out)
  t0 = time.time()

  # Bulk architecture, ranking in the query's MINIMUM-CELL frame to match
  # change_of_basis_op_to_nearest_setting's own internal selection criterion
  # (it picks its best setting by minimizing cell_distance(other's minimum
  # cell, setting['cell']) over self's cached nearly-reduced settings -- see
  # crystal.symmetry.change_of_basis_op_to_nearest_setting): (1) get the
  # query's S nearly-reduced settings once, as an (S,6) array of their cell
  # parameters; (2) get every database row's minimum-cell params via the
  # lean path (no change_of_basis_op object construction), as an (N,6)
  # array; (3) for each of the S settings, take the elementwise L1 distance
  # to all N rows and keep a running elementwise minimum -- no (S,N,...)
  # intermediate is ever built; (4) run the unmodified, exact
  # change_of_basis_op_to_nearest_setting API on only the top-M candidates
  # by this minimum-cell-frame distance for the final ranking.
  settings, cbi_near_ops = cs_ref.near_minimum_settings_and_cb_ops(
    length_tolerance=params.length_tolerance,
    angle_tolerance=params.angle_tolerance,
    test_multiples=False)
  query_settings = np.array([s['cell'] for s in settings], dtype=float)  # (S,6)

  # min-cell params are precomputed (see fetch_pdb_cells.py) and loaded
  # directly from the CSV; a row whose min-cell columns were blank at
  # precompute time (reduction failure) is NaN here and excluded via
  # valid_mask, mirroring the old per-row try/except loop's behavior.
  valid_mask = ~np.isnan(precomputed).any(axis=1)
  min_cell_params = precomputed[valid_mask]
  valid_indices = np.nonzero(valid_mask)[0]

  best_dist = np.full(min_cell_params.shape[0], np.inf)
  for setting_row in query_settings:
    d = np.abs(min_cell_params - setting_row).sum(axis=1)
    np.minimum(best_dist, d, out=best_dist)

  margin = max(100, 10 * params.n_results)
  top_idx = np.argpartition(best_dist, margin)[:margin] \
    if margin < len(best_dist) else np.arange(len(best_dist))
  top_idx = top_idx[np.argsort(best_dist[top_idx])]
  candidate_order = valid_indices[top_idx]
  candidate_dist = best_dist[top_idx]

  results = []
  for idx, dist in zip(candidate_order, candidate_dist):
    pdb_code, cs_row = rows[idx]
    try:
      cb_op = cs_ref.change_of_basis_op_to_nearest_setting(
        cs_row,
        length_tolerance=params.length_tolerance,
        angle_tolerance=params.angle_tolerance,
        test_multiples=False)
    except Exception:
      continue  # search failure, same bucket as parse failures
    transformed_uc = cs_row.unit_cell().change_basis(cb_op)
    reindexed = not cb_op.is_identity_op()
    results.append((dist, pdb_code, transformed_uc.parameters(), reindexed))
  elapsed = time.time() - t0
  print("done (%.1fs)" % elapsed, file=out)

  results.sort(key=lambda r: r[0])
  top = results[:params.n_results]

  print("", file=out)
  print("%-6s%-8s%-9s%9s%10s%10s%9s%9s%9s   %s" % (
    "Rank", "PDB", "Distance", "a", "b", "c", "alpha", "beta", "gamma",
    "Reindexed"), file=out)
  for rank, (dist, pdb_code, params_uc, reindexed) in enumerate(top, start=1):
    a, b, c, al, be, ga = params_uc
    print("%-6d%-8s%-9.3f%9.3f%10.3f%10.3f%9.3f%9.3f%9.3f   %s" % (
      rank, pdb_code, dist, a, b, c, al, be, ga,
      "yes" if reindexed else "no"), file=out)


if (__name__ == "__main__"):
  run(sys.argv[1:])
