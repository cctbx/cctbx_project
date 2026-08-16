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

import iotbx.phil
from cctbx import crystal
from cctbx.uctbx.near_minimum import cell_distance
from libtbx.utils import Sorry
from tqdm import tqdm

master_phil_str = """
unit_cell = None
  .type = unit_cell
  .help = "Query unit cell, e.g. unit_cell=39.741,183.767,140.649,90,90,90"
space_group = None
  .type = space_group
  .help = "Query space group, e.g. space_group=C222"
data_file = None
  .type = path
  .help = "PDB cell CSV (pdb_id,a,b,c,alpha,beta,gamma,Z,space_group), e.g. data/pdb_crystallography_data.csv"
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
  """Parse a PDB cell CSV (pdb_id,a,b,c,alpha,beta,gamma,Z,space_group) into
  a list of (pdb_code, crystal.symmetry) tuples.

  Rows whose space group symbol or cell parameters cannot be turned into a
  valid crystal.symmetry are silently skipped; the caller is responsible for
  reporting the total skip count (per-row warnings would flood the terminal
  over ~200k rows).
  """
  rows = []
  n_skipped = 0
  with open(path) as f:
    for i, line in enumerate(f):
      if i == 0:
        continue  # header row
      entry = line.strip().split(',')
      try:
        pdb_code = entry[0]
        a, b, c, al, be, ga = [float(x) for x in entry[1:7]]
        sg_str = entry[8]
        cs = crystal.symmetry(
          unit_cell=(a, b, c, al, be, ga),
          space_group_symbol=sg_str)
      except Exception:
        n_skipped += 1
        continue
      rows.append((pdb_code, cs))
  return rows, n_skipped


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

  rows, n_skipped = load_pdb_cells(params.data_file)
  print("Loaded %d cells from %s (%d skipped: unparsable space group or cell)" % (
    len(rows), params.data_file, n_skipped), file=out)

  print("Comparing nearly-reduced settings of the query against all database cells...",
        file=out)
  t0 = time.time()
  results = []
  for pdb_code, cs_row in tqdm(rows):
    try:
      cb_op = cs_ref.change_of_basis_op_to_nearest_setting(
        cs_row,
        length_tolerance=params.length_tolerance,
        angle_tolerance=params.angle_tolerance,
        test_multiples=False)
    except Exception:
      continue  # search failure, same bucket as parse failures
    transformed_uc = cs_row.unit_cell().change_basis(cb_op)
    dist = cell_distance(cs_ref.unit_cell().parameters(), transformed_uc.parameters())
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
