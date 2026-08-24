"""Fetch unit-cell and space-group data for every current PDB entry from the
RCSB Data API, and write an augmented CSV (the 9 legacy columns plus 6
precomputed minimum-cell columns) usable directly by
cctbx.search_nearest_cell's data_file parameter.

This is a stdlib-only port of the network logic prototyped in
/Users/dwpaley/daily/20251120/data/fetch_all.py (which used
requests/pandas/tqdm), matching the no-pandas/no-requests convention already
used elsewhere in cctbx/command_line/.
"""
# LIBTBX_SET_DISPATCHER_NAME cctbx.fetch_pdb_cells

from __future__ import absolute_import, division, print_function

import csv
import json
import sys
import time
import urllib.error
import urllib.request
from six.moves import cStringIO as StringIO

import iotbx.phil
from cctbx import crystal
from cctbx import sgtbx
from libtbx.utils import Sorry

master_phil_str = """
output_file = pdb_crystallography_data.csv
  .type = path
  .help = "Path to write the augmented 15-column PDB cell CSV"
limit = None
  .type = int
  .help = "Test/smoke mode: cap the number of PDB IDs fetched"
batch_size = 1000
  .type = int
  .help = "GraphQL batch size"
request_delay = 0.1
  .type = float
  .help = "Seconds to sleep between batches (politeness delay)"
holdings_url = "https://data.rcsb.org/rest/v1/holdings/current/entry_ids"
  .type = str
  .help = "RCSB holdings API endpoint listing all current PDB IDs"
graphql_url = "https://data.rcsb.org/graphql"
  .type = str
  .help = "RCSB GraphQL API endpoint"
"""

master_phil = iotbx.phil.parse(master_phil_str)

usage_string = """\
cctbx.fetch_pdb_cells output_file=pdb_crystallography_data.csv limit=100

Fetch unit-cell and space-group data for every current PDB entry from the
RCSB Data API, precompute each entry's minimum-cell parameters, and write an
augmented 15-column CSV (the 9 legacy columns plus 6 minimum-cell columns)
directly usable by cctbx.search_nearest_cell's data_file parameter."""


def get_all_pdb_ids(holdings_url):
  """Fetch the list of all current PDB IDs from the RCSB holdings API."""
  with urllib.request.urlopen(holdings_url) as response:
    return json.loads(response.read())


def create_graphql_query(pdb_ids):
  """Create a GraphQL query for unit cell and space group data."""
  ids_string = json.dumps(pdb_ids)
  query = f"""
    {{
      entries(entry_ids: {ids_string}) {{
        rcsb_id
        cell {{
          length_a
          length_b
          length_c
          angle_alpha
          angle_beta
          angle_gamma
          Z_PDB
        }}
        symmetry {{
          space_group_name_H_M
        }}
      }}
    }}
    """
  return query


def fetch_batch_data(graphql_url, pdb_ids):
  """Fetch crystallographic data for a batch of PDB IDs.

  Unlike the original requests-based version, urllib.request.urlopen raises
  urllib.error.HTTPError on a non-2xx response instead of returning it as a
  normal response object with a .status_code -- so the HTTP-error handling
  is a try/except around the urlopen call, not a status check after it.
  """
  query = create_graphql_query(pdb_ids)
  body = json.dumps({"query": query}).encode()
  request = urllib.request.Request(
    graphql_url, data=body, headers={"Content-Type": "application/json"},
    method="POST")
  try:
    with urllib.request.urlopen(request) as response:
      data = json.loads(response.read())
  except urllib.error.HTTPError as e:
    print(f"Error: HTTP {e.code}")
    print(e.read())
    return []

  if "errors" in data:
    print(f"GraphQL errors: {data['errors']}")
    return []

  return data.get("data", {}).get("entries", [])


def process_entry(entry):
  """Process a single entry into a flat dictionary for CSV.
  Returns None if unit cell or space group data is not available."""

  cell = entry.get("cell")
  if not cell:
    return None

  symmetry = entry.get("symmetry")
  if not symmetry:
    return None

  space_group = symmetry.get("space_group_name_H_M", "")
  if not space_group:
    return None

  result = {
    "pdb_id": entry.get("rcsb_id", ""),
    "unit_cell_a": cell.get("length_a"),
    "unit_cell_b": cell.get("length_b"),
    "unit_cell_c": cell.get("length_c"),
    "unit_cell_alpha": cell.get("angle_alpha"),
    "unit_cell_beta": cell.get("angle_beta"),
    "unit_cell_gamma": cell.get("angle_gamma"),
    "unit_cell_Z": cell.get("Z_PDB"),
    "space_group": space_group
  }

  return result


def fetch_all_data(pdb_ids, graphql_url, batch_size, request_delay, out=sys.stdout):
  """Fetch all crystallographic data in batches, returning a plain list of
  process_entry()-shaped dicts (no pandas.DataFrame)."""
  all_results = []
  skipped_count = 0

  total_batches = (len(pdb_ids) + batch_size - 1) // batch_size
  print(f"Fetching data in batches of {batch_size} ({total_batches} batches)...",
        file=out)

  for batch_num, i in enumerate(range(0, len(pdb_ids), batch_size), start=1):
    batch = pdb_ids[i:i + batch_size]

    try:
      entries = fetch_batch_data(graphql_url, batch)
      for entry in entries:
        processed = process_entry(entry)
        if processed is not None:
          all_results.append(processed)
        else:
          skipped_count += 1

      print(f"Progress: batch {batch_num}/{total_batches} "
            f"({len(all_results)} fetched, {skipped_count} skipped)", file=out)

      # Be nice to the API - small delay between batches
      time.sleep(request_delay)

    except Exception as e:
      print(f"Error processing batch starting at index {i}: {e}", file=out)
      continue

  print(f"Skipped {skipped_count} entries without crystallographic data", file=out)
  return all_results


def _crystal_symmetry_from_row(a, b, c, alpha, beta, gamma, space_group_symbol,
                                sg_info_cache):
  """Build a crystal.symmetry from 6 cell floats + a space-group symbol string,
  retrying against the explicit rhombohedral setting (symbol + ' :R') if the bare
  symbol is rejected. sg_info_cache is a caller-owned dict mutated in place to
  memoize sgtbx.space_group_info(symbol=...) across many rows sharing a symbol.
  Raises on failure (unparsable symbol / cell incompatible with either setting).
  Mirrors search_nearest_cell.py's load_pdb_cells inline block -- kept as a separate
  copy per the command_line/ no-shared-helper-modules convention; see
  tst_move_nearest.py's equivalence test."""
  sg_str = space_group_symbol
  if sg_str not in sg_info_cache:
    sg_info_cache[sg_str] = sgtbx.space_group_info(symbol=sg_str)
  try:
    return crystal.symmetry(
      unit_cell=(a, b, c, alpha, beta, gamma),
      space_group_info=sg_info_cache[sg_str])
  except Exception:
    sg_str_r = sg_str + ' :R'
    if sg_str_r not in sg_info_cache:
      sg_info_cache[sg_str_r] = sgtbx.space_group_info(symbol=sg_str_r)
    return crystal.symmetry(
      unit_cell=(a, b, c, alpha, beta, gamma),
      space_group_info=sg_info_cache[sg_str_r])


def write_cells_csv(base_rows, output_file, out=sys.stdout):
  """base_rows: iterable of dicts with keys pdb_id, unit_cell_a..gamma,
  unit_cell_Z, space_group (exactly process_entry()'s shape).
  For each row, computes the 6 min-cell columns via _crystal_symmetry_from_row +
  cs.minimum_cell().unit_cell().parameters() (blank fields on failure) and writes
  all 15 columns via csv.writer to output_file."""
  sg_info_cache = {}
  header = [
    "pdb_id", "unit_cell_a", "unit_cell_b", "unit_cell_c",
    "unit_cell_alpha", "unit_cell_beta", "unit_cell_gamma", "unit_cell_Z",
    "space_group", "min_cell_a", "min_cell_b", "min_cell_c",
    "min_cell_alpha", "min_cell_beta", "min_cell_gamma"]

  n_written = 0
  with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    for row in base_rows:
      base_fields = [
        row["pdb_id"], row["unit_cell_a"], row["unit_cell_b"], row["unit_cell_c"],
        row["unit_cell_alpha"], row["unit_cell_beta"], row["unit_cell_gamma"],
        row["unit_cell_Z"], row["space_group"]]
      try:
        cs = _crystal_symmetry_from_row(
          row["unit_cell_a"], row["unit_cell_b"], row["unit_cell_c"],
          row["unit_cell_alpha"], row["unit_cell_beta"], row["unit_cell_gamma"],
          row["space_group"], sg_info_cache)
        min_cell_fields = list(cs.minimum_cell().unit_cell().parameters())
      except Exception:
        min_cell_fields = [""] * 6
      writer.writerow(base_fields + min_cell_fields)
      n_written += 1
      if n_written % 5000 == 0:
        print(f"Wrote {n_written} rows...", file=out)

  print(f"Wrote {n_written} rows to {output_file}", file=out)


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
    usage_string=usage_string)
  params = cmdline.work.extract()

  if (params.output_file is None):
    raise Sorry("output_file not specified.")

  print(f"Fetching list of all PDB IDs from {params.holdings_url}...", file=out)
  pdb_ids = get_all_pdb_ids(params.holdings_url)
  print(f"Found {len(pdb_ids)} PDB entries", file=out)

  if (params.limit is not None):
    pdb_ids = pdb_ids[:params.limit]
    print(f"Limiting to first {len(pdb_ids)} entries", file=out)

  base_rows = fetch_all_data(
    pdb_ids, params.graphql_url, params.batch_size, params.request_delay, out)
  write_cells_csv(base_rows, params.output_file, out)


if (__name__ == "__main__"):
  run(sys.argv[1:])
