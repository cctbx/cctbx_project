"""
File metadata tracking for the AI agent (v113).

Each output file gets a metadata dict recording its
contents, quality, and provenance. Enables queries like
"find the latest model containing a ligand" instead of
filename pattern matching.

Entry points:
  build_file_metadata(...) -> dict
  find_latest_model_with_ligand(file_metadata) -> str
  find_best_model(file_metadata) -> str
  get_model_summary(file_metadata, path) -> str

2-space indentation, 80-char line width.
"""

from __future__ import absolute_import, division, print_function

import os


def build_file_metadata(
  file_path,
  validation_result=None,
  log_metrics=None,
  program_name=None,
  cycle_number=None,
  input_model=None,
  input_data=None,
):
  """Build metadata dict for a model file.

  Uses validation_result from Phase A if available.
  Falls back to minimal metadata (just provenance)
  if validation was not run or failed.

  Args:
    file_path: str, path to the output model.
    validation_result: dict from run_validation(),
      or None.
    log_metrics: dict with r_work, r_free, etc.
    program_name: str, program that produced this.
    cycle_number: int, cycle that produced this.
    input_model: str or None, input model path.
    input_data: str or None, input data path.

  Returns:
    dict with keys: path, basename, contents,
    quality, provenance.
  """
  meta = {
    "path": os.path.abspath(file_path),
    "basename": os.path.basename(file_path),
    "contents": {},
    "quality": {},
    "provenance": {
      "produced_by": program_name,
      "cycle": cycle_number,
      "input_model": input_model,
      "input_data": input_data,
    },
  }

  # Populate from validation result (Phase A)
  if validation_result:
    mc = validation_result.get("model_contents")
    if mc:
      meta["contents"] = mc

    geom = validation_result.get("geometry", {})
    if geom:
      meta["quality"]["rama_outliers"] = (
        geom.get("rama_outliers")
      )
      meta["quality"]["clashscore"] = (
        geom.get("clashscore")
      )
      meta["quality"]["rotamer_outliers"] = (
        geom.get("rotamer_outliers")
      )

  # R-factors from log_metrics (authoritative)
  if log_metrics:
    rf = log_metrics.get("r_free")
    if rf is not None:
      meta["quality"]["r_free"] = rf
    rw = log_metrics.get("r_work")
    if rw is not None:
      meta["quality"]["r_work"] = rw

  return meta


def find_latest_model_with_ligand(file_metadata):
  """Find the most recent model containing a ligand.

  Args:
    file_metadata: dict of {path: metadata_dict}

  Returns:
    str (file path) or None.
  """
  candidates = []
  for path, meta in file_metadata.items():
    ligs = meta.get(
      "contents", {}
    ).get("ligands", [])
    cycle = meta.get(
      "provenance", {}
    ).get("cycle", 0) or 0
    if ligs:
      candidates.append((cycle, path))
  if not candidates:
    return None
  candidates.sort(reverse=True)
  return candidates[0][1]


def find_best_model(file_metadata):
  """Find the model with the lowest R-free.

  Args:
    file_metadata: dict of {path: metadata_dict}

  Returns:
    str (file path) or None.
  """
  candidates = []
  for path, meta in file_metadata.items():
    rf = meta.get("quality", {}).get("r_free")
    if rf is not None:
      candidates.append((rf, path))
  if not candidates:
    return None
  candidates.sort()
  return candidates[0][1]


def get_model_summary(file_metadata, path):
  """One-line summary of a model's contents.

  Args:
    file_metadata: dict of {path: metadata_dict}
    path: str, key into file_metadata.

  Returns:
    str like "model_005.pdb (A,B; 1x ATP;
    187 waters; R-free=0.248)"
  """
  meta = file_metadata.get(path, {})
  name = meta.get(
    "basename", os.path.basename(path)
  )

  parts = []
  c = meta.get("contents", {})

  # Chains
  chains = c.get("chains", [])
  if chains:
    parts.append(",".join(chains))

  # Ligands
  ligs = c.get("ligands", [])
  if ligs:
    lig_names = [lig["name"] for lig in ligs[:3]]
    parts.append(
      "%dx %s" % (len(ligs), "+".join(lig_names))
    )

  # Waters
  wc = c.get("waters", 0)
  if wc:
    parts.append("%d waters" % wc)

  # Quality
  q = meta.get("quality", {})
  rf = q.get("r_free")
  if rf is not None:
    parts.append("R-free=%.3f" % rf)

  if parts:
    return "%s (%s)" % (name, "; ".join(parts))
  return name


def format_file_metadata_block(file_metadata):
  """Format a compact text block for the thinking prompt.

  Shows the best model summary. Designed to fit in
  ~150 chars of context window budget.

  Args:
    file_metadata: dict of {path: metadata_dict}

  Returns:
    str, e.g. "=== FILE METADATA ===\nBest model:
    model_005.pdb (A,B; 1x ATP; 187 waters;
    R-free=0.248)"
    or "" if no metadata.
  """
  if not file_metadata:
    return ""

  best = find_best_model(file_metadata)
  if not best:
    # No R-free data; pick most recent
    latest_cycle = -1
    latest_path = None
    for path, meta in file_metadata.items():
      cycle = meta.get(
        "provenance", {}
      ).get("cycle", 0) or 0
      if cycle > latest_cycle:
        latest_cycle = cycle
        latest_path = path
    best = latest_path

  if not best:
    return ""

  summary = get_model_summary(file_metadata, best)
  return "=== FILE METADATA ===\nBest model: %s" % (
    summary
  )
