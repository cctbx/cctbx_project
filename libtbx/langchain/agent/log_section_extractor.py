"""
Log Section Extractor for the Thinking Agent (v113).

Extracts scientifically informative sections from PHENIX program logs
within a character budget. Uses priority-ordered keyword tables so the
most important information is always included first.

Pure function — no LLM, no I/O, no external dependencies.
"""

from __future__ import absolute_import, division, print_function


# =========================================================================
# Per-program keyword tables
# =========================================================================
# Each entry: (section_name, [keywords])
# Listed in PRIORITY ORDER — first section is extracted first.
# If budget runs out, later sections are omitted.

SECTION_MARKERS = {

  "phenix.xtriage": [
    ("Twinning", [
      "L-test", "H-test", "Twinning fraction", "twin_fraction"]),
    ("Anomalous signal", [
      "Anomalous signal", "measurability", "d_min_anom",
      "anomalous_signal"]),
    ("Wilson statistics", [
      "Wilson plot", "Wilson", "ice_ring", "ice ring"]),
    ("Space group", [
      "Space group", "Possible space groups", "space_group_info"]),
  ],

  "phenix.phaser": [
    ("MR scores", [
      "TFZ", "LLG", "RFZ", "SOLU SET"]),
    ("Packing", [
      "Packing", "clashes"]),
    ("Search strategy", [
      "Rotation Function", "Search strategy", "copies placed"]),
  ],

  "phenix.autosol": [
    ("Phasing stats", [
      "FOM", "figure_of_merit", "sites found", "BAYES-CC"]),
    ("Density", [
      "Connectivity", "Density modification", "density_modification"]),
  ],

  "phenix.autobuild": [
    ("Building progress", [
      "Residues built", "placed residues", "Completeness"]),
    ("Map quality", [
      "Map CC", "map_correlation", "R-free"]),
  ],

  "phenix.refine": [
    ("R-factors", [
      "Start: r_work", "Final: r_work", "r_free", "r_work",
      "Refinement statistics"]),
    ("Geometry", [
      "Ramachandran", "clashscore", "RMS(bonds)", "RMS(angles)"]),
    ("Twinning", [
      "Twinning law", "twin_law"]),
    ("Difference map", [
      "Unexplained peaks", "mFo-DFc", "difference_map"]),
  ],

}

# Window sizes around keyword matches
_BEFORE = 5   # lines before match
_AFTER = 15   # lines after match
_TAIL_LINES = 100  # fallback: last N lines


def extract_sections(log_text, program_name, max_chars=3500):
  """Extract informative sections from a program log.

  Args:
    log_text:     Full program output (string).
    program_name: e.g. "phenix.xtriage", "phenix.phaser".
    max_chars:    Character budget (default 3500).

  Returns:
    String containing extracted sections with headers.
    Empty string if log_text is empty.
  """
  if not log_text or not log_text.strip():
    return ""

  lines = log_text.splitlines()
  markers = SECTION_MARKERS.get(program_name)

  if not markers:
    return _fallback_tail(lines, program_name, max_chars)

  output_parts = []
  total_chars = 0
  any_match = False

  for section_name, keywords in markers:
    # Find all matching line indices
    match_indices = set()
    for i, line in enumerate(lines):
      for kw in keywords:
        if kw in line:
          match_indices.add(i)
          break

    if not match_indices:
      continue

    any_match = True

    # Build windows around matches
    windows = []
    for idx in sorted(match_indices):
      start = max(0, idx - _BEFORE)
      end = min(len(lines), idx + _AFTER + 1)
      windows.append((start, end))

    # Merge overlapping windows
    merged = _merge_windows(windows)

    # Build section text
    header = "--- %s (%s) ---" % (section_name, program_name)
    section_lines = [header]
    for start, end in merged:
      section_lines.extend(lines[start:end])
    section_text = "\n".join(section_lines) + "\n"

    # Budget check
    if total_chars + len(section_text) > max_chars:
      # Try to fit a truncated version
      remaining = max_chars - total_chars
      if remaining > len(header) + 20:
        truncated = _truncate_at_line(section_text, remaining)
        output_parts.append(truncated)
        total_chars += len(truncated)
      if output_parts:
        output_parts.append("[remaining sections omitted]\n")
      break

    output_parts.append(section_text)
    total_chars += len(section_text)

  if not any_match:
    return _fallback_tail(lines, program_name, max_chars)

  return "\n".join(output_parts)


def _merge_windows(windows):
  """Merge overlapping (start, end) intervals."""
  if not windows:
    return []
  sorted_w = sorted(windows)
  merged = [sorted_w[0]]
  for start, end in sorted_w[1:]:
    if start <= merged[-1][1]:
      merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    else:
      merged.append((start, end))
  return merged


def _truncate_at_line(text, max_chars):
  """Truncate text to fit within max_chars, ending at a line boundary."""
  if len(text) <= max_chars:
    return text
  cut = text[:max_chars]
  last_nl = cut.rfind("\n")
  if last_nl > 0:
    return cut[:last_nl + 1]
  return cut


def _fallback_tail(lines, program_name, max_chars):
  """Return last _TAIL_LINES lines as fallback."""
  tail = lines[-_TAIL_LINES:] if len(lines) > _TAIL_LINES else lines
  header = "--- Last %d lines (%s) ---" % (len(tail), program_name)
  text = header + "\n" + "\n".join(tail) + "\n"
  if len(text) > max_chars:
    text = _truncate_at_line(text, max_chars)
  return text
