"""
Placement Checker for PHENIX AI Agent.

Determines whether a model is already placed in the unit cell by comparing
unit cell dimensions between the model and the reflection data / map.

This is Tier 1 of the three-tier placement detection system:

    Tier 1 (this module) — unit cell comparison, free and immediate
    Tier 2 — existing _has_placed_model() heuristics
    Tier 3 — diagnostic probe (phenix.model_vs_data / phenix.map_correlations)

Fail-safe contract
------------------
Every public function returns False / None on any error.  A parse failure
never triggers a mismatch — the workflow simply falls through to Tier 2/3.

Usage::

    from libtbx.langchain.agent.placement_checker import (
        check_xray_cell_mismatch,
        check_cryoem_cell_mismatch,
    )

    # X-ray
    if check_xray_cell_mismatch(pdb_path, mtz_path):
        # route to molecular replacement immediately

    # Cryo-EM
    if check_cryoem_cell_mismatch(pdb_path, map_path):
        # route to docking immediately
"""

from __future__ import absolute_import, division, print_function

import os
import re


# =============================================================================
# CONSTANTS
# =============================================================================

# Fraction tolerance for cell parameter comparison (5 %).
# Loose enough to accommodate cryo-EM sub-box padding; tight enough to catch
# genuine crystal-form mismatches that would require MR.
DEFAULT_TOLERANCE = 0.05


# =============================================================================
# UNIT CELL READERS
# =============================================================================

def read_pdb_unit_cell(pdb_path):
    """
    Parse the CRYST1 record from a PDB file.

    Format::

        CRYST1   57.230   57.230  146.770  90.00  90.00  90.00 P 41 21 2

    Args:
        pdb_path (str): Path to PDB or mmCIF file.

    Returns:
        tuple | None: (a, b, c, alpha, beta, gamma) as floats, or None if
        the record cannot be found or parsed.
    """
    try:
        with open(pdb_path) as fh:
            for line in fh:
                if line.startswith("CRYST1"):
                    # Columns 7-54 contain a b c alpha beta gamma
                    parts = line[6:54].split()
                    if len(parts) >= 6:
                        return tuple(float(p) for p in parts[:6])
        return None
    except Exception:
        return None


def read_mtz_unit_cell(mtz_path):
    """
    Read the crystal symmetry unit cell from an MTZ file.

    Tries two strategies in order:
      1. iotbx.mtz (preferred — no subprocess)
      2. mtzdump subprocess stdout (fallback for environments without iotbx)

    Args:
        mtz_path (str): Path to MTZ file.

    Returns:
        tuple | None: (a, b, c, alpha, beta, gamma) as floats, or None.
    """
    # Strategy 1: iotbx.mtz
    try:
        import iotbx.mtz
        mtz_obj = iotbx.mtz.object(file_name=mtz_path)
        cs = mtz_obj.crystals()[0].crystal_symmetry()
        uc = cs.unit_cell().parameters()
        return tuple(uc[:6])
    except Exception:
        pass

    # Strategy 2: mtzdump subprocess
    try:
        from libtbx import easy_run
        result = easy_run.fully_buffered(
            command=["mtzdump", "hklin", mtz_path, "END"],
            use_shell_in_subprocess=False,
        )
        for line in result.stdout_lines:
            # "* Cell Dimensions : (a=57.23, b=57.23, c=146.77, ..."
            # or plain "Cell dimensions ..."
            m = re.search(
                r'[Cc]ell\s+[Dd]imensions?\s*[:\*]?\s*'
                r'([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+'
                r'([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)',
                line,
            )
            if m:
                return tuple(float(m.group(i)) for i in range(1, 7))
    except Exception:
        pass

    return None


def read_map_unit_cells(map_path):
    """
    Run ``phenix.show_map_info`` and parse both unit cells it reports.

    The program outputs two blocks::

        Unit cell and grid for full map (one full unit cell):
           Unit cell:       (320.000, 320.000, 320.000, 90.00, 90.00, 90.00) A

        Unit cell and map for the part of map that is present:
           Map Unit cell:   (59.000, 67.000, 106.000, 90.00, 90.00, 90.00) A

    Args:
        map_path (str): Path to CCP4/MRC map file.

    Returns:
        tuple: ``(full_cell, present_cell)`` where each element is either a
        6-float tuple or None.  Both may be None if the program fails.
    """
    full_cell = None
    present_cell = None
    try:
        import shlex as _shlex
        from libtbx import easy_run
        cmd = _shlex.split("phenix.show_map_info %s" % map_path)
        result = easy_run.fully_buffered(command=cmd, use_shell_in_subprocess=False)
        lines = result.stdout_lines

        # Pattern for a parenthesised 6-float tuple
        _CELL_PAT = re.compile(
            r'\(\s*([0-9.]+)\s*,\s*([0-9.]+)\s*,\s*([0-9.]+)\s*,'
            r'\s*([0-9.]+)\s*,\s*([0-9.]+)\s*,\s*([0-9.]+)\s*\)'
        )

        for line in lines:
            m = _CELL_PAT.search(line)
            if not m:
                continue
            cell = tuple(float(m.group(i)) for i in range(1, 7))
            if "Map Unit cell" in line or "Map unit cell" in line:
                present_cell = cell
            elif "Unit cell" in line and full_cell is None:
                full_cell = cell
    except Exception:
        pass

    return full_cell, present_cell


# =============================================================================
# COMPARISON
# =============================================================================

def cells_are_compatible(cell_a, cell_b, tolerance=DEFAULT_TOLERANCE):
    """
    Compare two unit cells for compatibility.

    Returns True if all six parameters agree within ``tolerance`` (fractional).
    The comparison is symmetric and uses the larger value as the denominator
    to avoid division by zero.

    Args:
        cell_a (tuple): (a, b, c, alpha, beta, gamma)
        cell_b (tuple): (a, b, c, alpha, beta, gamma)
        tolerance (float): Fractional tolerance (default 0.05 = 5 %).

    Returns:
        bool: True if cells are compatible; False if any parameter differs
        by more than ``tolerance``.
    """
    if cell_a is None or cell_b is None:
        return True   # Can't compare → assume compatible (fail-safe)
    if len(cell_a) != 6 or len(cell_b) != 6:
        return True

    for v_a, v_b in zip(cell_a, cell_b):
        denom = max(abs(v_a), abs(v_b), 1e-6)
        if abs(v_a - v_b) / denom > tolerance:
            return False
    return True


# =============================================================================
# HIGH-LEVEL MISMATCH CHECKS
# =============================================================================

def check_xray_cell_mismatch(pdb_path, mtz_path):
    """
    Check whether model and MTZ unit cells are incompatible for X-ray.

    Returns True **only** when both cells can be read AND they are definitively
    incompatible (outside 5 % tolerance on any parameter).

    Args:
        pdb_path (str): Path to PDB/CIF model.
        mtz_path (str): Path to MTZ reflection file.

    Returns:
        bool: True → definite mismatch (model needs MR).
              False → compatible, or data unavailable (fail-safe).
    """
    if not pdb_path or not mtz_path:
        return False
    if not os.path.exists(pdb_path) or not os.path.exists(mtz_path):
        return False

    model_cell = read_pdb_unit_cell(pdb_path)
    mtz_cell = read_mtz_unit_cell(mtz_path)

    if model_cell is None or mtz_cell is None:
        return False   # Couldn't read one or both → fail-safe

    return not cells_are_compatible(model_cell, mtz_cell)


def check_cryoem_cell_mismatch(pdb_path, map_path):
    """
    Check whether model and map unit cells are incompatible for cryo-EM.

    The model is tested against **both** the full-map cell and the
    present-portion cell.  Matching either → compatible (model may have been
    placed in a sub-box extracted from the full reconstruction).

    Returns True only when both map cells can be read AND the model matches
    neither.

    Args:
        pdb_path (str): Path to PDB/CIF model.
        map_path (str): Path to CCP4/MRC map file.

    Returns:
        bool: True → definite mismatch (model needs docking).
              False → compatible, or data unavailable (fail-safe).
    """
    if not pdb_path or not map_path:
        return False
    if not os.path.exists(pdb_path) or not os.path.exists(map_path):
        return False

    model_cell = read_pdb_unit_cell(pdb_path)
    if model_cell is None:
        return False   # Can't read model cell → fail-safe

    full_cell, present_cell = read_map_unit_cells(map_path)

    # If we couldn't read either map cell → fail-safe
    if full_cell is None and present_cell is None:
        return False

    # Compatible if model matches the full map cell
    if full_cell is not None and cells_are_compatible(model_cell, full_cell):
        return False

    # Compatible if model matches the present-portion cell
    if present_cell is not None and cells_are_compatible(model_cell, present_cell):
        return False

    # Both cells readable and model matches neither → definite mismatch
    if full_cell is not None or present_cell is not None:
        # At least one cell was readable; model matched none
        return True

    return False
