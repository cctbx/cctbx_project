"""
Shared file utility functions for PHENIX AI Agent.

This module provides common file classification and utility functions
used across multiple modules to avoid code duplication.
"""

from __future__ import absolute_import, division, print_function

import os
import re


# =============================================================================
# MTZ CLASSIFICATION
# =============================================================================

# Patterns that indicate map coefficients (phases) rather than data
MAP_COEFFS_PATTERNS = [
    'map_coeffs',      # Explicit map coefficients
    'denmod',          # Density-modified maps
    'density_mod',     # Alternative naming
]


def classify_mtz_type(filepath):
    """
    Classify MTZ file as data_mtz or map_coeffs_mtz based on filename patterns.

    MTZ files fall into two categories:
    - data_mtz: Contains measured Fobs and R-free flags (for refinement)
    - map_coeffs_mtz: Contains calculated phases (for ligand fitting, visualization)

    Classification rules:
    1. refine_NNN.mtz or refine_NNN_NNN.mtz -> map_coeffs_mtz (refine map output)
    2. *_001.mtz pattern -> map_coeffs_mtz (other numbered map outputs)
    3. Contains 'map_coeffs', 'denmod', 'density_mod' -> map_coeffs_mtz
    4. Contains '_data.mtz' or 'refinement_data' -> data_mtz
    5. Default -> data_mtz

    Args:
        filepath: Path to MTZ file (string)

    Returns:
        str: "data_mtz" or "map_coeffs_mtz"

    Examples:
        >>> classify_mtz_type("/path/to/refine_001.mtz")
        'map_coeffs_mtz'
        >>> classify_mtz_type("/path/to/refine_001_001.mtz")
        'map_coeffs_mtz'
        >>> classify_mtz_type("/path/to/7qz0_refine_001.mtz")
        'map_coeffs_mtz'
        >>> classify_mtz_type("/path/to/data.mtz")
        'data_mtz'
        >>> classify_mtz_type("/path/to/denmod_map.mtz")
        'map_coeffs_mtz'
    """
    basename = os.path.basename(filepath).lower()

    # Pattern 1: Refine output MTZ — all variants:
    #   refine_001.mtz, refine_001_001.mtz,
    #   7qz0_refine_001.mtz, 7qz0_refine_001_001.mtz
    if re.match(r'(?:.*_)?refine_\d{3}(?:_\d{3})?\.mtz$', basename):
        return "map_coeffs_mtz"

    # Pattern 2: Other numbered outputs like model_refine_001.mtz
    if re.match(r'.*_001\.mtz$', basename):
        return "map_coeffs_mtz"

    # Pattern 3: Explicit map coefficients or density-modified
    for pattern in MAP_COEFFS_PATTERNS:
        if pattern in basename:
            return "map_coeffs_mtz"

    # Pattern 4: Explicit data files (keep as data)
    if '_data.mtz' in basename or 'refinement_data' in basename:
        return "data_mtz"

    # Default: treat as data MTZ
    return "data_mtz"


def get_mtz_stage(filepath, category):
    """
    Determine the specific stage/subcategory for an MTZ file.

    Args:
        filepath: Path to MTZ file
        category: Parent category ("data_mtz" or "map_coeffs_mtz")

    Returns:
        str: Specific stage like "refine_map_coeffs", "denmod_map_coeffs", etc.
    """
    basename = os.path.basename(filepath).lower()

    if category == "data_mtz":
        if '_data.mtz' in basename:
            return "original_data_mtz"
        if 'phased' in basename:
            return "phased_data_mtz"
        return "data_mtz"

    elif category == "map_coeffs_mtz":
        if 'denmod' in basename or 'density_mod' in basename:
            return "denmod_map_coeffs"
        if 'map_coeffs' in basename and 'predict' in basename:
            return "predict_build_map_coeffs"
        if 'overall_best_map_coeffs' in basename:
            return "predict_build_map_coeffs"
        # Default for map coeffs is refine output
        return "refine_map_coeffs"

    return category


# =============================================================================
# FILE EXTENSION UTILITIES
# =============================================================================

# Common file extensions by category
EXTENSION_CATEGORIES = {
    # Structure files
    '.pdb': 'model',
    '.cif': None,  # Could be model, ligand, or data - needs context
    '.ent': 'model',

    # Reflection data
    '.mtz': None,  # Needs classify_mtz_type()
    '.sca': 'data_mtz',
    '.hkl': 'data_mtz',

    # Maps
    '.map': 'map',
    '.ccp4': 'map',
    '.mrc': 'map',

    # Sequences
    '.fa': 'sequence',
    '.fasta': 'sequence',
    '.seq': 'sequence',
    '.dat': 'sequence',  # Common for PHENIX sequence files
}


def get_category_for_extension(filepath):
    """
    Get the likely file category based on extension.

    Note: Some extensions (like .cif, .mtz) need additional context
    to determine the exact category. This function returns None for
    those cases.

    Args:
        filepath: Path to file

    Returns:
        str or None: Category name, or None if extension needs context
    """
    ext = os.path.splitext(filepath)[1].lower()
    return EXTENSION_CATEGORIES.get(ext)


def is_mtz_file(filepath):
    """Check if file is an MTZ file."""
    return filepath.lower().endswith('.mtz')


def is_model_file(filepath):
    """Check if file is a model file (PDB/mmCIF)."""
    lower = filepath.lower()
    return lower.endswith('.pdb') or lower.endswith('.ent') or \
           (lower.endswith('.cif') and 'ligand' not in lower)


def is_map_file(filepath):
    """Check if file is a map file."""
    lower = filepath.lower()
    return lower.endswith('.map') or lower.endswith('.ccp4') or \
           lower.endswith('.mrc')


def is_sequence_file(filepath):
    """Check if file is a sequence file."""
    lower = filepath.lower()
    return lower.endswith('.fa') or lower.endswith('.fasta') or \
           lower.endswith('.seq') or lower.endswith('.dat')


# =============================================================================
# PATTERN MATCHING
# =============================================================================

def matches_exclude_pattern(basename, patterns):
    """Check if basename matches any exclude/prefer pattern (word-boundary aware).

    All patterns use word-boundary semantics on the filename stem:
      'ligand' matches 'ligand.pdb', 'my_ligand.pdb', 'ligand_001.pdb'
      'ligand' does NOT match 'noligand.pdb' (no word boundary before 'ligand')
      'lig.pdb' matches 'lig.pdb' but NOT 'nolig.pdb'

    Word boundaries are: start-of-string, end-of-string, or separators
    (underscore, hyphen, dot).

    Used by:
      - CommandBuilder._find_file_for_slot (extension fallback)
      - ai_agent._find_candidate_for_slot (safety-net injection)
    """
    bn_lower = basename.lower()
    stem = bn_lower.rsplit('.', 1)[0] if '.' in bn_lower else bn_lower
    for pat in patterns:
        pat_lower = pat.lower()
        # Extract pattern stem (strip extension if present)
        pat_stem = pat_lower.rsplit('.', 1)[0] if '.' in pat_lower else pat_lower
        # If pattern has extension, require that the basename ends with it
        if '.' in pat_lower:
            pat_ext = '.' + pat_lower.rsplit('.', 1)[1]
            if not bn_lower.endswith(pat_ext):
                continue
        # Word-boundary match on stems
        if re.search(
            r'(?:^|[_\-\.])' + re.escape(pat_stem) + r'(?=[_\-\.]|$)',
            stem):
            return True
    return False
