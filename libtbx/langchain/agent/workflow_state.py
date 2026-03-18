"""
Workflow State Utilities for PHENIX AI Agent.

This module provides:
- File categorization by type and purpose
- History analysis (what programs have been run)
- Experiment type detection (X-ray vs Cryo-EM)
- Workflow state detection (delegating to WorkflowEngine)
- Prompt formatting for LLM

The actual workflow logic is defined in:
- knowledge/workflows.yaml (state machine)
- knowledge/file_categories.yaml (file categorization rules)
- agent/workflow_engine.py (YAML interpreter)
"""

from __future__ import absolute_import, division, print_function
import os
import re
import fnmatch
import logging

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Safe MTZ utility imports.
#
# classify_mtz_type and get_mtz_stage live in file_utils.py.  In some
# deployments file_utils.py may be missing or may have classify_mtz_type but
# NOT get_mtz_stage (added later).  Importing them in a single statement
# means a missing get_mtz_stage kills the import of classify_mtz_type too.
#
# Solution: import each function separately, and provide inline fallbacks
# for BOTH functions if the import fails.  This means _import_mtz_utils()
# ALWAYS returns working functions — callers never need to check for None.
# ---------------------------------------------------------------------------

def _import_mtz_utils():
    """Import classify_mtz_type and get_mtz_stage with inline fallbacks.

    ALWAYS returns two callable functions.  Priority:
    1. Import from libtbx.langchain.agent.file_utils
    2. Import from agent.file_utils
    3. Inline implementation (minimal regex-based classification)

    Returns (classify_mtz_type, get_mtz_stage) — both always callable.
    """
    classify_fn = None
    stage_fn = None

    # Import classify_mtz_type
    try:
        from libtbx.langchain.agent.file_utils import classify_mtz_type
        classify_fn = classify_mtz_type
    except ImportError:
        try:
            from agent.file_utils import classify_mtz_type
            classify_fn = classify_mtz_type
        except ImportError:
            pass

    # Inline fallback for classify_mtz_type
    if classify_fn is None:
        _refine_mtz_re = re.compile(r'(?:.*_)?refine_\d{3}(?:_\d{3})?\.mtz$')
        _data_mtz_re = re.compile(r'(?:.*_)?refine_\d{3}_data\.mtz$')
        _map_markers = ('map_coeffs', 'denmod', 'density_mod', 'overall_best')
        def _classify_fallback(filepath):
            bn = os.path.basename(filepath).lower()
            if _data_mtz_re.match(bn):
                return "data_mtz"
            if _refine_mtz_re.match(bn):
                return "map_coeffs_mtz"
            # Files ending in _data.mtz are always
            # data (not map coefficients), even when
            # the name contains a map marker like
            # "overall_best".  Check this BEFORE the
            # marker scan.
            if bn.endswith('_data.mtz'):
                return "data_mtz"
            if any(m in bn for m in _map_markers):
                return "map_coeffs_mtz"
            return "data_mtz"
        classify_fn = _classify_fallback

    # Import get_mtz_stage
    try:
        from libtbx.langchain.agent.file_utils import get_mtz_stage
        stage_fn = get_mtz_stage
    except ImportError:
        try:
            from agent.file_utils import get_mtz_stage
            stage_fn = get_mtz_stage
        except ImportError:
            pass

    # Inline fallback for get_mtz_stage
    if stage_fn is None:
        def _stage_fallback(filepath, category):
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
                return "refine_map_coeffs"
            return category
        stage_fn = _stage_fallback

    return classify_fn, stage_fn



# =============================================================================
# FILE VALIDITY
# =============================================================================

# CCP4/MRC format: bytes 208-211 (0-indexed) must be the ASCII string 'MAP '
# This is the standard magic-byte check for CCP4-format density maps.
_CCP4_MAGIC_OFFSET = 208
_CCP4_MAGIC_BYTES  = b'MAP '

_MAP_EXTENSIONS    = frozenset(['.ccp4', '.mrc', '.map'])
_MODEL_EXTENSIONS  = frozenset(['.pdb', '.cif'])
# Records that must start the first non-blank line of a valid PDB/CIF file
_PDB_VALID_RECORDS = frozenset([
    # Title / citation section
    'HEADER', 'OBSLTE', 'TITLE ', 'SPLIT ', 'CAVEAT', 'COMPND', 'SOURCE',
    'KEYWDS', 'EXPDTA', 'NUMMDL', 'MDLTYP', 'AUTHOR', 'REVDAT', 'SPRSDE',
    'JRNL  ', 'REMARK',
    # Primary structure
    'DBREF ', 'DBREF1', 'DBREF2', 'SEQADV', 'SEQRES', 'MODRES',
    # Heterogen (ligands, small molecules)
    'HET   ', 'HETNAM', 'HETSYN', 'FORMUL',
    # Secondary structure
    'HELIX ', 'SHEET ',
    # Connectivity annotation
    'SSBOND', 'LINK  ', 'CISPEP',
    # Miscellaneous features
    'SITE  ',
    # Crystallographic transformation
    'CRYST1', 'ORIGX1', 'ORIGX2', 'ORIGX3',
    'SCALE1', 'SCALE2', 'SCALE3',
    'MTRIX1', 'MTRIX2', 'MTRIX3',
    # Coordinate records (the most common first lines for small-molecule files)
    'MODEL ', 'ATOM  ', 'ANISOU', 'TER   ', 'HETATM', 'ENDMDL',
    # Bookkeeping
    'CONECT', 'MASTER', 'END   ',
])
_CIF_VALID_STARTS  = ('data_', '_', 'loop_', '#')


def _is_valid_file(path):
    """
    Lightweight structural validity check before categorizing a file.

    Non-existent files are passed through (return True) — the three layers
    below only apply to files that are physically on disk.  A missing file is
    not our concern here; the downstream program that tries to open it will
    fail with a clear error.  This also ensures that callers who pass
    hypothetical filename strings (e.g. in tests) are not penalised.

    Layer 1 — Size:  zero-byte files are always invalid (crashed write).
    Layer 2 — Header (map files):  CCP4/MRC files must have the 'MAP ' magic
               at byte offset 208.  A non-zero but header-corrupt map is
               rejected so that the agent never passes a broken file downstream.
    Layer 3 — First record (model files):  PDB/CIF files must begin with a
               recognised record keyword.  Catches binary data masquerading as
               coordinates.

    Returns True only when the file passes all applicable layers (or doesn't exist).
    Errors during the check are caught and treated as invalid (conservative).
    """
    # Defensive: if path is not a string, it cannot be a file path.
    # This catches client bugs where e.g. a list of half-map paths ends up
    # as a single element in available_files.
    if not isinstance(path, (str, bytes, os.PathLike)):
        return False

    # If the file doesn't exist at all, pass it through.  The validity checks
    # below are for files that ARE on disk: they catch zero-byte crashed writes
    # and corrupt headers.  Tests frequently pass hypothetical filenames; a
    # missing file isn't our problem here — the downstream program that tries to
    # open it will fail with a clear error.
    if not os.path.exists(path):
        return True

    try:
        size = os.path.getsize(path)
    except OSError:
        # Can't read size — treat as invalid (conservative)
        return False

    # Layer 1: size — zero bytes means a crashed/incomplete write
    if size == 0:
        return False

    _, ext = os.path.splitext(path.lower())

    # Layer 2: CCP4/MRC header magic bytes
    if ext in _MAP_EXTENSIONS:
        if size < _CCP4_MAGIC_OFFSET + 4:
            # File too small to contain a proper CCP4 header (1024 bytes minimum)
            return False
        try:
            with open(path, 'rb') as fh:
                fh.seek(_CCP4_MAGIC_OFFSET)
                magic = fh.read(4)
            if magic != _CCP4_MAGIC_BYTES:
                print("WARNING: %s has invalid CCP4/MRC header (magic bytes mismatch) "
                      "— excluded from categorization" % path)
                return False
        except (OSError, IOError):
            return False

    # Layer 3: PDB/CIF structural validity check
    elif ext in _MODEL_EXTENSIONS:
        try:
            with open(path, 'r', errors='replace') as fh:
                content = fh.read()
            if ext == '.cif':
                # CIF must have at least one data_ block or loop
                first_line = next(
                    (l.strip() for l in content.splitlines() if l.strip()), '')
                if not first_line.startswith(_CIF_VALID_STARTS):
                    print("WARNING: %s does not start with a valid mmCIF record "
                          "— excluded from categorization" % path)
                    return False
            else:  # .pdb
                # PDB must contain at least one ATOM or HETATM record.
                # We scan up to 2000 lines rather than checking only the first
                # line, because PDB files can have very long headers (REMARK,
                # COMPND, SEQRES, etc.) before the coordinate records.
                # 3dnd.pdb has 546 header lines before the first ATOM record.
                lines = content.splitlines()
                has_coords = any(
                    len(l) >= 6 and l[:6].upper() in ('ATOM  ', 'HETATM')
                    for l in lines[:2000]
                )
                if not has_coords:
                    print("WARNING: %s contains no ATOM/HETATM records "
                          "— excluded from categorization" % path)
                    return False
        except (OSError, IOError):
            return False

    return True


# =============================================================================
# FILE CATEGORIZATION
# =============================================================================

def _load_category_rules():
    """Load file category rules from YAML."""
    try:
        from libtbx.langchain.knowledge.yaml_loader \
          import load_file_categories
        return load_file_categories()
    except ImportError:
        try:
            from knowledge.yaml_loader \
              import load_file_categories
            return load_file_categories()
        except ImportError:
            return None


def _match_pattern(filename, pattern):
    """Check if filename matches a pattern (supports * wildcards)."""
    # Convert pattern to regex-friendly fnmatch pattern
    return fnmatch.fnmatch(filename.lower(), pattern.lower())


def _pdb_is_small_molecule(path, max_bytes=32768):
    """
    Return True if the PDB file is a small-molecule coordinate file
    (ligand, cofactor, ion, etc.) rather than a macromolecular model.

    Detection strategy (applied to the first ``max_bytes``):
      1. Count total coordinate records (ATOM + HETATM).
      2. Very small files (≤ 150 coordinate records) → small molecule.
         Real proteins have hundreds to thousands of atoms; ligands
         typically have 10–100.  Some ligand files (e.g. atp.pdb) use
         ATOM records instead of HETATM, so record type alone is not
         reliable for small files.
      3. Larger files → small molecule only if HETATM-only (no ATOM chains).

    This catches ligand files named after their PDB hetcode, e.g. atp.pdb,
    gdp.pdb, hem.pdb — names that have no 'lig' or 'ligand' substring and so
    escape all pattern-based categorization.

    Args:
        path:      Full path to the PDB file.
        max_bytes: How much of the file to read (default 32 KB — enough to
                   see ~400 lines and reliably distinguish small ligands
                   from protein models).

    Returns:
        True  → file is a small molecule
        False → file has enough ATOM records to be a polymer, or is
                unreadable, or empty
    """
    try:
        atom_count = 0
        hetatm_count = 0
        with open(path, 'r', errors='replace') as fh:
            for line in fh.read(max_bytes).splitlines():
                if line.startswith('ATOM  ') or line.startswith('ATOM '):
                    atom_count += 1
                elif line.startswith('HETATM'):
                    hetatm_count += 1
        total = atom_count + hetatm_count
        if total == 0:
            return False              # No coordinates — can't tell

        # Small files are ligands regardless of record type.
        # ATP has ~31 atoms; smallest crystallographic protein ~500+ atoms.
        if total <= 150:
            return True

        # Larger files: small molecule only if HETATM-only
        return hetatm_count > 0 and atom_count == 0
    except Exception:
        return False                  # Be conservative on read errors


def _pdb_is_protein_model(path, max_bytes=32768):
    """
    Return True if the PDB file is predominantly a macromolecular model
    (protein, DNA, RNA) as opposed to a small-molecule ligand.

    Heuristics (applied to the first ``max_bytes``):
      1. Count total coordinate records (ATOM + HETATM)
      2. Very small files (≤ 150 coordinate records) → NOT protein.
         Real proteins have hundreds to thousands of atoms; ligands
         typically have 10–100.  Some ligand files (e.g. atp.pdb) use
         ATOM records instead of HETATM, so record type alone is not
         reliable for small files.
      3. Larger files → protein if majority ATOM records.

    Unlike ``not _pdb_is_small_molecule()``, this function returns False for
    unreadable or non-existent files rather than True, making it safe for use
    as a rejection filter in ligand-slot guards where a false positive would
    incorrectly discard a valid candidate.

    Args:
        path:      Full path to the PDB file.
        max_bytes: How much of the file to read (default 32 KB — enough to
                   see ~400 ATOM lines and reliably distinguish small ligands
                   from protein models).

    Returns:
        True  → file is positively identified as a protein model
        False → file is a small molecule, unreadable, or non-existent
    """
    try:
        atom_count = 0
        hetatm_count = 0
        with open(path, 'r', errors='replace') as fh:
            for line in fh.read(max_bytes).splitlines():
                if line.startswith('ATOM  ') or line.startswith('ATOM '):
                    atom_count += 1
                elif line.startswith('HETATM'):
                    hetatm_count += 1
        total = atom_count + hetatm_count
        if total == 0:
            return False          # No coordinates at all
        # Small files are NEVER protein models.  Ligands (ATP, GDP, heme,
        # etc.) typically have 10-100 atoms and may use ATOM records.
        # The smallest protein domain used in crystallography has ~500 atoms.
        if total <= 150:
            return False
        # Larger files: protein if majority ATOM records
        return atom_count > hetatm_count
    except Exception:
        return False                # Can't read → don't reject



def _categorize_files(available_files, ligand_hints=None, files_local=True):
    """
    Categorize files by type and purpose.

    Uses rules from knowledge/file_categories.yaml.

    Args:
        available_files: List of file paths
        ligand_hints: Optional set of basenames known to be ligands
            (from client-side best_files tracker).  Used as a fallback
            when _pdb_is_small_molecule() cannot read file content
            (e.g., on a remote server where client files are not on disk).

    Returns dict with keys for BOTH:
    - Subcategories: refined, phaser_output, predicted, etc.
    - Parent categories: model, search_model, ligand, map, mtz, sequence

    Files in subcategories are automatically "bubbled up" to their parent
    semantic categories. For example:
        - refined -> also in model
        - phaser_output -> also in model
        - predicted -> also in search_model
        - processed_predicted -> also in search_model
    """
    # Filter out zero-byte, corrupt, or structurally invalid files before categorizing.
    # This ensures has_full_map / has_model etc. only reflect genuinely usable files,
    # which in turn makes not_has conditions (Category A) reliable.
    # Filter out zero-byte, corrupt, or structurally invalid files before categorizing.
    # Skip on server where client files aren't on disk (would pass everything through
    # anyway since _is_valid_file returns True for non-existent files).
    if files_local:
        available_files = [f for f in available_files if _is_valid_file(f)]

    # Try to load YAML rules
    category_rules = _load_category_rules()

    if category_rules:
        files = _categorize_files_yaml(available_files, category_rules)
    else:
        # Fallback to hardcoded rules if YAML not available
        files = _categorize_files_hardcoded(available_files, ligand_hints=ligand_hints,
                                            files_local=files_local)

    # Bubble up subcategories to their parent semantic categories
    files = _bubble_up_to_parents(files, category_rules)

    # Post-processing: Rescue "no_ligand" / "noligand" PDB files that
    # were misclassified as ligand_pdb.  Filenames like 1J4R_no_ligand.pdb
    # contain the word "ligand" but are protein models with the ligand
    # omitted — the opposite of a ligand file.  Move them from ligand
    # categories to model.  Runs after both YAML and hardcoded paths.
    #
    # Also handles the case where the hardcoded exclusion prevented the
    # file from entering ligand_pdb at all — in that case the file is
    # in pdb but not in model (no subcategory was matched).  Promote
    # it to model so BUILD can find it.
    _anti_ligand_patterns = [
        'no_ligand', 'noligand', 'sans_ligand',
        'without_ligand', 'apo_ligand',
    ]
    # Path 1: file got INTO ligand_pdb (YAML path) — rescue it
    for f in list(files.get("ligand_pdb", [])):
        bn = os.path.basename(f).lower()
        if any(pat in bn for pat in _anti_ligand_patterns):
            files["ligand_pdb"].remove(f)
            if f in files.get("ligand", []):
                files["ligand"].remove(f)
            for cat in ("model", "pdb"):
                if cat not in files:
                    files[cat] = []
                if f not in files[cat]:
                    files[cat].append(f)

    # Path 2: Promote orphaned PDB files in "pdb" to "model".
    #
    # The hardcoded categorizer puts all PDB files in "pdb", then
    # sorts them into subcategories (refined, phaser_output, etc.)
    # which bubble up to "model" or "search_model".  User-supplied
    # input files like 1aba.pdb, myprotein.pdb don't match any
    # subcategory — they stay in "pdb" alone and has_model=False
    # in PERCEIVE, so refine is never offered.
    #
    # The YAML categorizer has a "*" catch-all for unclassified_pdb
    # that bubbles to model, but the hardcoded path lacks this.
    #
    # Fix: any PDB file in "pdb" that's NOT in a model/search_model
    # subcategory AND NOT in ligand_pdb should be promoted to model.
    _all_model_subcats = set()
    for subcat, parent in SUBCATEGORY_TO_PARENT.items():
        if parent in ("model", "search_model"):
            _all_model_subcats.add(subcat)
    _ligand_set = set(files.get("ligand_pdb", []))
    _intermediate_set = set(files.get("intermediate", [])
                            + files.get("intermediate_mr", []))
    for f in list(files.get("pdb", [])):
        bn = os.path.basename(f).lower()
        if not bn.endswith('.pdb'):
            continue
        if f in _ligand_set:
            continue  # genuinely a ligand — don't promote
        if f in _intermediate_set:
            continue  # intermediate output — don't promote
        _in_subcat = any(
            f in files.get(sc, [])
            for sc in _all_model_subcats
        )
        if not _in_subcat:
            if "model" not in files:
                files["model"] = []
            if f not in files["model"]:
                files["model"].append(f)

    # Post-processing: Cross-check MTZ categorization against file_utils.
    #
    # The YAML pattern-based categorizer can misclassify refine output MTZ files
    # (e.g. refine_001_001.mtz) as data_mtz if the patterns are stale or
    # incomplete.  The regex in file_utils.classify_mtz_type() is the canonical
    # authority — it handles all known refine output naming variants.
    #
    # When a file is in data_mtz but classify_mtz_type says map_coeffs_mtz,
    # move it: data_mtz → map_coeffs_mtz + refine_map_coeffs.  Without this,
    # ligandfit's exclude_categories=[data_mtz] guard blocks the file even
    # from the best_files fallback path.
    try:
        classify_mtz_type, get_mtz_stage = _import_mtz_utils()

        for f in list(files.get("data_mtz", [])):
            if not f.lower().endswith('.mtz'):
                continue
            canonical = classify_mtz_type(f)
            if canonical == "map_coeffs_mtz":
                # Mis-categorized: move from data_mtz to map_coeffs_mtz
                files["data_mtz"].remove(f)
                # Add to parent category
                if "map_coeffs_mtz" not in files:
                    files["map_coeffs_mtz"] = []
                if f not in files["map_coeffs_mtz"]:
                    files["map_coeffs_mtz"].append(f)
                # Add to appropriate subcategory
                stage = get_mtz_stage(f, "map_coeffs_mtz")
                if stage and stage not in ("map_coeffs_mtz",):
                    if stage not in files:
                        files[stage] = []
                    if f not in files[stage]:
                        files[stage].append(f)
    except Exception:
        pass  # Non-fatal: worst case is the old behaviour

    # Post-processing: Reclassify HETATM-only PDB files that pattern matching
    # placed in 'unclassified_pdb' (and therefore 'model').
    #
    # Files named after PDB het-codes (atp.pdb, gdp.pdb, hem.pdb, …) have no
    # 'lig' or 'ligand' in their name, so they slip through to unclassified_pdb
    # → model.  Content inspection is the reliable tiebreaker: if a PDB file
    # contains only HETATM records and no ATOM chain records, it is a small
    # molecule — not a macromolecular model.
    #
    # Only files that are SOLELY in 'unclassified_pdb' are inspected here.
    # Files that have already matched a more specific model subcategory (refined,
    # phaser_output, docked, etc.) are left alone regardless of content.

    # Pre-step: rescue known program outputs that ended up in unclassified_pdb.
    # map_to_model*.pdb is the final output of phenix.map_to_model.  A partial
    # cryo-EM build can have ≤150 atoms and trip _pdb_is_small_molecule.
    # Reclassify as autobuild_output BEFORE the small-molecule check.
    for f in list(files.get("unclassified_pdb", [])):
        bn = os.path.basename(f).lower()
        if 'map_to_model' in bn:
            if "autobuild_output" not in files:
                files["autobuild_output"] = []
            if f not in files["autobuild_output"]:
                files["autobuild_output"].append(f)
            # Also ensure it's in "model" parent (may already be via bubble-up)
            if "model" not in files:
                files["model"] = []
            if f not in files["model"]:
                files["model"].append(f)

    _model_subcats = {
        "refined", "rsr_output", "phaser_output", "autobuild_output",
        "docked", "with_ligand", "ligand_fit_output", "model_cif",
    }
    for f in list(files.get("unclassified_pdb", [])):
        # Skip if already in a specific model subcategory (content check not needed)
        if any(f in files.get(sc, []) for sc in _model_subcats):
            continue
        # Primary: content inspection (works on client where files are on disk)
        is_small = _pdb_is_small_molecule(f) if files_local else False
        # Fallback: client-side hints (works on server where files are NOT on disk)
        if not is_small and ligand_hints:
            is_small = os.path.basename(f) in ligand_hints
        if is_small:
            # Move: unclassified_pdb → ligand_pdb, model → ligand
            files["unclassified_pdb"].remove(f)
            if f in files.get("model", []):
                files["model"].remove(f)
            if f in files.get("pdb", []):
                files["pdb"].remove(f)
            for lst_key in ("ligand_pdb", "ligand"):
                if lst_key not in files:
                    files[lst_key] = []
                if f not in files[lst_key]:
                    files[lst_key].append(f)

    # Post-processing: Validate ligand_pdb classification.
    #
    # The YAML categorizer (and potentially the hardcoded categorizer) may
    # misclassify protein PDB files as ligands when broad patterns match
    # filenames like 1aba.pdb, 3gx5.pdb, etc.  A real protein with a few
    # HETATM ligand/cofactor atoms is still a macromolecular model, not a
    # ligand coordinate file.
    #
    # Content inspection is the reliable tiebreaker: if a PDB file in
    # ligand_pdb is predominantly ATOM records (i.e. _pdb_is_protein_model
    # returns True), it is a false positive and should be reclassified as
    # unclassified_pdb → model.
    #
    # Only runs when files are on local disk (files_local=True).
    if files_local and files.get("ligand_pdb"):
        ligand_keep = []
        ligand_rescued = []
        for f in files["ligand_pdb"]:
            if not f.lower().endswith('.pdb'):
                ligand_keep.append(f)
                continue
            if _pdb_is_protein_model(f):
                ligand_rescued.append(f)
            else:
                # File-size fallback: small-molecule ligands are typically
                # <5 KB (10-100 atoms × ~80 bytes/line).  A PDB file >10 KB
                # in ligand_pdb is almost certainly a protein model that
                # _pdb_is_protein_model failed to detect (e.g., unusual
                # formatting, large header, or read error).
                try:
                    fsize = os.path.getsize(f)
                    if fsize > 10000:
                        print("INFO: Rescuing '%s' from ligand_pdb "
                              "(file size %d > 10 KB, likely protein model)"
                              % (os.path.basename(f), fsize))
                        ligand_rescued.append(f)
                    else:
                        ligand_keep.append(f)
                except OSError:
                    ligand_keep.append(f)
        if ligand_rescued:
            files["ligand_pdb"] = ligand_keep
            print("INFO: Rescued %d protein model(s) from ligand_pdb: %s"
                  % (len(ligand_rescued),
                     ", ".join(os.path.basename(f) for f in ligand_rescued)))
            # Remove from ligand parent
            for f in ligand_rescued:
                if f in files.get("ligand", []):
                    files["ligand"].remove(f)
            # Add to unclassified_pdb + pdb + model
            for dest in ("unclassified_pdb", "pdb", "model"):
                if dest not in files:
                    files[dest] = []
                for f in ligand_rescued:
                    if f not in files[dest]:
                        files[dest].append(f)

    # Post-processing: Detect reference model PDB files among models.
    #
    # When 2+ PDB files are categorized as "model", one may be a
    # high-resolution reference model for restraints (not a model to
    # refine).  If auto-fill picks both as primary model inputs,
    # phenix.refine crashes with "Wrong number of models".
    #
    # Heuristic: if a model file's basename (excluding agent outputs)
    # matches keywords like "reference", "homolog", "template",
    # "restraint", "high_res", reclassify it to "reference_model".
    # The command builder won't use it as a primary model, but it
    # remains available for reference_model.file= via strategy flags.
    #
    # Guard: skip filenames that look like agent outputs (refine_*,
    # autobuild_*, etc.) to avoid false positives.
    model_files = files.get("model", [])
    if len(model_files) >= 2:
        _REF_KEYWORDS = re.compile(
            r'(reference|homolog|template|restraint|high.res)',
            re.IGNORECASE)
        _AGENT_OUTPUT_PREFIXES = (
            'refine_', 'autobuild_', 'autosol_', 'phaser_',
            'resolve_', 'predict_', 'real_space_', 'dock_',
            'map_to_model_', 'pdbtools_', 'molprobity_', 'rsr_',
        )
        for f in list(model_files):  # copy for safe removal
            bn = os.path.basename(f).lower()
            if any(bn.startswith(p) for p in _AGENT_OUTPUT_PREFIXES):
                continue
            if _REF_KEYWORDS.search(bn):
                model_files.remove(f)
                files.setdefault("reference_model", []).append(f)
                print("INFO: Reclassified '%s' from model "
                      "to reference_model (filename keyword)"
                      % os.path.basename(f))
                break  # Only reclassify one file

    # Post-processing: Validate half-map classification.
    #
    # Both the YAML and hardcoded categorizers can misclassify sequentially
    # numbered maps (map_1.ccp4, map_2.ccp4 from resolve_cryo_em segmentation)
    # as half-maps if broad patterns like *_[12].* are used.
    #
    # Half-maps from EMDB, RELION, CryoSPARC, and other cryo-EM software
    # always contain 'half' in their filename.  Any file in half_map without
    # 'half' in the name is reclassified as full_map.
    if "half_map" in files and files["half_map"]:
        reclassified = []
        keep = []
        for f in files["half_map"]:
            if 'half' in os.path.basename(f).lower():
                keep.append(f)
            else:
                reclassified.append(f)
        if reclassified:
            files["half_map"] = keep
            if "full_map" not in files:
                files["full_map"] = []
            for f in reclassified:
                if f not in files["full_map"]:
                    files["full_map"].append(f)
            # Also ensure reclassified files are in 'map' parent
            if "map" not in files:
                files["map"] = []
            for f in reclassified:
                if f not in files["map"]:
                    files["map"].append(f)

    # Post-processing: Detect half-map pairs among full_map files.
    #
    # Some cryo-EM half-maps are named with _1/_2 suffixes instead of
    # containing 'half' (e.g. 7mjs_23883_H_1.ccp4, 7mjs_23883_H_2.ccp4,
    # 7n8i_24237_box_1.ccp4, 7n8i_24237_box_2.ccp4).
    # The strict 'half'-in-name check (above) correctly avoids false
    # positives on segmented maps, but misses these legitimate pairs.
    #
    # Two-tier heuristic:
    # Tier 1: exactly 2 full_map files whose basenames differ only by
    #   a trailing _1/_2 before the extension, AND no half_maps.
    #   When these are the ONLY map files, they're almost certainly
    #   half-maps (no companion full map to lose).
    # Tier 2: >= 3 full_map files with a _1/_2 pair AND at least one
    #   other map remaining after removing the pair (the companion
    #   full map guard prevents false positives on segmented outputs).
    if not files.get("half_map"):
        _pair_re = re.compile(
            r'^(.+)[_-]([12])\.(\w+)$', re.IGNORECASE)
        _by_prefix = {}
        for f in files.get("full_map", []):
            m = _pair_re.match(os.path.basename(f))
            if m:
                prefix = m.group(1).lower()
                _by_prefix.setdefault(prefix, []).append(f)

        n_full = len(files.get("full_map", []))
        for prefix, pair in _by_prefix.items():
            if len(pair) != 2:
                continue
            # Tier 1: exactly 2 full_map files that form a pair
            # (no companion full map needed — these ARE the only maps)
            if n_full == 2:
                for f in pair:
                    files["full_map"].remove(f)
                    files["half_map"].append(f)
                break  # Only one pair possible with 2 files
            # Tier 2: pair + companion full map
            elif n_full >= 3:
                remaining = [
                    f for f in files["full_map"]
                    if f not in pair]
                if remaining:
                    for f in pair:
                        files["full_map"].remove(f)
                        files["half_map"].append(f)

    # Post-processing: If we have exactly one half-map and no full maps,
    # treat it as a full map. Half-maps only make sense in pairs for FSC.
    # A user providing a single map (even if named like a half-map) wants to use it.
    if "half_map" in files and "full_map" in files:
        if len(files["half_map"]) == 1 and len(files["full_map"]) == 0:
            files["full_map"].append(files["half_map"][0])
            files["half_map"] = []

    # Post-processing: MTZ classification safety net.
    #
    # The authoritative MTZ classifier (file_utils.classify_mtz_type) uses
    # well-tested regexes to distinguish data_mtz from map_coeffs_mtz.
    # YAML pattern matching and the hardcoded categorizer may disagree if
    # file_categories.yaml has stale/missing patterns on the production server.
    #
    # This safety net cross-checks every MTZ file against classify_mtz_type().
    # If a file is misclassified (e.g. refine_001_001.mtz ends up in data_mtz
    # but NOT in map_coeffs_mtz), it is moved to the correct category.
    # Without this, ligandfit cannot find map coefficients after refinement.
    try:
        classify_mtz_type, get_mtz_stage = _import_mtz_utils()

        # Ensure subcategory and parent category lists exist
        for cat in ["map_coeffs_mtz", "data_mtz", "refine_map_coeffs",
                     "denmod_map_coeffs", "predict_build_map_coeffs"]:
            if cat not in files:
                files[cat] = []

        # Collect all MTZ files from all categories
        all_mtz = set()
        for cat_list in files.values():
            if isinstance(cat_list, list):
                for f in cat_list:
                    if isinstance(f, str) and f.lower().endswith('.mtz'):
                        all_mtz.add(f)

        for f in all_mtz:
            correct_type = classify_mtz_type(f)
            in_data = f in files.get("data_mtz", [])
            in_map_coeffs = f in files.get("map_coeffs_mtz", [])

            if correct_type == "map_coeffs_mtz" and not in_map_coeffs:
                # File should be in map_coeffs_mtz but isn't — fix it
                # Determine the correct subcategory
                stage = get_mtz_stage(f, "map_coeffs_mtz")
                if stage and stage in files:
                    if f not in files[stage]:
                        files[stage].append(f)
                # Add to parent category
                if f not in files["map_coeffs_mtz"]:
                    files["map_coeffs_mtz"].append(f)
                # Remove from data_mtz if it was misplaced there
                if in_data:
                    files["data_mtz"].remove(f)
                logger.warning(
                    "MTZ safety net: moved %s from data_mtz to %s/%s",
                    os.path.basename(f), stage or "map_coeffs_mtz", "map_coeffs_mtz")

            elif correct_type == "map_coeffs_mtz" and in_map_coeffs and in_data:
                # File is correctly in map_coeffs_mtz but ALSO in data_mtz.
                # This happens when the YAML categorizer's Step 1 extension
                # match puts it in data_mtz (exclude patterns didn't fire)
                # AND Step 2 pattern match puts it in a map_coeffs subcategory
                # that bubbles up to map_coeffs_mtz.
                # The command builder's exclude_categories: [data_mtz] check
                # will REJECT the file if it's in both.  Remove from data_mtz.
                files["data_mtz"].remove(f)
                logger.warning(
                    "MTZ safety net: removed %s from data_mtz (was in both "
                    "data_mtz and map_coeffs_mtz)", os.path.basename(f))

            elif correct_type == "data_mtz" and not in_data:
                # File should be in data_mtz but isn't — fix it
                if f not in files["data_mtz"]:
                    files["data_mtz"].append(f)
                # Remove from map_coeffs_mtz if it was misplaced there
                if in_map_coeffs:
                    files["map_coeffs_mtz"].remove(f)
                    # Also remove from any map_coeffs subcategories
                    for subcat in ["refine_map_coeffs", "denmod_map_coeffs",
                                   "predict_build_map_coeffs"]:
                        if f in files.get(subcat, []):
                            files[subcat].remove(f)

    except Exception:
        pass  # Safety net must not break categorization

    # Post-processing: Map extension safety net (Fix 6, v115).
    #
    # Ensure ALL .ccp4/.mrc/.map files appear in at least one map
    # category.  Some tutorials (actin_sharpen_local_resolution)
    # have map files that fall through the YAML and hardcoded
    # categorizers, producing 0 files in original_files.
    #
    # Categorize by filename pattern:
    #   'half' in name → half_map
    #   everything else → full_map
    try:
        all_maps = set()
        for cat in ("full_map", "half_map",
                    "optimized_full_map", "sharpened",
                    "map"):
            for f in files.get(cat, []):
                all_maps.add(f)

        for f in available_files:
            if not f:
                continue
            _, ext = os.path.splitext(f.lower())
            if ext not in _MAP_EXTENSIONS:
                continue
            if f in all_maps:
                continue
            # Uncategorized map file — rescue it
            bn = os.path.basename(f).lower()
            if "half" in bn:
                dest = "half_map"
            else:
                dest = "full_map"
            if dest not in files:
                files[dest] = []
            files[dest].append(f)
            # Also add to parent 'map' category
            if "map" not in files:
                files["map"] = []
            if f not in files["map"]:
                files["map"].append(f)
    except Exception:
        pass  # Safety net must not break categorization

    # Post-processing: Promote orphan map files to full_map.
    #
    # Mirrors the orphan-PDB → model promotion above.  Map files that
    # end up in the "map" parent category but NOT in any subcategory
    # (full_map, half_map, optimized_full_map) are invisible to the
    # workflow engine's has_full_map check, so real_space_refine and
    # other programs with requires_full_map=true are never offered.
    #
    # Root cause: the YAML excludes list for full_map has "*_a.*" and
    # "*_b.*" (intended for half-map suffixes) which false-positive on
    # filenames like "emd-20026_auto_sharpen_A.ccp4" (the "_A" matches
    # "*_a.*" case-insensitively).  Rather than fix every exclude
    # pattern, promote any orphan to full_map — if it's not a half-map
    # or optimized map, it's a full map by elimination.
    _map_subcats = {"full_map", "half_map", "optimized_full_map"}
    _in_subcat = set()
    for sc in _map_subcats:
        for f in files.get(sc, []):
            _in_subcat.add(f)
    for f in list(files.get("map", [])):
        if f not in _in_subcat:
            if "full_map" not in files:
                files["full_map"] = []
            if f not in files["full_map"]:
                files["full_map"].append(f)

    # ── Post-categorization: phased_data_mtz (v115.08) ─────
    # Content-based promotion: check ALL data_mtz files for phase
    # columns (iotbx Tier 1 + ASCII heuristic Tier 2).  Promoted
    # files are moved from data_mtz to phased_data_mtz atomically
    # (no file exists in both categories).  Runs after BOTH the
    # YAML and hardcoded paths to handle:
    #   - YAML path: pattern-matched files need removal from data_mtz
    #   - Hardcoded path: content-based detection + promotion
    _resolve_phased_promotions(files)

    # ── Post-categorization: ignored_formats (v115.08) ─────────
    # Detect known file formats that the agent recognizes but
    # cannot use for auto-fill.  Surfaces a non-actionable WARNING
    # so the LLM and user are aware without triggering false
    # recovery attempts.
    _ignored = []
    for f in available_files:
        ext = os.path.splitext(f)[1].lower()
        if ext in _IGNORED_FORMATS:
            _ignored.append({
                'file': os.path.basename(f),
                'extension': ext,
                'note': _IGNORED_FORMATS[ext],
            })
    if _ignored:
        files['ignored_formats'] = _ignored

    return files


# ── Known but unsupported file formats (v115.08) ──────────────
# Extensions that the agent recognizes but cannot use for auto-fill.
# Kept at module level (like _PHASE_LABEL_TOKENS) to avoid
# re-creating the dict on every call.
_IGNORED_FORMATS = {
    '.cv': ('CNS cross-validation file — not supported for '
            'auto-fill in current version. Proceeding with '
            'other available data.'),
}


# Mapping from subcategory to parent semantic category
# This is the source of truth for bubbling up
SUBCATEGORY_TO_PARENT = {
    # Model subcategories (positioned, ready for refinement)
    "refined": "model",
    "rsr_output": "model",
    "phaser_output": "model",
    "autobuild_output": "model",
    "docked": "model",
    "with_ligand": "model",
    "ligand_fit_output": "ligand",  # Ligand fragment from ligandfit, not a full model
    "model_cif": "model",
    "unclassified_pdb": "model",  # Generic PDB files bubble to model category
                                  # but _has_placed_model checks subcategories to determine placement

    # Search model subcategories (templates, NOT positioned)
    "predicted": "search_model",
    "processed_predicted": "search_model",
    "pdb_template": "search_model",

    # Ligand subcategories
    "ligand_pdb": "ligand",
    "ligand_cif": "ligand",

    # Map subcategories
    "full_map": "map",
    "half_map": "map",
    "optimized_full_map": "map",
    "sharpened": "map",

    # Data MTZ subcategories (measured Fobs, R-free)
    "original_data_mtz": "data_mtz",
    "phased_data_mtz": "data_mtz",

    # Map coefficients MTZ subcategories (calculated phases)
    "refine_map_coeffs": "map_coeffs_mtz",
    "denmod_map_coeffs": "map_coeffs_mtz",
    "predict_build_map_coeffs": "map_coeffs_mtz",

    # Intermediate - these should NOT be bubbled up or tracked
    # Set to "intermediate" parent so they're excluded from model/search_model
    "intermediate_mr": "intermediate",
    "autobuild_temp": "intermediate",
    "carryover_temp": "intermediate",
}


def _bubble_up_to_parents(files, category_rules=None):
    """
    Ensure files in subcategories also appear in their parent semantic categories.

    This enables programs to request by parent category (e.g., "model")
    while files are categorized into specific subcategories (e.g., "refined").

    Args:
        files: Dict of category -> list of files
        category_rules: Optional YAML rules (for dynamic parent lookup)

    Returns:
        Updated files dict with parent categories populated
    """
    # Ensure parent categories exist
    for parent in ["model", "search_model", "ligand", "intermediate", "map", "data_mtz", "map_coeffs_mtz", "sequence"]:
        if parent not in files:
            files[parent] = []

    # Build parent mapping from YAML if available
    parent_map = dict(SUBCATEGORY_TO_PARENT)  # Start with hardcoded

    if category_rules:
        for cat_name, cat_def in category_rules.items():
            parent = cat_def.get("parent_category")
            if parent:
                parent_map[cat_name] = parent

    # Bubble up each subcategory to its parent
    for subcat, parent in parent_map.items():
        if parent is None:
            continue  # Don't bubble up if no parent
        if subcat not in files:
            continue
        if parent not in files:
            files[parent] = []

        for f in files[subcat]:
            if f not in files[parent]:
                files[parent].append(f)

    # Also ensure backward compatibility: "pdb" category contains all model and search_model
    if "pdb" not in files:
        files["pdb"] = []
    for f in files.get("model", []):
        if f not in files["pdb"]:
            files["pdb"].append(f)
    for f in files.get("search_model", []):
        if f not in files["pdb"]:
            files["pdb"].append(f)

    return files


def _categorize_files_yaml(available_files, rules):
    """
    Categorize files using YAML-defined rules.

    This handles two types of categories:
    1. Extension-based primary categories (data_mtz, map_coeffs_mtz, map, sequence)
    2. Pattern-based subcategories with semantic parents (refined->model, predicted->search_model)

    Files are first matched to subcategories by patterns, then bubbled up to parent categories.
    """
    # Initialize all categories from YAML
    files = {cat: [] for cat in rules.keys()}

    # Ensure semantic parent categories exist
    for parent in ["model", "search_model", "ligand", "map", "data_mtz", "map_coeffs_mtz", "sequence", "intermediate"]:
        if parent not in files:
            files[parent] = []

    # Also ensure pdb exists for backward compatibility
    if "pdb" not in files:
        files["pdb"] = []

    # Group categories by extension for primary matching
    # Only include categories that are PURELY extension-based (no patterns)
    # Categories with patterns are handled in Step 2
    ext_to_categories = {}
    ext_to_excludes = {}  # Track excludes for each category
    for cat_name, cat_def in rules.items():
        # Skip semantic parent categories - they don't match by extension
        if cat_def.get("is_semantic_parent"):
            continue
        # Skip categories that have patterns - they need pattern matching in Step 2
        if cat_def.get("patterns"):
            continue
        for ext in cat_def.get("extensions", []):
            if ext not in ext_to_categories:
                ext_to_categories[ext] = []
            ext_to_categories[ext].append(cat_name)
            # Store excludes for this category
            if cat_def.get("excludes"):
                ext_to_excludes[cat_name] = cat_def.get("excludes", [])

    for f in available_files:
        f_lower = f.lower()
        basename = os.path.basename(f_lower)
        _, ext = os.path.splitext(f_lower)

        # Step 1: Primary categorization by extension (for non-PDB files)
        primary_categories = ext_to_categories.get(ext, [])
        for cat in primary_categories:
            # Check excludes before adding
            excludes = ext_to_excludes.get(cat, [])
            excluded = False
            for exc_pattern in excludes:
                if _match_pattern(basename, exc_pattern):
                    excluded = True
                    break
            if not excluded and f not in files[cat]:
                files[cat].append(f)

        # Step 2: Subcategorization by patterns
        # This is where we match PDB files to specific subcategories like "refined", "predicted"
        for cat_name, cat_def in rules.items():
            # Skip semantic parent categories
            if cat_def.get("is_semantic_parent"):
                continue

            # Skip deprecated categories
            if cat_def.get("is_deprecated"):
                continue

            # Check if this category uses subcategory_of (old style) or parent_category (new semantic style)
            old_parent = cat_def.get("subcategory_of")  # e.g., refined -> pdb
            semantic_parent = cat_def.get("parent_category")  # e.g., refined -> model

            # For old-style subcategories, check if file is in parent
            if old_parent and f not in files.get(old_parent, []):
                continue

            # For new-style semantic subcategories with extension requirements
            if semantic_parent and "extensions" in cat_def:
                cat_extensions = cat_def.get("extensions", [])
                if not any(f_lower.endswith(e) for e in cat_extensions):
                    continue

            # Check patterns
            patterns = cat_def.get("patterns", [])
            excludes = cat_def.get("excludes", [])
            max_len = cat_def.get("max_basename_length")

            # Check excludes first
            excluded = False
            for exc_pattern in excludes:
                if _match_pattern(basename, exc_pattern):
                    excluded = True
                    break

            if excluded:
                continue

            # Check max basename length
            if max_len and len(os.path.splitext(basename)[0]) > max_len:
                continue

            # Check patterns
            matched = False
            if not patterns and semantic_parent:
                # No patterns = use extension matching only (already checked above)
                matched = any(f_lower.endswith(e) for e in cat_def.get("extensions", []))
            else:
                for pattern in patterns:
                    if pattern == "*":
                        # Wildcard matches all (used with excludes)
                        matched = True
                        break
                    if _match_pattern(basename, pattern):
                        matched = True
                        break

            if matched and f not in files[cat_name]:
                files[cat_name].append(f)

                # Bubble up to semantic parent category
                # e.g. refined -> model, ligand_cif -> ligand,
                #      original_data_mtz -> data_mtz, refine_map_coeffs -> map_coeffs_mtz
                parent = cat_def.get("parent_category")
                if parent and parent in files and f not in files[parent]:
                    files[parent].append(f)

                # Also add to "also_in" categories
                for also_cat in cat_def.get("also_in", []):
                    if also_cat in files and f not in files[also_cat]:
                        files[also_cat].append(f)

    # ── POST-PROCESSING ─────────────────────────────────────────────────────
    # Remove intermediate files from model and search_model categories.
    # A file may match BOTH an intermediate subcategory (e.g. *EDITED*) and
    # the unclassified_pdb catch-all (*), which bubbles up to model.  The
    # intermediate classification should always take priority — these files
    # must never be used as program inputs.
    intermediate_files = set(files.get("intermediate", []))
    if intermediate_files:
        for safe_cat in ("model", "search_model", "pdb", "unclassified_pdb"):
            if safe_cat in files:
                before = len(files[safe_cat])
                files[safe_cat] = [f for f in files[safe_cat]
                                   if f not in intermediate_files]

    return files


def _categorize_files_hardcoded(available_files, ligand_hints=None, files_local=True):
    """
    Categorize files by type and purpose.

    Args:
        available_files: List of file paths
        ligand_hints: Optional set of basenames known to be ligands

    Returns dict with keys:
        data_mtz, map_coeffs_mtz, pdb, sequence, map, ligand_cif, ligand_pdb,
        phaser_output, refined, rsr_output, with_ligand, ligand_fit, predicted,
        full_map, half_map (for cryo-EM half-maps)
    """
    files = {
        "data_mtz": [],  # Reflection data with Fobs, R-free (for refinement)
        "map_coeffs_mtz": [],  # Map coefficients with phases (for ligand fitting)
        "pdb": [],
        "sequence": [],
        "map": [],  # All map files (for backward compatibility)
        "full_map": [],  # Full cryo-EM maps
        "half_map": [],  # Half maps (usually come in pairs)
        "ligand_cif": [],
        "ligand_pdb": [],
        "phaser_output": [],
        "refined": [],  # X-ray refinement output
        "rsr_output": [],  # real_space_refine output (cryo-EM)
        "with_ligand": [],
        "ligand_fit_output": [],
        "predicted": [],
        "processed_predicted": [],
        "autobuild_output": [],
        "docked": [],  # dock_in_map output
        "intermediate_mr": [],  # Intermediate MR files (never use for refinement)
    }

    # X-ray data file extensions (PHENIX can read these)
    xray_data_extensions = ('.mtz', '.sca', '.hkl', '.sdf')

    def is_half_map(basename):
        """Detect if a map file is a half-map based on naming conventions.

        Only matches files with 'half' in the name (e.g. half_map_1.mrc,
        half1.mrc, emd_12345_half_map_2.map).

        The old regexes [_-]?[12]$ and [_-][ab]$ were removed because
        they false-positive on sequentially numbered maps like map_1.ccp4,
        map_2.ccp4 (segmented maps from resolve_cryo_em), segment_a.mrc,
        etc.  Real half-maps from EMDB, RELION, CryoSPARC, and other
        cryo-EM software all contain 'half' in their filename.
        """
        return 'half' in basename.lower()

    # Import shared MTZ classification
    classify_mtz_type, get_mtz_stage = _import_mtz_utils()

    for f in available_files:
        f_lower = f.lower()
        basename = os.path.basename(f_lower)

        # Primary type categorization
        if f_lower.endswith(xray_data_extensions):
            # Classify into data_mtz or map_coeffs_mtz
            mtz_type = classify_mtz_type(f)  # Pass full path
            files[mtz_type].append(f)
            # Also populate the specific subcategory (refine_map_coeffs, etc.)
            # so that programs requesting subcategories can find the file.
            if mtz_type == "map_coeffs_mtz":
                stage = get_mtz_stage(f, "map_coeffs_mtz")
                if stage and stage != "map_coeffs_mtz":
                    if stage not in files:
                        files[stage] = []
                    if f not in files[stage]:
                        files[stage].append(f)
        elif f_lower.endswith('.pdb'):
            files["pdb"].append(f)

            # Subcategorize PDBs by their origin/purpose
            # Track whether this PDB matched any known program-output subcategory.
            # Files that match are known outputs and must NOT be reclassified as
            # small molecules later (a small cryo-EM build can have ≤150 atoms
            # and trip _pdb_is_small_molecule).
            _is_program_output = False

            if 'phaser' in basename or basename.startswith('phaser'):
                files["phaser_output"].append(f)
                _is_program_output = True

            # MR solution files: names like mup_mr_solution.pdb
            # indicate an already-placed MR model.  Classify as
            # phaser_output so _has_placed_model recognizes them.
            # Boundary-aware: require _ before mr_solution (or at
            # start of name) so nmr_solution.pdb doesn't match.
            if ('_mr_solution' in basename or
                basename.startswith('mr_solution')):
                if f not in files["phaser_output"]:
                    files["phaser_output"].append(f)
                _is_program_output = True

            if 'refine' in basename and 'real_space' not in basename and 'rsr' not in basename:
                files["refined"].append(f)
                _is_program_output = True

            # RSR output detection - real_space_refine outputs contain 'real_space_refined'
            # e.g., model_real_space_refined_000.pdb
            if 'real_space_refined' in basename or 'rsr_' in basename or '_rsr' in basename:
                files["rsr_output"].append(f)
                _is_program_output = True

            if 'with_ligand' in basename:
                files["with_ligand"].append(f)
                _is_program_output = True
            elif '_modified' in basename and basename.endswith('.pdb'):
                # pdbtools default output: {input}_modified.pdb
                # In the agent workflow, pdbtools is used to combine
                # protein + ligand, so _modified.pdb = with_ligand
                files["with_ligand"].append(f)
                _is_program_output = True
            if 'ligand_fit' in basename or 'ligandfit' in basename:
                files["ligand_fit_output"].append(f)
                _is_program_output = True

            if 'predict' in basename or 'alphafold' in basename or 'colabfold' in basename:
                files["predicted"].append(f)
                _is_program_output = True

            if 'processed' in basename:
                files["processed_predicted"].append(f)
                _is_program_output = True

            is_autobuild = (
                ('autobuild' in basename or 'auto_build' in basename) or
                ('AutoBuild' in os.path.basename(f)) or
                ('overall_best' in basename) or
                ('build' in basename or 'built' in basename) or
                'buccaneer' in basename or
                'arp_warp' in basename or
                ('shelxe' in basename and 'trace' in basename)
            )
            if is_autobuild and 'predict' not in basename:
                files["autobuild_output"].append(f)
                _is_program_output = True

            # map_to_model output — the final positioned model from cryo-EM
            # model building.  Matched by "map_to_model" in basename.  Must
            # come AFTER the autobuild check (which catches generic 'build')
            # and BEFORE the small-molecule check (a partial build with few
            # atoms would otherwise be misclassified as a ligand).
            if 'map_to_model' in basename:
                files["autobuild_output"].append(f)
                _is_program_output = True

            if 'dock' in basename and 'map' in basename:
                files["docked"].append(f)
                _is_program_output = True
            # Also match placed_model* from dock_in_map output
            if basename.startswith('placed_model') or '_placed' in basename:
                files["docked"].append(f)
                _is_program_output = True

            # Intermediate MR files from dock_in_map - never use for refinement
            if basename.startswith('run_mr') or fnmatch.fnmatch(basename, '*mr.[0-9]*'):
                files["intermediate_mr"].append(f)

            # Ligand coordinate files — require 'lig' or 'ligand' at a word boundary
            # to avoid false positives on names like 'nsf-d2_noligand.pdb' where
            # 'ligand' appears as a suffix of a different word.
            _is_ligand_name = (
                (re.search(r'(^|[_\-\.])lig([_\-\.]|\.pdb$|$)', basename) and
                 len(basename) < 20) or
                re.search(r'(^|[_\-\.])ligand([_\-\.]|\.pdb$|$)', basename)
            )
            if _is_ligand_name:
                if not any(x in basename for x in [
                    'ligand_fit', 'ligandfit', 'with_ligand',
                    'no_ligand', 'noligand', 'sans_ligand',
                    'without_ligand', 'apo_ligand',
                ]):
                    files["ligand_pdb"].append(f)
            elif not _is_program_output:
                # Content-based detection: HETATM-only files whose names don't
                # contain 'lig' or 'ligand' (e.g. atp.pdb, gdp.pdb, hem.pdb).
                # Primary: read PDB content (works on client where files are on disk).
                # Fallback: client-side ligand_hints (works on server).
                #
                # CRITICAL: Only check files that didn't match any program-output
                # subcategory above.  A small cryo-EM map_to_model output (≤150
                # atoms) would otherwise be misclassified as a ligand.
                is_small = _pdb_is_small_molecule(f) if files_local else False
                if not is_small and ligand_hints:
                    is_small = os.path.basename(f) in ligand_hints
                if is_small:
                    # Move: pdb → ligand_pdb
                    files["pdb"].remove(f)
                    files["ligand_pdb"].append(f)

        elif f_lower.endswith(('.fa', '.fasta', '.seq', '.dat')):
            files["sequence"].append(f)
        elif f_lower.endswith(('.mrc', '.ccp4', '.map')):
            files["map"].append(f)
            if is_half_map(basename):
                files["half_map"].append(f)
            else:
                files["full_map"].append(f)
        elif f_lower.endswith('.cif'):
            if 'refine' in basename:
                files["pdb"].append(f)
                files["refined"].append(f)
            elif any(pat in basename for pat in (
                    'overall_best', 'autobuild', 'predict_and_build',
                    'autosol', 'model_cif')):
                # Whole-model mmCIF from autobuild/predict_and_build:
                # geometry restraints for the full protein, NOT a ligand.
                # Putting these in ligand_cif causes them to be injected
                # positionally into phenix.refine → "wrong number of models".
                files["pdb"].append(f)
            else:
                files["ligand_cif"].append(f)

    # ── Post-categorization: MTZ array detection (Fix I1) ──
    # Detect MTZ files with multiple observation arrays.
    # When found, record the available labels so the command
    # builder can inject explicit labels= and avoid the
    # "ambiguous data" crash in xtriage/autosol/phaser.
    files["mtz_array_info"] = _detect_mtz_arrays(
        files.get("data_mtz", []))

    return files


def _mtz_has_phase_columns(filepath):
    """Return True if the MTZ file contains experimental phase columns.

    Checks for PHIB (type 'P') or Hendrickson-Lattman coefficients
    (type 'A') using iotbx.  Falls back to label-name heuristics if
    iotbx is unavailable (e.g. in unit-test environments).

    Used by the phased_data_mtz categorisation to prevent a correctly-
    named but wrong-content file from replacing one crash (missing file)
    with another (missing labels inside autobuild).

    Never raises.
    """
    if not filepath or not os.path.isfile(filepath):
        return False
    try:
        # Primary path: iotbx miller array column types
        try:
            from iotbx.reflection_file_reader import any_reflection_file
        except ImportError:
            from iotbx import reflection_file_reader
            any_reflection_file = reflection_file_reader.any_reflection_file

        reader = any_reflection_file(filepath)
        if reader is None:
            return False
        for ma in reader.as_miller_arrays():
            # type_string() returns one char: P=phase, A=HL coefficients
            type_str = ""
            try:
                type_str = ma.info().type_string()
            except Exception:
                pass
            if type_str in ('P', 'A'):
                return True
        return False

    except Exception:
        pass

    # Fallback: direct mtz object column types
    try:
        from iotbx import mtz as _iotbx_mtz
        obj = _iotbx_mtz.object(filepath)
        for col in obj.columns():
            if col.type() in ('P', 'A'):
                return True
        return False
    except Exception:
        pass

    # Last resort: raw label-name heuristic on ALL mtz columns.
    # _read_mtz_array_labels() returns only F/I arrays, so it would
    # never find PHIB/FOM (which are phase columns, not obs arrays).
    # Instead, iterate every column via iotbx.mtz.object directly.
    _PHASE_LABEL_SUBSTRINGS = (
        'PHIB', 'PHIM', 'FOM', 'FOMM',
        'HLA', 'HLB', 'HLC', 'HLD',
    )
    try:
        from iotbx import mtz as _iotbx_mtz2
        obj2 = _iotbx_mtz2.object(filepath)
        for col in obj2.columns():
            label_up = col.label().upper()
            if any(p in label_up for p in _PHASE_LABEL_SUBSTRINGS):
                return True
    except Exception:
        pass

    return False


# ── Module-level cache for phase-column detection (v115.08) ────────
# Persists across cycles within the same process.  Same pattern as
# _load_done_tracking_configs() which caches YAML parsing at module
# level.  Key: (absolute_path, mtime) — absolute path prevents
# cross-directory collisions in long-running processes; mtime ensures
# re-check when programs overwrite files (e.g. phenix.refine produces
# a new MTZ each cycle).
_PHASE_COLUMN_CACHE = {}


def _has_phase_columns_cached(filepath):
    """Cached wrapper for _mtz_has_phase_columns().

    Cache key: (os.path.abspath(filepath), mtime).  Returns cached
    result on hit; calls _mtz_has_phase_columns() on miss and stores
    the result.

    Never raises.
    """
    try:
        abs_path = os.path.abspath(filepath)
        mtime = os.path.getmtime(abs_path)
    except (OSError, TypeError):
        return False
    key = (abs_path, mtime)
    if key in _PHASE_COLUMN_CACHE:
        return _PHASE_COLUMN_CACHE[key]
    result = _mtz_has_phase_columns(abs_path)
    _PHASE_COLUMN_CACHE[key] = result
    return result


# ── Tier 2: ASCII header heuristic for .hkl/.sca files (v115.08) ──
# When iotbx cannot resolve column types (e.g. .hkl files without a
# .def file), fall back to scanning the first 100 lines for known
# phase column label tokens.

_PHASE_LABEL_TOKENS = frozenset({
    'PHIB', 'PHIM', 'FOM', 'FOMM',
    'HLA', 'HLB', 'HLC', 'HLD',
})

_COMMENT_CHARS = ('#', '!', '*', ';')


def _ascii_phase_heuristic(filepath):
    """Check ASCII reflection file header for phase column labels.

    Reads first 100 lines.  Skips blank and comment lines (starting
    with #, !, *, ;).  Tokenizes remaining lines and checks for
    known phase column labels as exact whitespace-separated tokens.
    This avoids false-positives from comments containing phase-related
    words.

    Never raises.
    """
    if not filepath or not os.path.isfile(filepath):
        return False
    try:
        with open(filepath, 'r') as fh:
            for i, line in enumerate(fh):
                if i >= 100:
                    break
                stripped = line.strip()
                if not stripped:
                    continue
                # Skip comment lines
                if stripped[0] in _COMMENT_CHARS:
                    continue
                tokens = stripped.upper().split()
                if _PHASE_LABEL_TOKENS & set(tokens):
                    return True
    except Exception:
        pass  # Binary file or encoding issue — not ASCII phases
    return False


def _resolve_phased_promotions(files):
    """Promote data_mtz files with phase columns to phased_data_mtz.

    Replaces the old marker-gated promotion loop (v115.07 Fix I3).
    Now content-based: every file in data_mtz is checked for actual
    phase columns, regardless of filename.

    Two-tier detection:
      Tier 1: iotbx column-type check (P or A types) — cached
      Tier 2: ASCII header heuristic (for .hkl without iotbx headers)

    Atomic: each promoted file is added to phased_data_mtz AND removed
    from data_mtz in a single pass.  No file exists in both categories
    after this function returns.

    Never raises.
    """
    if "phased_data_mtz" not in files:
        files["phased_data_mtz"] = []

    to_promote = []
    for f in files.get("data_mtz", []):
        # Tier 1: iotbx content check (cached)
        if _has_phase_columns_cached(f):
            to_promote.append(f)
            continue
        # Tier 2: ASCII header heuristic (for .hkl/.sca files)
        if _ascii_phase_heuristic(f):
            to_promote.append(f)

    # Atomic: promote + remove in one pass.
    # IMPORTANT: only remove from data_mtz when OTHER data files
    # remain.  A phased file (e.g. nsf-d2_phases.hkl) contains
    # BOTH Fobs AND phases — it IS valid data for refinement.
    # Programs like phenix.refine/phaser/autosol do NOT check
    # phased_data_mtz (only data_mtz/original_data_mtz).  If we
    # remove the only data file, those programs can't find data.
    # When a separate scale/data file exists (like rab3a_scale.hkl),
    # removal is correct — it prevents ambiguous file selection.
    _non_promoted_data = [
        f for f in files.get("data_mtz", [])
        if f not in to_promote
    ]
    _can_remove = len(_non_promoted_data) > 0

    for f in to_promote:
        if f not in files["phased_data_mtz"]:
            files["phased_data_mtz"].append(f)
        if _can_remove and f in files["data_mtz"]:
            files["data_mtz"].remove(f)
            try:
                logger.debug(
                    "CATEGORIZE: Promoted %s to phased_data_mtz "
                    "(removed from data_mtz; other data files "
                    "remain)", os.path.basename(f))
            except Exception:
                pass  # Logging unavailable in test environments
        elif not _can_remove and f in files["data_mtz"]:
            try:
                logger.debug(
                    "CATEGORIZE: Promoted %s to phased_data_mtz "
                    "(kept in data_mtz — only data file, needed "
                    "by refine/phaser/autosol)",
                    os.path.basename(f))
            except Exception:
                pass  # Logging unavailable in test environments

    # ── Category exclusivity enforcement ──────────────────────
    # When the YAML categorizer placed a file in BOTH data_mtz
    # (by extension) and phased_data_mtz (by *phases* pattern),
    # and the content check didn't fire (e.g. iotbx can't parse
    # the .hkl), enforce exclusivity — but ONLY when other data
    # files remain.  Same rationale as above: programs like
    # phenix.refine only look in data_mtz.
    _phased_set = set(files.get("phased_data_mtz", []))
    if _phased_set:
        _data_after_removal = [
            f for f in files.get("data_mtz", [])
            if f not in _phased_set
        ]
        if len(_data_after_removal) > 0:
            _before = len(files.get("data_mtz", []))
            files["data_mtz"] = _data_after_removal
            _removed = _before - len(files["data_mtz"])
            if _removed > 0:
                try:
                    logger.debug(
                        "CATEGORIZE: Exclusivity enforcement "
                        "removed %d file(s) from data_mtz "
                        "(other data files remain)", _removed)
                except Exception:
                    pass  # Logging unavailable in test envs
        # else: all data_mtz files are phased — keep them in
        # both categories so refine/phaser can find data.

    return files


def _detect_mtz_arrays(data_mtz_files):
    """Detect MTZ files with multiple observation arrays (Fix I1).

    When an MTZ has both merged (Iobs) and anomalous (I(+)/I(-))
    arrays, programs like xtriage crash without explicit labels.

    Returns:
        dict: {filepath: {"merged": label_str,
                          "anomalous": label_str,
                          "preferred": label_str}}
        Empty dict if no multi-array files found or iotbx
        unavailable.

    Never raises.
    """
    if not data_mtz_files:
        return {}

    result = {}
    for filepath in data_mtz_files:
        if not filepath.lower().endswith('.mtz'):
            continue
        arrays = _read_mtz_array_labels(filepath)
        if arrays and len(arrays) > 1:
            info = _classify_mtz_arrays(arrays)
            if info.get("merged") and info.get("anomalous"):
                result[filepath] = info
    return result


def _read_mtz_array_labels(filepath):
    """Read observation array labels from an MTZ file.

    Returns list of (label_string, is_anomalous) tuples,
    or empty list if iotbx unavailable or file unreadable.
    """
    try:
        try:
            from iotbx.reflection_file_reader import (
                any_reflection_file)
        except ImportError:
            from iotbx import reflection_file_reader
            any_reflection_file = (
                reflection_file_reader
                .any_reflection_file)

        reader = any_reflection_file(filepath)
        if reader is None:
            return []
        arrays = reader.as_miller_arrays()
        obs_arrays = []
        for ma in arrays:
            if ma.is_xray_intensity_array():
                label = ma.info().label_string()
                is_anom = ma.anomalous_flag()
                obs_arrays.append((label, is_anom))
            elif ma.is_xray_amplitude_array():
                label = ma.info().label_string()
                is_anom = ma.anomalous_flag()
                obs_arrays.append((label, is_anom))
        return obs_arrays
    except Exception:
        return []


def _classify_mtz_arrays(arrays):
    """Classify MTZ arrays into merged/anomalous with
    a preferred default.

    Args:
        arrays: list of (label_string, is_anomalous)

    Returns:
        dict with "merged", "anomalous", "preferred" keys.
    """
    merged = None
    anomalous = None
    for label, is_anom in arrays:
        if is_anom and anomalous is None:
            anomalous = label
        elif not is_anom and merged is None:
            merged = label

    # Ranking rule: default to merged (safer for MR/refine).
    # SAD/MAD workflows will override to anomalous based on
    # experiment context in the command builder.
    preferred = merged or anomalous
    return {
        "merged": merged,
        "anomalous": anomalous,
        "preferred": preferred,
    }


def _detect_experiment_type(files, history_info=None):
    """
    Determine if this is X-ray crystallography or Cryo-EM.

    Logic:
    - If mtriage has been run → cryo-EM (definitive)
    - Has map (full or half) but no MTZ → cryo-EM
    - Has MTZ → X-ray (even if map also present)
    - Neither → unknown (default to X-ray)
    """
    if history_info and history_info.get("mtriage_done"):
        return "cryoem"

    if history_info and history_info.get("rsr_done"):
        return "cryoem"

    has_data_mtz = bool(files.get("data_mtz")) or bool(files.get("map_coeffs_mtz"))
    has_map = bool(files.get("map")) or bool(files.get("full_map")) or bool(files.get("half_map"))

    if has_map and not has_data_mtz:
        return "cryoem"
    else:
        return "xray"


# =============================================================================
# HISTORY ANALYSIS
# =============================================================================

def _is_failed_result(result):
    """
    Check if a result string indicates failure.

    Uses specific Phenix terminal-error phrases to avoid false positives.
    Bare "ERROR" and ": ERROR" are intentionally excluded because Phenix logs
    contain these in non-fatal contexts (parameter descriptions, help text,
    "No ERROR detected", "Expected errors: 0", etc.)

    Priority order per J2 audit spec:
      1. Exit-code check (done at call sites before calling this function)
      2. Output-file check (done by _clear_zombie_done_flags / J5)
      3. Log-text check (this function) — uses specific terminal phrases only

    Args:
        result: Result string from history entry

    Returns:
        bool: True if result indicates a definitive failure
    """
    if not result:
        return False

    result_upper = result.upper()

    # Simple substring patterns — each chosen to be specific enough to avoid
    # matching non-fatal Phenix log content.
    # 'ERROR:' in any form is intentionally excluded: "No error: all checks passed"
    # uppercases to "NO ERROR: ALL CHECKS PASSED" which would match, and
    # Phenix-level errors are always reported via "Sorry:", "*** Error:", or exit code.
    simple_patterns = [
        'FAILED',           # Common explicit failure indicator
        'SORRY:',           # Phenix-specific error prefix (e.g. "Sorry: bad input")
        'SORRY ',           # Phenix error prefix followed by text
        '*** ERROR',        # Phenix fatal error banner (*** Error: ...)
        'FATAL:',           # Fatal error prefix
        'TRACEBACK',        # Python exception traceback header
        'EXCEPTION',        # Python exception indicator
    ]

    return any(p in result_upper for p in simple_patterns)


# Cache for done tracking config loaded from YAML
_DONE_TRACKING_CACHE = None

# Allowed count_field values — rejects typos at load time
ALLOWED_COUNT_FIELDS = {"refine_count", "rsr_count", "phaser_count"}


def _load_done_tracking_configs():
    """Load done_tracking configuration from programs.yaml.

    Returns list of dicts for all programs with history_detection, each with:
      - flag: done flag name (e.g., 'dock_done')
      - strategy: 'set_flag' (default), 'run_once', or 'count'
      - count_field: counter name for strategy='count' (e.g., 'refine_count')
      - markers: list of strings to match in combined text (OR logic,
                 substring matching)
      - exclude_markers: list of strings that reject a match (checked FIRST)
      - alt_markers: optional list of alternate marker strings
      - alt_requires: optional list of strings that must ALL be present
                      alongside alt_markers (AND logic)
      - success_flag: optional additional flag set on success

    Validates count_field against ALLOWED_COUNT_FIELDS at load time to
    prevent typos from silently creating garbage attributes.
    """
    global _DONE_TRACKING_CACHE
    if _DONE_TRACKING_CACHE is not None:
        return _DONE_TRACKING_CACHE

    configs = []
    try:
        try:
            from libtbx.langchain.knowledge.yaml_loader import load_programs
        except ImportError:
            from knowledge.yaml_loader import load_programs
        programs = load_programs()
        for name, defn in programs.items():
            if not isinstance(defn, dict):
                continue
            tracking = defn.get("done_tracking", {})
            detection = tracking.get("history_detection")
            if not detection or not isinstance(detection, dict):
                continue

            strategy = tracking.get("strategy", "set_flag")
            count_field = tracking.get("count_field")

            # Validate count_field at load time
            if strategy == "count":
                if not count_field:
                    print("  [workflow_state] Warning: %s has strategy='count' "
                          "but no count_field" % name)
                    continue
                if count_field not in ALLOWED_COUNT_FIELDS:
                    raise ValueError(
                        "Unknown count_field %r in done_tracking for %s. "
                        "Allowed: %s" % (count_field, name, ALLOWED_COUNT_FIELDS))

            configs.append({
                "program": name,
                "flag": tracking.get("flag"),
                "strategy": strategy,
                "count_field": count_field,
                "markers": detection.get("markers", []),
                "exclude_markers": detection.get("exclude_markers", []),
                "alt_markers": detection.get("alt_markers", []),
                "alt_requires": detection.get("alt_requires", []),
                "success_flag": detection.get("success_flag"),
            })
    except Exception as e:
        print("  [workflow_state] Warning: could not load done_tracking "
              "from programs.yaml: %s" % str(e))

    _DONE_TRACKING_CACHE = configs
    return configs


def _set_done_flags(info, combined, result):
    """Set done flags from YAML history_detection for all strategies.

    Handles: set_flag, run_once, and count strategies.
    The run_once filtering itself happens in program_registration;
    here we just set the flag and increment counts.

    Args:
        info: The info dict being built by _analyze_history
        combined: Lowercase string of program + command
        result: Result string from history entry
    """
    if _is_failed_result(result):
        return  # All strategies require success

    configs = _load_done_tracking_configs()

    for config in configs:
        flag = config["flag"]
        if not flag:
            continue

        # Exclude markers take precedence — checked FIRST
        if any(m in combined for m in config["exclude_markers"]):
            continue

        # Check primary markers (OR logic, substring matching)
        matched = any(m in combined for m in config["markers"])

        # Check alt_markers with alt_requires (AND logic)
        if not matched and config["alt_markers"] and config["alt_requires"]:
            if (any(m in combined for m in config["alt_markers"]) and
                    all(r in combined for r in config["alt_requires"])):
                matched = True

        if matched:
            info[flag] = True

            # Count strategy: increment count field
            if config["strategy"] == "count" and config["count_field"]:
                info[config["count_field"]] = info.get(config["count_field"], 0) + 1

            # Optional success flag
            if config["success_flag"]:
                info[config["success_flag"]] = True


def _analyze_history(history):
    """
    Extract information about what has been done from history.

    Returns dict with program completion flags and counts.

    Note: Simple done flags for programs with run_once: true are auto-generated
    from programs.yaml via program_registration. Complex flags (counts, success
    conditions) are still handled manually below.
    """
    # =========================================================================
    # Initialize flags from YAML done_tracking configs
    # =========================================================================
    info = {
        "programs_run": set(),
        # predict_and_build flags — no history_detection (Python-only cascade)
        "predict_done": False,
        "predict_full_done": False,
        # Post-ligandfit refinement tracking
        "needs_post_ligandfit_refine": False,
        # Metrics extracted from history
        "last_program": None,
        "last_r_free": None,
        "last_map_cc": None,
        "last_clashscore": None,
        "last_tfz": None,
        "resolution": None,
        "anomalous_resolution": None,
        "anomalous_measurability": None,
        "has_anomalous": False,
        "strong_anomalous": False,
        "has_twinning": False,
        "twin_law": None,
        "twin_fraction": None,
        "has_ncs": False,  # NCS detected in data

        # Placement probe results (Tier 3)
        "placement_probed": False,
        "placement_probe_result": None,
    }

    # Initialize done flags and count fields from YAML (single source of truth)
    configs = _load_done_tracking_configs()
    for config in configs:
        if config["flag"]:
            info[config["flag"]] = False
        if config["strategy"] == "count" and config["count_field"]:
            info[config["count_field"]] = 0
        if config.get("success_flag"):
            info[config["success_flag"]] = False

    if not history:
        return info

    # =========================================================================
    # Process history for flags and metrics
    # =========================================================================
    for entry in history:
        prog = ""
        cmd = ""

        if isinstance(entry, str):
            prog = entry.lower()
            cmd = entry.lower()
        elif isinstance(entry, dict):
            prog = (entry.get("program") or "").lower()
            cmd = (entry.get("command") or "").lower()
            # Extract metrics
            # NOTE: history from session has 'analysis' key, but after transport has 'metrics'
            analysis = entry.get("analysis", entry.get("metrics", {}))
            if isinstance(analysis, dict):
                # Coerce numeric values — JSON round-tripping can
                # turn floats into strings (e.g. 0.385 → "0.385").
                def _sf(v):
                    if v is None:
                        return None
                    try:
                        return float(v)
                    except (ValueError, TypeError):
                        return None
                if analysis.get("r_free"):
                    info["last_r_free"] = _sf(analysis["r_free"])
                if analysis.get("map_cc"):
                    info["last_map_cc"] = _sf(analysis["map_cc"])
                if analysis.get("clashscore"):
                    info["last_clashscore"] = _sf(analysis["clashscore"])
                if analysis.get("tfz"):
                    info["last_tfz"] = _sf(analysis["tfz"])
                if analysis.get("resolution"):
                    info["resolution"] = _sf(analysis["resolution"])
                if analysis.get("anomalous_resolution"):
                    _ar = _sf(analysis["anomalous_resolution"])
                    info["anomalous_resolution"] = _ar
                    # anomalous_resolution is the estimated d-spacing limit of
                    # useful anomalous signal from xtriage.  A value > 6 Å means
                    # signal is only present at very low angles (essentially noise
                    # from Friedel pair differences) and is NOT useful for phasing.
                    # Example: anomalous_resolution=9.8 with measurability=0.032
                    # is negligible and must NOT gate autosol availability.
                    if _ar is not None and _ar < 6.0:
                        info["has_anomalous"] = True
                # Independent 'if' (not elif): has_anomalous from analysis
                # must be checked even when anomalous_resolution is present
                # but >= 6.0.  The old 'elif' silently skipped this when
                # anomalous_resolution was set.  (v115.03 fix)
                if analysis.get("has_anomalous"):
                    info["has_anomalous"] = analysis["has_anomalous"]
                # Store anomalous measurability for decision making
                if analysis.get("anomalous_measurability") is not None:
                    info["anomalous_measurability"] = analysis["anomalous_measurability"]
                    # Strong anomalous signal if measurability > 0.10
                    if analysis["anomalous_measurability"] > 0.10:
                        info["has_anomalous"] = True
                        info["strong_anomalous"] = True
                    elif analysis["anomalous_measurability"] >= 0.06:
                        info["has_anomalous"] = True  # weak signal (v115.03)
                    elif analysis["anomalous_measurability"] < 0.06:
                        # Negligible signal — clear unless protected
                        if info.get("strong_anomalous"):
                            pass  # Don't clear strong signal
                        elif analysis.get("has_anomalous"):
                            pass  # Trust explicit flag (v115.03)
                        else:
                            info["has_anomalous"] = False
                # Twinning threshold (0.20) from workflows.yaml shared section
                if analysis.get("twin_law") and analysis.get("twin_fraction"):
                    twin_frac = analysis["twin_fraction"]
                    if twin_frac > 0.20 and not analysis.get("no_twinning_suspected"):
                        info["has_twinning"] = True
                        info["twin_law"] = analysis["twin_law"]
                        info["twin_fraction"] = twin_frac
                # NCS detection (from map_symmetry or similar)
                if analysis.get("ncs_found") or analysis.get("has_ncs"):
                    info["has_ncs"] = True
        else:
            continue

        combined = prog + " " + cmd

        info["programs_run"].add(prog)
        info["last_program"] = prog

        # =====================================================================
        # Program done flags — YAML-driven detection for all strategies
        # =====================================================================
        result = entry.get("result", "") if isinstance(entry, dict) else ""

        # Handles: set_flag, run_once, and count strategies.
        # Covers all programs with history_detection in programs.yaml.
        _set_done_flags(info, combined, result)

        # The ONE remaining Python-only case: predict_and_build cascade.
        # When predict runs fully, it sets predict_full_done.
        # We intentionally do NOT set refine_done here: predict_and_build
        # performs internal refinement passes, but those are not the same as
        # a standalone phenix.refine run on the final output model.  Setting
        # refine_done=True here caused the LOOP WARNING (refine_done=True but
        # phenix.refine still in valid_programs) and also triggered the file
        # selector to pick intermediate multi-model PDBs from PredictAndBuild
        # subdirectories instead of the correct overall_best output.
        if "predict_and_build" in combined:
            if not _is_failed_result(result):
                info["predict_done"] = True
                if "stop_after_predict=true" not in combined and "stop_after_predict=True" not in combined:
                    info["predict_full_done"] = True
                    # refine_count is NOT incremented here; the agent should
                    # run a separate phenix.refine step on the p&b output.

        # Track whether refinement is needed after ligandfit.
        # After ligandfit adds a ligand, the complex always needs re-refinement.
        # This flag is True when ligandfit succeeded but no refine happened after.
        if not _is_failed_result(result):
            if "ligandfit" in combined:
                info["needs_post_ligandfit_refine"] = True
            elif "refine" in combined and "real_space" not in combined:
                info["needs_post_ligandfit_refine"] = False

    # =========================================================================
    # Placement probe detection (Tier 3 post-processing)
    # =========================================================================
    # Identify whether model_vs_data or map_correlations ran as a *placement
    # probe* (i.e. before any refinement or docking) vs. as normal validation.
    #
    # Strategy: iterate history a second time, tracking whether we have seen
    # any refine/dock cycle yet.  The FIRST occurrence of model_vs_data or
    # map_correlations that appears BEFORE the first refine/dock is the probe.
    #
    # This requires no schema change — phase is inferred from position.
    if history:
        _seen_refine_or_dock = False
        for _entry in history:
            if not isinstance(_entry, dict):
                continue
            _eprog  = (_entry.get("program") or "").lower()
            _ecmd   = (_entry.get("command") or "").lower()
            _ecomb  = _eprog + " " + _ecmd
            _result = _entry.get("result", "")

            if _is_failed_result(_result):
                # ── S2L: probe failures that carry placement information ──────
                # A failed map_correlations run before any refine/dock is the
                # cryo-EM placement probe.  A failed model_vs_data run before
                # any refine is the X-ray placement probe.  Some crashes are
                # themselves definitive answers — extract the signal before
                # discarding the entry.
                #
                # "entirely outside map": model coordinates are completely
                #   outside the map box → model is definitively not placed →
                #   needs_dock.
                #
                # "crystal symmetry mismatch" from model_vs_data:
                #   The probe couldn't even compare — INCONCLUSIVE.  We do NOT
                #   route to MR, because phenix.refine is often more permissive
                #   (it can handle small cell differences) and may succeed where
                #   model_vs_data refuses.  Mark as probed/inconclusive so the
                #   probe is not retried, then fall through to normal refine.
                #
                # Any other failure: mark placement_probed=True with result=None
                #   (inconclusive) so the probe never runs a third time.
                #   Routing: placement_uncertain becomes False, falls through to
                #   obtain_model → predict/dock fallback.
                _is_xray_probe = ("model_vs_data" in _ecomb and not _seen_refine_or_dock)
                _is_cryoem_probe = ("map_correlations" in _ecomb and not _seen_refine_or_dock)
                if _is_cryoem_probe or _is_xray_probe:
                    _rl = (_result or "").lower()
                    _outside_signals = [
                        "entirely outside map",
                        "outside map",
                        "model is outside",
                        "model entirely outside",
                        "stopping as model",
                    ]
                    if any(s in _rl for s in _outside_signals):
                        # Hard evidence: model is not in the map at all
                        info["placement_probed"] = True
                        info["placement_probe_result"] = "needs_dock"
                    elif not info.get("placement_probed"):
                        # Unknown failure (including crystal_symmetry_mismatch
                        # from model_vs_data) — prevent infinite probe retry
                        # but treat as inconclusive so refine can still run.
                        info["placement_probed"] = True
                        # Leave placement_probe_result as None (inconclusive)
                    continue   # Ignore failed cycles for all other probe detection

                # ── Failed RSR/refine with symmetry mismatch ─────────
                # "Symmetry and/or box dimensions mismatch" from
                # real_space_refine is definitive evidence that the
                # model has NOT been docked into this map.  The same
                # error from phenix.refine (X-ray) means the model
                # cell doesn't match the data cell → needs MR.
                # Set placement_probed so the next cycle routes to
                # dock_in_map (cryo-EM) or phaser (X-ray).
                _rl = (_result or "").lower()
                _sym_mismatch_signals = [
                    "symmetry and/or box",
                    "dimensions mismatch",
                    "unit cell mismatch",
                    "unit_cell_mismatch",
                    "crystal symmetry mismatch",
                ]
                _is_rsr = "real_space_refine" in _ecomb
                _is_refine = (
                    "refine" in _ecomb
                    and "real_space" not in _ecomb
                    and "model_vs_data" not in _ecomb)
                if (_is_rsr or _is_refine) and any(
                    s in _rl for s in _sym_mismatch_signals
                ):
                    if not info.get("placement_probed"):
                        info["placement_probed"] = True
                        info["placement_probe_result"] = (
                            "needs_dock" if _is_rsr
                            else "needs_mr")
                continue   # Ignore failed cycles for remaining detection

            # Once refinement or docking has run, anything after is validation
            if ("refine" in _ecomb and "real_space" not in _ecomb and
                    "model_vs_data" not in _ecomb):
                _seen_refine_or_dock = True
            if "dock_in_map" in _ecomb or "real_space_refine" in _ecomb:
                _seen_refine_or_dock = True

            # model_vs_data before any refine → X-ray placement probe
            if "model_vs_data" in _ecomb and not _seen_refine_or_dock:
                _analysis = _entry.get("analysis", _entry.get("metrics", {}))

                # ── Try R-free first (primary metric) ────────────────────────
                _rfree = None
                if isinstance(_analysis, dict):
                    _rfree = _analysis.get("r_free")
                # Also try to parse from result text if not in analysis
                if _rfree is None:
                    _m = re.search(r'r_free\s*[=:]\s*([0-9.]+)', _result)
                    if _m:
                        try:
                            _rfree = float(_m.group(1))
                        except ValueError:
                            pass
                if _rfree is not None:
                    try:
                        _rfree = float(_rfree)
                        info["placement_probed"] = True
                        info["placement_probe_result"] = (
                            "placed" if _rfree < 0.50 else "needs_mr"
                        )
                    except (ValueError, TypeError):
                        pass

                # ── R-work fallback when R-free is absent or "None" ──────────
                # model_vs_data sometimes outputs "r_free: None" when the
                # reflection file has no free-set flags.  R-work > 0.45 is
                # equally unambiguous evidence that the model is not placed
                # (a correctly placed model has R-work ≈ 0.20–0.35).
                if not info.get("placement_probed"):
                    _rwork = None
                    if isinstance(_analysis, dict):
                        _rwork = _analysis.get("r_work")
                    if _rwork is None:
                        # Match both 'r_work: 0.58' and 'R Work: 0.58' (display form)
                        _mw = re.search(r'r.?work\s*[=:]\s*([0-9.]+)', _result,
                                        re.IGNORECASE)
                        if _mw:
                            try:
                                _rwork = float(_mw.group(1))
                            except ValueError:
                                pass
                    if _rwork is not None:
                        try:
                            _rwork = float(_rwork)
                            info["placement_probed"] = True
                            # R-work < 0.45 → model is placed (conservative threshold)
                            # R-work ≥ 0.45 → model is not placed → needs MR
                            info["placement_probe_result"] = (
                                "placed" if _rwork < 0.45 else "needs_mr"
                            )
                        except (ValueError, TypeError):
                            pass

            # map_correlations before any refine/dock → cryo-EM placement probe
            if "map_correlations" in _ecomb and not _seen_refine_or_dock:
                _analysis = _entry.get("analysis", _entry.get("metrics", {}))
                _cc = None
                if isinstance(_analysis, dict):
                    _cc = (_analysis.get("cc_mask") or _analysis.get("map_cc")
                           or _analysis.get("cc_volume"))
                if _cc is not None:
                    try:
                        _cc = float(_cc)
                        info["placement_probed"] = True
                        info["placement_probe_result"] = (
                            "placed" if _cc > 0.15 else "needs_dock"
                        )
                    except (ValueError, TypeError):
                        pass

    # ── Post-probe correction: unset validation_done if it was set by a ─────
    # placement probe rather than actual post-refinement validation.
    # model_vs_data has done_tracking.flag = "validation_done" in YAML,
    # so _set_done_flags() always sets it.  But if model_vs_data ran BEFORE
    # any refinement (i.e. as a placement probe), it's not real validation.
    # Leaving it set causes the workflow to think validation is complete and
    # skip to "complete" phase prematurely.
    if info.get("placement_probed") and info.get("validation_done"):
        # Check if any validation program ran AFTER refinement
        _validation_after_refine = False
        _seen_refine = False
        for _entry in (history or []):
            if not isinstance(_entry, dict):
                continue
            _ep = (_entry.get("program") or "").lower()
            if "refine" in _ep and "model_vs_data" not in _ep:
                _seen_refine = True
            if _seen_refine and _ep in ("phenix.molprobity", "phenix.model_vs_data",
                                         "phenix.map_correlations", "phenix.validation_cryoem"):
                _validation_after_refine = True
                break
        if not _validation_after_refine:
            info["validation_done"] = False

    return info


# =============================================================================
# J5: ZOMBIE STATE DETECTION
# =============================================================================
# When a session restarts, history may record a program as "done" but its key
# output file may be missing (crashed write, user deleted it, partial run).
# Without this check the phase detector sees done_flag=True and skips the
# program, but the dependent file flags (has_full_map, has_placed_model, etc.)
# are False because _categorize_files found nothing.  The workflow is stuck.
#
# Strategy (from audit spec, Trap 2):
#   1. Check each done_flag against its expected output pattern in available_files.
#   2. If the history says done but NO matching file exists → clear the done_flag
#      AND the associated file-based flag (in-memory only, do NOT rewrite history).
#   3. Log a clear diagnostic so the user can see what happened.
#   4. The program becomes re-eligible on the next cycle.

import re as _re

# Map: done_flag → (file_pattern, file_flag_to_clear)
# Zombie check table entries:
#   (done_flag, output_pattern, file_flag_to_clear, accept_any_pdb)
#
# done_flag:         key in history_info to check (str)
# output_pattern:    regex matched against basenames of available_files
# file_flag_to_clear: context flag also cleared when zombie detected (or None)
# accept_any_pdb:    if True, the presence of ANY .pdb file counts as sufficient
#                    evidence that a model exists — the pattern is tried first but
#                    any .pdb is an acceptable fallback.  This prevents false
#                    zombie detections when tests (or users) name output files
#                    generically rather than using phenix's default naming.
#
# Note: resolve_cryo_em and dock_in_map use accept_any_pdb=False because they
# produce specifically-named map/model files that the workflow depends on by name.
# Model-refining programs (predict, refine, rsr) use accept_any_pdb=True because
# any PDB on disk is plausible evidence of their output.
_ZOMBIE_CHECK_TABLE = [
    # resolve_cryo_em produces denmod_map.ccp4 (or similar .ccp4/.mrc output)
    ("resolve_cryo_em_done", _re.compile(r"denmod.*\.(ccp4|mrc)$", _re.IGNORECASE),
     "has_full_map", False),
    # predict_and_build full run produces *_overall_best.pdb;
    # accept_any_pdb=True: generic names like "predicted_model.pdb" are valid evidence
    ("predict_full_done",    _re.compile(r"overall_best.*\.pdb$|.*_best.*\.pdb$", _re.IGNORECASE),
     "has_placed_model", True),
    # dock_in_map produces *_docked.pdb, *dock*map*.pdb, placed_model*.pdb, or *_placed*.pdb
    # Pattern mirrors the "docked" file category in file_categories.yaml
    ("dock_done",            _re.compile(r"docked.*\.pdb$|.*_docked.*\.pdb$|dock.*map.*\.pdb$|placed_model.*\.pdb$|.*_placed.*\.pdb$", _re.IGNORECASE),
     "has_placed_model", True),
    # phenix.refine produces *_refine_001.pdb (or higher numbered);
    # accept_any_pdb=True: generic names like "refined_model.pdb" are valid evidence
    ("refine_done",          _re.compile(r"refine_\d+\.pdb$", _re.IGNORECASE),
     None, True),
    # real_space_refine produces *_real_space_refined*.pdb;
    # accept_any_pdb=True: generic names like "refined_model.pdb" are valid evidence
    ("rsr_done",             _re.compile(r"real_space_refin.*\.pdb$", _re.IGNORECASE),
     None, True),
    # map_sharpening: output filenames vary: *_sharpened*.ccp4, *auto_sharpen*.ccp4,
    # *sharpen_A*.ccp4, etc.  Match "sharpen" (covers "sharpened" too).
    # The 'Sorry:' failure mode IS caught by _is_failed_result (has PHIB/FOM
    # pattern matches SORRY:).  This entry is a backstop for any graceful
    # failure mode that does not produce SORRY:/FAILED in the result string
    # (P5 fix, log-confirmed 2026-03-09).
    ("map_sharpening_done",  _re.compile(
        r"sharpen.*\.(ccp4|mrc)$|.*sharpen.*\.(ccp4|mrc)$",
        _re.IGNORECASE),
     None, False),
]


def _clear_zombie_done_flags(history_info, available_files, log_func=None):
    """
    Detect and resolve zombie states: done_flag=True but output file missing.

    Modifies history_info in-place.  Does NOT rewrite persisted history.
    Returns list of diagnostic messages (empty if no zombies found).

    Args:
        history_info: Dict produced by _analyze_history (modified in-place)
        available_files: List of file paths currently on disk (already filtered
                         for zero-byte files by _validate_file())
        log_func: Optional callable(msg) for logging
    """
    diagnostics = []

    if not available_files:
        available_files = []

    # Build a set of basenames for fast matching.
    # Filter out any non-string elements (defensive against malformed available_files).
    basenames = {os.path.basename(f).lower() for f in available_files
                 if isinstance(f, (str, bytes, os.PathLike))}
    any_pdb_on_disk = any(bn.endswith(".pdb") for bn in basenames)

    for done_flag, pattern, file_flag, accept_any_pdb in _ZOMBIE_CHECK_TABLE:
        if not history_info.get(done_flag):
            continue  # Not marked done → no zombie possible

        # Check if any available file matches the expected output pattern
        output_found = any(pattern.search(bn) for bn in basenames)

        # Fallback: for model-producing programs, any .pdb on disk is acceptable
        # evidence that something was produced. This prevents false zombie detections
        # when files are named generically (e.g. "refined_model.pdb" instead of
        # "model_real_space_refined_000.pdb").
        if not output_found and accept_any_pdb and any_pdb_on_disk:
            output_found = True

        if not output_found:
            msg = ("PERCEIVE: %s=True but output file (%s) not found "
                   "— clearing to allow re-run" % (done_flag, pattern.pattern))
            diagnostics.append(msg)
            if log_func:
                log_func(msg)

            # Clear the done flag
            history_info[done_flag] = False

            # Clear the associated file-based flag (if applicable)
            if file_flag and history_info.get(file_flag):
                msg2 = "PERCEIVE: Also clearing %s (output missing)" % file_flag
                diagnostics.append(msg2)
                if log_func:
                    log_func(msg2)
                history_info[file_flag] = False

            # For count-based flags: refine_done / rsr_done imply count > 0
            # If the output is missing, decrement the count by 1 (minimum 0)
            if done_flag == "refine_done":
                old = history_info.get("refine_count", 0)
                history_info["refine_count"] = max(0, old - 1)
                if old > 0:
                    diagnostics.append(
                        "PERCEIVE: refine_count decremented %d→%d (refine output missing)" % (
                            old, history_info["refine_count"]))
            elif done_flag == "rsr_done":
                old = history_info.get("rsr_count", 0)
                history_info["rsr_count"] = max(0, old - 1)
                if old > 0:
                    diagnostics.append(
                        "PERCEIVE: rsr_count decremented %d→%d (RSR output missing)" % (
                            old, history_info["rsr_count"]))

    return diagnostics


# =============================================================================
# WORKFLOW STATE DETECTION
# =============================================================================


def detect_workflow_state(history, available_files, analysis=None, maximum_automation=True,
                         use_yaml_engine=True, directives=None, session_info=None,
                         files_local=True):
    """
    Determine current workflow state based on history and files.

    This function delegates to the YAML-driven WorkflowEngine for state detection.

    Args:
        history: List of cycle records from client
        available_files: List of available file paths
        analysis: Current log analysis dict (optional)
        maximum_automation: If True, use fully automated cryo-EM path
        use_yaml_engine: If True, use YAML-driven WorkflowEngine (default: True)
        directives: Optional user directives dict

    Returns:
        dict: {
            state: str,              # State name
            experiment_type: str,    # "xray" or "cryoem"
            valid_programs: list,    # Programs allowed in this state
            reason: str,             # Human-readable explanation
            conditions: dict,        # Conditional program availability
            automation_path: str,    # "stepwise" or "automated" (cryo-EM only)
            categorized_files: dict, # Pre-categorized files (full_map, half_map, etc.)
        }
    """
    # Categorize files
    # Extract ligand hints from client-side best_files so the server can
    # correctly classify small-molecule PDBs even without disk access.
    # On the client, _pdb_is_small_molecule() reads file content directly.
    # On the server, it cannot read client files, so we use best_files
    # (computed on the client) as an authoritative hint.
    ligand_hints = None
    if session_info:
        best_files = session_info.get("best_files", {})
        if best_files:
            ligand_hints = set()
            for key in ("ligand", "ligand_cif", "ligand_pdb"):
                val = best_files.get(key)
                if val:
                    if isinstance(val, list):
                        for v in val:
                            if v:
                                ligand_hints.add(os.path.basename(v))
                    else:
                        ligand_hints.add(os.path.basename(val))
    files = _categorize_files(available_files, ligand_hints=ligand_hints,
                              files_local=files_local)

    # Analyze history
    history_info = _analyze_history(history)

    # J5: Zombie state detection — clear done flags whose output files are missing.
    # This handles crashed/killed runs where history recorded success but the
    # output file was never fully written (or was subsequently deleted).
    # Diagnostic messages are passed back via the state dict so PERCEIVE can log
    # them; silently discarding them would make re-runs confusing ("why did that
    # program run again?").
    zombie_diagnostics = _clear_zombie_done_flags(history_info, available_files)

    # Determine experiment type
    experiment_type = _detect_experiment_type(files, history_info)

    # Use YAML-driven workflow engine
    if use_yaml_engine:
        try:
            # Lazy import to avoid circular dependencies
            from libtbx.langchain.agent.workflow_engine import WorkflowEngine

            engine = WorkflowEngine()
            state = engine.get_workflow_state(experiment_type, files, history_info, analysis,
                                             directives, maximum_automation, session_info,
                                             files_local=files_local)

            state["categorized_files"] = state.get("categorized_files", files)  # S2c: respect promoted files
            # Set automation_path for both experiment types
            state["automation_path"] = "automated" if maximum_automation else "stepwise"
            # Surface zombie diagnostics for PERCEIVE logging (non-empty only)
            if zombie_diagnostics:
                state["zombie_diagnostics"] = zombie_diagnostics
            return state
        except Exception as e:
            import sys
            print("Warning: YAML workflow engine failed: %s" % e, file=sys.stderr)

    # Fallback: return minimal state (should not happen if YAML is properly configured)
    return {
        "state": "unknown",
        "experiment_type": experiment_type,
        "valid_programs": ["STOP"],
        "reason": "Workflow engine unavailable",
        "conditions": {},
        "automation_path": "automated" if maximum_automation else "stepwise",
        "categorized_files": files,
    }


# =============================================================================
# VALIDATION
# =============================================================================

def validate_program_choice(chosen_program, workflow_state):
    """
    Validate that a program choice is allowed in the current state.

    Args:
        chosen_program: Program the LLM chose
        workflow_state: Dict from detect_workflow_state()

    Returns:
        tuple: (is_valid: bool, error_message: str or None)
    """
    if chosen_program is None:
        return True, None

    if chosen_program == "STOP":
        return True, None

    valid = workflow_state["valid_programs"]

    if chosen_program in valid:
        # P3 fix: autobuild prerequisite check for X-ray experiments.
        # Even when autobuild is in valid_programs (some analyze-step YAML
        # variants include it), it requires phasing to have been done.
        # Without a phased MTZ, autobuild fails with 'MTZ lacks phase/FOM
        # columns'. CASP7-T0283-mr confirmed autobuild ran before phaser
        # on cycle 1 because the 5 ab-initio PDBs + MTZ satisfied file
        # checks; phasing was never done.
        if chosen_program == "phenix.autobuild":
            context = workflow_state.get("context", {})
            exp_type = workflow_state.get("experiment_type", "")
            if exp_type == "xray":
                phasing_done = (
                    context.get("phaser_done") or
                    context.get("autosol_done")
                )
                placed_from_history = context.get(
                    "has_placed_model_from_history")
                if not phasing_done and not placed_from_history:
                    error = (
                        "phenix.autobuild requires a phased MTZ "                        "(PHIB/FOM columns). Run phenix.phaser or "                        "phenix.autosol first to obtain phase information. "                        "Current state: phaser_done=%s, autosol_done=%s."
                    ) % (
                        context.get("phaser_done", False),
                        context.get("autosol_done", False),
                    )
                    return False, error
        # P6 fix: predict_and_build guard for autosol.
        # If a sequence file + model are present and predict_and_build has not
        # yet run, the correct workflow is predict_and_build (MR with a predicted
        # model), not autosol (SAD/MAD experimental phasing).  Block autosol
        # unless the user explicitly requested MR-SAD via use_mr_sad.
        if chosen_program == "phenix.autosol":
            context = workflow_state.get("context", {})
            if (context.get("has_sequence") and
                    context.get("has_model_for_mr") and
                    not context.get("predict_full_done") and
                    not context.get("use_mr_sad")):
                return False, (
                    "predict_and_build workflow detected: a sequence file and "
                    "model are present and predict_and_build has not yet run. "
                    "Use phenix.predict_and_build (MR with predicted model) "
                    "instead of phenix.autosol (SAD/MAD experimental phasing). "
                    "To override, set use_mr_sad=true in workflow_preferences."
                )
        return True, None

    try:
        from libtbx.langchain.knowledge.yaml_loader import get_all_programs
        all_known_programs = get_all_programs()
    except Exception:
        all_known_programs = []  # Graceful degradation

    if chosen_program in all_known_programs:
        error = (
            "Program '%s' is not valid in state '%s'. "
            "Valid programs: %s. Reason: %s"
        ) % (
            chosen_program,
            workflow_state["state"],
            ", ".join(valid),
            workflow_state["reason"]
        )
    else:
        error = "Unknown program '%s'. Valid programs: %s" % (chosen_program, ", ".join(valid))

    return False, error


# =============================================================================
# PROMPT FORMATTING
# =============================================================================

def format_workflow_for_prompt(workflow_state):
    """
    Format workflow state for inclusion in LLM prompt.

    Args:
        workflow_state: Output from detect_workflow_state()

    Returns:
        str: Formatted text for prompt
    """
    lines = []

    lines.append("### WORKFLOW STATE: %s" % workflow_state["state"])
    lines.append("Experiment type: %s" % workflow_state["experiment_type"])

    if workflow_state.get("automation_path"):
        lines.append("Automation path: %s" % workflow_state["automation_path"])

    lines.append("")
    lines.append(workflow_state["reason"])
    lines.append("")
    lines.append("**VALID PROGRAMS FOR THIS STATE:**")
    lines.append(", ".join(workflow_state["valid_programs"]))
    lines.append("")
    lines.append("⚠️ You MUST choose a program from the list above, or set \"stop\": true.")
    lines.append("Choosing an invalid program will cause a validation error.")

    if workflow_state.get("conditions"):
        lines.append("")
        lines.append("Conditional availability:")
        for prog, condition in workflow_state["conditions"].items():
            lines.append("  - %s: requires %s" % (prog, condition))

    # Show stepwise mode hint for both cryo-EM and X-ray
    if workflow_state.get("automation_path") == "stepwise":
        stepwise_states = ["cryoem_analyzed", "xray_initial", "xray_placed"]
        if workflow_state["state"] in stepwise_states:
            lines.append("")
            lines.append("NOTE (Stepwise mode): predict_and_build will use stop_after_predict=true")

    return "\n".join(lines)
