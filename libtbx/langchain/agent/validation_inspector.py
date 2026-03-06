"""
Headless structural validation for the AI agent (v113).

Runs the same validation the Refine GUI uses, but without
any GUI imports. Returns a plain Python dict of structured
results that the thinking agent can reason about.

Entry point: run_validation(model_path, ...) -> dict or None.

Design principle: this function NEVER raises. If validation
fails for any reason, it returns None and the agent falls
back to its existing log-only analysis. This ensures that
wiring it in cannot break the agent.

2-space indentation, 80-char line width.
"""

from __future__ import absolute_import, division, print_function

import os
import traceback
import logging

logger = logging.getLogger(__name__)


def run_validation(model_path, data_path=None,
                   map_coeffs_path=None,
                   sequence_path=None,
                   experiment_type="xray"):
  """Run headless validation on a model.

  Args:
    model_path: str, path to PDB/mmCIF model file.
    data_path: str or None, path to MTZ with Fobs/R-free.
    map_coeffs_path: str or None, path to MTZ with map
      coefficients (for difference density peak search).
    sequence_path: str or None, path to sequence file.
    experiment_type: "xray" or "cryoem".

  Returns:
    dict with structured validation results, or None if
    validation fails. Never raises — logs errors and
    returns None so the agent continues with log-only
    analysis (the current fallback behavior).
  """
  if not model_path or not os.path.isfile(model_path):
    logger.warning(
      "Validation skipped: model not found: %s",
      model_path,
    )
    return None
  try:
    return _run_validation_inner(
      model_path=model_path,
      data_path=data_path,
      map_coeffs_path=map_coeffs_path,
      sequence_path=sequence_path,
      experiment_type=experiment_type,
    )
  except Exception:
    logger.warning(
      "Validation failed, continuing with log-only "
      "analysis:\n%s", traceback.format_exc()
    )
    return None


def _run_validation_inner(model_path, data_path,
                          map_coeffs_path,
                          sequence_path,
                          experiment_type):
  """Inner validation. May raise."""
  from iotbx.data_manager import DataManager
  dm = DataManager()
  dm.process_model_file(model_path)
  model = dm.get_model()

  result = {}

  # Model contents inventory (chains, ligands, waters)
  result["model_contents"] = _inspect_model(model)

  # Geometry (Ramachandran, rotamers, clashscore)
  result["geometry"] = _validate_geometry(model)

  # Per-chain completeness (Step 1.3)
  result["chain_completeness"] = (
    _chain_completeness(model, sequence_path)
  )

  # Data-dependent validation (ligand RSCC, etc.)
  if data_path and experiment_type == "xray":
    result["data_model"] = _validate_data_model(
      model, data_path, map_coeffs_path
    )

  # Difference density peaks
  if map_coeffs_path:
    result["diff_peaks"] = _find_diff_peaks(
      model, map_coeffs_path
    )

  return result


# --- Stubs for subsequent steps (A2-A5) ---

def _inspect_model(model):
  """Extract model contents inventory.

  Walks the PDB hierarchy to identify chains, standard
  residues, ligands (HETATM with >3 atoms), waters,
  and common ions. Uses the atom.hetero flag as the
  primary discriminator for ligand detection.

  Returns dict:
    chains: list of chain IDs
    residue_count: int (protein + nucleic acid)
    ligands: list of {name, chain, resid, n_atoms}
    waters: int
    ions: list of {name, chain, resid}
    has_hetatm: bool
  """
  hierarchy = model.get_hierarchy()
  chains = []
  ligands = []
  ions = []
  water_count = 0
  residue_count = 0

  _WATER = frozenset(["HOH", "WAT", "DOD"])
  _IONS = frozenset([
    "ZN", "MG", "CA", "FE", "MN", "CU", "NI",
    "CO", "NA", "K", "CL", "CD", "HG", "PT",
    "AU", "AG", "PB", "SR", "BA", "RB", "CS",
  ])
  _STD_AA = frozenset([
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN",
    "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
    "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "MSE", "UNK", "SEC", "PYL",
  ])
  _STD_NA = frozenset([
    "DA", "DT", "DG", "DC", "DU",
    "A", "T", "G", "C", "U",
  ])

  # Track ligand (name, chain, resid) to avoid dupes
  # from alt confs
  seen_ligands = set()

  for chain_obj in hierarchy.chains():
    cid = chain_obj.id.strip()
    if cid and cid not in chains:
      chains.append(cid)
    for rg in chain_obj.residue_groups():
      for ag in rg.atom_groups():
        resname = ag.resname.strip()
        resid = rg.resseq_as_int()
        atoms = ag.atoms()
        n_atoms = len(atoms)

        if resname in _WATER:
          water_count += 1
        elif resname in _IONS and n_atoms <= 2:
          ions.append({
            "name": resname,
            "chain": cid,
            "resid": resid,
          })
        elif resname in _STD_AA or resname in _STD_NA:
          residue_count += 1
        elif n_atoms > 3:
          # Non-standard, multi-atom -> ligand candidate
          is_het = any(a.hetero for a in atoms)
          if is_het:
            key = (resname, cid, resid)
            if key not in seen_ligands:
              seen_ligands.add(key)
              ligands.append({
                "name": resname,
                "chain": cid,
                "resid": resid,
                "n_atoms": n_atoms,
              })
        else:
          # Small non-standard (<=3 atoms, not ion)
          # e.g. ACE, NH2 caps, buffer artifacts
          # Count as residue if not HETATM
          is_het = any(a.hetero for a in atoms)
          if not is_het:
            residue_count += 1

  return {
    "chains": chains,
    "residue_count": residue_count,
    "ligands": ligands,
    "waters": water_count,
    "ions": ions,
    "has_hetatm": len(ligands) > 0 or len(ions) > 0,
  }


def _validate_geometry(model):
  """Run geometry validation (Ramachandran, rotamers, clashscore).

  Calls the same mmtbx validation that MolProbity uses.
  The Refine GUI displays these metrics (Output.py line
  549-551), confirming they're available headlessly.

  Bonds/angles RMSD require a restraints manager, which
  may not be available for all model types. Falls back
  to None if unavailable (log_analysis already extracts
  them from the refinement log).

  Returns dict:
    rama_favored: float (fraction, 0.0-1.0)
    rama_outliers: float (fraction)
    rama_outlier_list: list of "chain/resname resid"
    rotamer_outliers: float (fraction)
    rotamer_outlier_list: list of "chain/resname resid"
    clashscore: float
    bonds_rmsd: float or None
    angles_rmsd: float or None
  """
  from mmtbx.validation import ramalyze
  from mmtbx.validation import rotalyze
  from mmtbx.validation import clashscore as cs_mod

  hierarchy = model.get_hierarchy()

  # --- Ramachandran ---
  rama = ramalyze.ramalyze(pdb_hierarchy=hierarchy)
  rama_outlier_list = []
  for r in rama.results:
    if r.is_outlier():
      rama_outlier_list.append(
        "%s/%s%s" % (
          r.chain_id.strip(),
          r.resname.strip(),
          r.resseq.strip(),
        )
      )

  # --- Rotamers ---
  rota = rotalyze.rotalyze(pdb_hierarchy=hierarchy)
  rota_outlier_list = []
  for r in rota.results:
    if r.is_outlier():
      rota_outlier_list.append(
        "%s/%s%s" % (
          r.chain_id.strip(),
          r.resname.strip(),
          r.resseq.strip(),
        )
      )

  # Get counts with version-safe attribute access
  n_rota = getattr(rota, 'n_total', None)
  if n_rota is None:
    n_rota = len(rota.results)
  n_rota_out = getattr(rota, 'n_outliers', None)
  if n_rota_out is None:
    n_rota_out = len(rota_outlier_list)

  # --- Clashscore ---
  clash = cs_mod.clashscore(pdb_hierarchy=hierarchy)
  clashscore_val = clash.get_clashscore()

  # --- Bonds and angles (optional) ---
  bonds_rmsd = None
  angles_rmsd = None
  try:
    grm = model.get_restraints_manager()
    if grm is not None:
      energies = grm.energies_sites(
        sites_cart=model.get_sites_cart(),
        compute_gradients=False,
      )
      bonds_rmsd = energies.bond_deviations()[2]
      angles_rmsd = energies.angle_deviations()[2]
  except Exception:
    # Not critical — log_analysis has these from the
    # refinement log. Silently skip.
    pass

  return {
    "rama_favored": rama.percent_favored / 100.0,
    "rama_outliers": rama.percent_outliers / 100.0,
    "rama_outlier_list": rama_outlier_list[:10],
    "rotamer_outliers": (
      n_rota_out / max(n_rota, 1)
    ),
    "rotamer_outlier_list": rota_outlier_list[:10],
    "clashscore": clashscore_val,
    "bonds_rmsd": bonds_rmsd,
    "angles_rmsd": angles_rmsd,
  }


def _chain_completeness(model, sequence_path=None):
  """Per-chain completeness assessment (Step 1.3).

  Walks the model hierarchy to report built residues,
  gaps in numbering, average B-factor, and disordered
  fraction per chain. If a sequence file is provided,
  compares against expected residue counts.

  "Disordered fraction" = fraction of residues with
  average B-factor > 2x the chain average, OR with
  missing main-chain atoms (no CA).

  Args:
    model: cctbx model object.
    sequence_path: str or None. Path to FASTA/PIR
      sequence file for expected residue counts.

  Returns:
    list of {chain_id, residues_built,
      residues_expected, completeness_fraction,
      gaps: [{start, end}], avg_b_factor,
      disordered_fraction}
  """
  try:
    return _chain_completeness_inner(
      model, sequence_path
    )
  except Exception:
    logger.debug(
      "_chain_completeness failed",
      exc_info=True,
    )
    return []


def _chain_completeness_inner(model, sequence_path):
  """Inner chain completeness. May raise."""
  hierarchy = model.get_hierarchy()

  # --- Read expected chain lengths from sequence ---
  expected_lengths = {}
  if sequence_path and os.path.isfile(sequence_path):
    expected_lengths = _read_sequence_lengths(
      sequence_path
    )

  # Standard residue names (protein + nucleic acid)
  _STD_AA = frozenset([
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN",
    "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
    "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "MSE", "UNK", "SEC", "PYL",
  ])
  _STD_NA = frozenset([
    "DA", "DT", "DG", "DC", "DU",
    "A", "T", "G", "C", "U",
  ])

  results = []
  for chain_obj in hierarchy.chains():
    cid = chain_obj.id.strip()
    if not cid:
      continue

    residue_numbers = []
    b_factors = []
    disordered_indices = set()  # track by index
    n_residues = 0

    for rg in chain_obj.residue_groups():
      # Use first atom_group (skip altlocs)
      ags = rg.atom_groups()
      if not ags:
        continue
      ag = ags[0]
      resname = ag.resname.strip()
      if resname not in _STD_AA and \
          resname not in _STD_NA:
        continue

      res_idx = n_residues
      n_residues += 1
      resid = rg.resseq_as_int()
      residue_numbers.append(resid)

      # Per-residue average B-factor
      atoms = ag.atoms()
      if atoms.size() > 0:
        b_vals = atoms.extract_b()
        avg_b = sum(b_vals) / len(b_vals)
        b_factors.append(avg_b)

        # Check for missing main-chain atoms
        atom_names = set(
          a.name.strip() for a in atoms
        )
        has_ca = "CA" in atom_names
        if not has_ca:
          disordered_indices.add(res_idx)
      else:
        b_factors.append(0.0)

    if n_residues == 0:
      continue

    # --- Compute chain-level statistics ---
    chain_avg_b = (
      sum(b_factors) / len(b_factors)
      if b_factors else 0.0
    )

    # Also flag residues with high B-factor
    if chain_avg_b > 0:
      for idx, b in enumerate(b_factors):
        if b > 2.0 * chain_avg_b:
          disordered_indices.add(idx)
    disordered_frac = (
      len(disordered_indices) / max(n_residues, 1)
    )

    # --- Find gaps in residue numbering ---
    gaps = []
    if residue_numbers:
      residue_numbers.sort()
      for i in range(len(residue_numbers) - 1):
        curr = residue_numbers[i]
        nxt = residue_numbers[i + 1]
        if nxt - curr > 1:
          gaps.append({
            "start": curr + 1,
            "end": nxt - 1,
          })

    # --- Expected residues ---
    expected = expected_lengths.get(cid)
    if expected is None and residue_numbers:
      # Estimate from numbering range
      expected = (
        residue_numbers[-1]
        - residue_numbers[0] + 1
      )
    completeness = (
      n_residues / max(expected, 1)
      if expected else None
    )

    results.append({
      "chain_id": cid,
      "residues_built": n_residues,
      "residues_expected": expected,
      "completeness_fraction": (
        round(completeness, 3)
        if completeness is not None else None
      ),
      "gaps": gaps[:10],  # cap for sanity
      "avg_b_factor": round(chain_avg_b, 1),
      "disordered_fraction": round(
        disordered_frac, 3
      ),
    })

  return results


def _read_sequence_lengths(sequence_path):
  """Read chain lengths from a FASTA/PIR sequence file.

  Returns dict of {chain_id: residue_count}.
  Very forgiving parser — skips malformed entries.
  """
  result = {}
  try:
    from iotbx import bioinformatics
    seqs = bioinformatics.any_sequence_format(
      file_name=sequence_path
    )
    if seqs and hasattr(seqs, 'sequence'):
      # Single sequence
      result["A"] = len(seqs.sequence)
    elif seqs:
      # Multiple sequences
      for i, seq_obj in enumerate(seqs):
        seq = getattr(seq_obj, 'sequence', '')
        name = getattr(seq_obj, 'name', '')
        # Try to extract chain ID from name
        cid = name.strip()[:1].upper() if name \
          else chr(65 + i)  # A, B, C...
        if seq:
          result[cid] = len(seq)
  except Exception:
    # Non-critical — completeness will be estimated
    # from residue numbering instead.
    pass
  return result


def _validate_data_model(model, data_path,
                         map_coeffs_path=None):
  """Compute per-ligand real-space correlation (Step 1.2).

  Uses the 2mFo-DFc map from the refinement output MTZ
  (preferred) or computes it via mmtbx.f_model.manager
  from Fobs+model (fallback). Computes per-residue CC
  and filters to ligands.

  Also computes a Z-score relative to the protein CC
  distribution for resolution-independent assessment.

  R-factors are NOT computed here — they come from
  log_analysis which is more authoritative.

  Args:
    model: cctbx model object.
    data_path: str, path to MTZ with Fobs/R-free.
    map_coeffs_path: str or None, path to MTZ with
      2mFo-DFc and mFo-DFc map coefficients.

  Returns dict:
    ligand_rscc: list of {name, chain, resid, rscc,
      rscc_z, b_mean, occupancy}
  """
  result = {"ligand_rscc": []}
  try:
    return _validate_data_model_inner(
      model, data_path, map_coeffs_path,
    )
  except Exception:
    logger.debug(
      "_validate_data_model failed",
      exc_info=True,
    )
    return result


def _validate_data_model_inner(model, data_path,
                               map_coeffs_path):
  """Inner ligand RSCC. May raise."""
  # Import with fallback — module path varies
  # across PHENIX versions.
  try:
    from mmtbx.real_space import \
      real_space_correlation
  except ImportError:
    from mmtbx import real_space_correlation

  hierarchy = model.get_hierarchy()
  xrs = model.get_xray_structure()

  # --- Strategy 1: Read 2mFo-DFc from map coeffs ---
  # Preferred: uses the refinement's own map, which
  # includes proper bulk solvent and scaling.
  map_data = None
  d_min = None
  if map_coeffs_path and os.path.isfile(
    map_coeffs_path
  ):
    map_data, d_min = _read_2fofc_map(
      map_coeffs_path
    )

  # --- Strategy 2: Compute via f_model.manager ---
  # Fallback when map_coeffs_path is not available.
  if map_data is None:
    map_data, d_min = _compute_2fofc_map(
      model, data_path
    )

  if map_data is None:
    return {"ligand_rscc": []}

  # --- Per-residue real-space CC ---
  rsc_result = real_space_correlation.simple(
    fmodel=None,
    pdb_hierarchy=hierarchy,
    map_data=map_data,
    unit_cell=xrs.unit_cell(),
    d_min=d_min,
  )

  # Collect protein CCs for Z-score computation
  protein_ccs = []
  _STD_AA = frozenset([
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN",
    "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
    "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "MSE",
  ])

  for r in rsc_result:
    resname = getattr(r, 'resname', '').strip()
    if resname in _STD_AA:
      cc = getattr(r, 'cc', None)
      if cc is not None:
        protein_ccs.append(float(cc))

  # Compute mean/std for Z-score
  if len(protein_ccs) >= 10:
    mean_cc = sum(protein_ccs) / len(protein_ccs)
    variance = sum(
      (x - mean_cc) ** 2 for x in protein_ccs
    ) / len(protein_ccs)
    std_cc = variance ** 0.5
  else:
    mean_cc = None
    std_cc = None

  # --- Identify ligand residues ---
  _WATER = frozenset(["HOH", "WAT", "DOD"])
  ligand_results = []
  for r in rsc_result:
    resname = getattr(r, 'resname', '').strip()
    chain = getattr(r, 'chain_id', '').strip()
    resid = getattr(r, 'resseq', 0)
    if isinstance(resid, str):
      try:
        resid = int(resid.strip())
      except (ValueError, TypeError):
        resid = 0
    cc = getattr(r, 'cc', None)
    b = getattr(r, 'b', None)
    occ = getattr(r, 'occupancy', None)
    n_atoms = getattr(r, 'n_atoms', 0)

    # Skip standard residues, waters, tiny groups
    if resname in _STD_AA or resname in _WATER:
      continue
    if n_atoms is not None and n_atoms <= 3:
      continue
    if cc is None:
      continue

    cc_float = float(cc)

    # Z-score relative to protein distribution
    z_score = None
    if mean_cc is not None and std_cc and \
        std_cc > 0.001:
      z_score = (cc_float - mean_cc) / std_cc

    ligand_results.append({
      "name": resname,
      "chain": chain,
      "resid": resid,
      "rscc": round(cc_float, 3),
      "rscc_z": (
        round(z_score, 2)
        if z_score is not None else None
      ),
      "b_mean": (
        round(float(b), 1) if b is not None
        else None
      ),
      "occupancy": (
        round(float(occ), 2) if occ is not None
        else None
      ),
    })

  return {"ligand_rscc": ligand_results}


def _read_2fofc_map(map_coeffs_path):
  """Read 2mFo-DFc map from a refinement output MTZ.

  Looks for FWT/PHFWT, 2FOFCWT/PH2FOFCWT, or similar
  column names. Excludes mFo-DFc columns (DELFWT,
  FOFCWT).

  Returns:
    (map_data, d_min) or (None, None) if not found.
  """
  from iotbx import mtz as iotbx_mtz

  mtz_obj = iotbx_mtz.object(
    file_name=map_coeffs_path
  )
  # 2mFo-DFc column names by priority.
  # Check individual labels, not substrings.
  _2FOFC_NAMES = frozenset([
    "FWT", "2FOFCWT", "2FOFC", "2MFO-DFC",
  ])
  # mFo-DFc column names to EXCLUDE.
  _FOFC_NAMES = frozenset([
    "DELFWT", "FOFCWT", "FOFC", "MFO-DFC",
  ])

  two_fofc = None
  for ma in mtz_obj.as_miller_arrays():
    if not ma.is_complex_array():
      continue
    col_names = set(
      str(l).strip().upper()
      for l in ma.info().labels
    )
    # Skip if any column is an mFo-DFc name
    if col_names & _FOFC_NAMES:
      continue
    # Match if any column is a 2mFo-DFc name
    if col_names & _2FOFC_NAMES:
      two_fofc = ma
      break

  if two_fofc is None:
    return (None, None)

  fft_map = two_fofc.fft_map(
    resolution_factor=0.25,
  )
  fft_map.apply_sigma_scaling()
  return (
    fft_map.real_map_unpadded(),
    two_fofc.d_min(),
  )


def _compute_2fofc_map(model, data_path):
  """Compute 2mFo-DFc map via f_model.manager.

  Fallback when refinement map coefficients are not
  available. Handles bulk solvent and scaling properly.

  Returns:
    (map_data, d_min) or (None, None) on failure.
  """
  try:
    from mmtbx import f_model as f_model_mod
    from iotbx import mtz as iotbx_mtz

    # Read Fobs
    mtz_obj = iotbx_mtz.object(file_name=data_path)
    fobs = None
    for ma in mtz_obj.as_miller_arrays():
      labels = ",".join(
        str(l) for l in ma.info().labels
      ).upper()
      if ma.is_xray_amplitude_array() or \
          ma.is_xray_intensity_array():
        if "FREE" not in labels:
          fobs = ma.as_amplitude_array()
          break
    if fobs is None:
      return (None, None)

    xrs = model.get_xray_structure()
    fmodel = f_model_mod.manager(
      f_obs=fobs,
      xray_structure=xrs,
    )
    fmodel.update_all_scales()

    # Get 2mFo-DFc map coefficients
    map_coeffs = fmodel.map_coefficients(
      map_type="2mFo-DFc"
    )
    fft_map = map_coeffs.fft_map(
      resolution_factor=0.25,
    )
    fft_map.apply_sigma_scaling()
    return (
      fft_map.real_map_unpadded(),
      fobs.d_min(),
    )
  except Exception:
    logger.debug(
      "_compute_2fofc_map failed", exc_info=True,
    )
    return (None, None)


def _find_diff_peaks(model, map_coeffs_path):
  """Find significant peaks in mFo-DFc map (Step 1.1).

  Uses the mFo-DFc map coefficients from a refinement
  output MTZ to find positive peaks (unmodeled features)
  and negative holes (wrongly-placed atoms).

  For each peak, finds the nearest residue using
  brute-force distance search on model atom positions.

  Args:
    model: cctbx model object.
    map_coeffs_path: str, path to MTZ with map
      coefficients (FOFCWT/PHFOFCWT columns).

  Returns dict:
    positive: list of {height, xyz, near_residue,
      near_chain, distance} — peaks > 4.0 sigma,
      sorted by height descending. Capped at 20.
    negative: list of same format — holes < -4.0
      sigma, sorted by depth. Capped at 10.
    peak_count: int (total positive peaks > 3.5 sigma,
      before the cap).
  """
  empty = {"positive": [], "negative": [],
           "peak_count": 0}
  try:
    return _find_diff_peaks_inner(
      model, map_coeffs_path
    )
  except Exception:
    logger.debug(
      "_find_diff_peaks failed", exc_info=True,
    )
    return empty


# Module-level cache for KD-tree / atom lookup.
# Keyed by (path, mtime, size, space_group, unit_cell).
_atom_lookup_cache = {}
_atom_lookup_cache_key = None


def _find_diff_peaks_inner(model, map_coeffs_path):
  """Inner diff peaks. May raise."""
  from iotbx import mtz as iotbx_mtz
  from cctbx import maptbx
  from scitbx.array_family import flex

  if not os.path.isfile(map_coeffs_path):
    return {"positive": [], "negative": [],
            "peak_count": 0}

  # --- Read mFo-DFc map coefficients from MTZ ---
  mtz_obj = iotbx_mtz.object(
    file_name=map_coeffs_path
  )

  # mFo-DFc column names by priority.
  _FOFC_NAMES = frozenset([
    "DELFWT", "FOFCWT", "FOFC", "MFO-DFC",
  ])
  # 2mFo-DFc column names to EXCLUDE (we want the
  # difference map, not the weighted map).
  _2FOFC_NAMES = frozenset([
    "FWT", "2FOFCWT", "2FOFC", "2MFO-DFC",
  ])

  fofc_coeffs = None
  for ma in mtz_obj.as_miller_arrays():
    if not ma.is_complex_array():
      continue
    col_names = set(
      str(l).strip().upper()
      for l in ma.info().labels
    )
    # Skip 2mFo-DFc columns
    if col_names & _2FOFC_NAMES:
      continue
    # Match mFo-DFc columns
    if col_names & _FOFC_NAMES:
      fofc_coeffs = ma
      break

  if fofc_coeffs is None:
    # Fallback: any complex array with "DEL" in a
    # column name (e.g. DELFWT)
    for ma in mtz_obj.as_miller_arrays():
      if not ma.is_complex_array():
        continue
      col_names = set(
        str(l).strip().upper()
        for l in ma.info().labels
      )
      if any("DEL" in c for c in col_names):
        fofc_coeffs = ma
        break

  if fofc_coeffs is None:
    logger.debug(
      "_find_diff_peaks: no mFo-DFc columns found "
      "in %s", map_coeffs_path,
    )
    return {"positive": [], "negative": [],
            "peak_count": 0}

  # Ensure complex array for FFT
  if not fofc_coeffs.is_complex_array():
    try:
      fofc_coeffs = fofc_coeffs.as_complex_array()
    except Exception:
      return {"positive": [], "negative": [],
              "peak_count": 0}

  # --- FFT to real-space map ---
  fft_map = fofc_coeffs.fft_map(
    resolution_factor=0.25,
  )
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()

  crystal_symmetry = fofc_coeffs.crystal_symmetry()
  uc = crystal_symmetry.unit_cell()

  # Tags array: all zeros = search everywhere.
  # Required by maptbx.peak_list (not optional).
  tags = flex.int(map_data.accessor(), 0)

  # Positive peaks (sorted by height descending)
  pos_peaks_raw = maptbx.peak_list(
    data=map_data,
    tags=tags,
    peak_search_level=1,
    max_peaks=100,
    interpolate=True,
  )

  # --- Build atom lookup for nearest-residue ---
  atom_info = _get_atom_lookup(
    model, map_coeffs_path, crystal_symmetry
  )
  sites_cart = atom_info["sites_cart"]
  atom_labels = atom_info["labels"]

  # --- Process positive peaks ---
  positive = []
  peak_count_35 = 0  # total > 3.5σ
  heights = pos_peaks_raw.heights()
  sites = pos_peaks_raw.sites()
  n_peaks = len(heights)
  for i in range(n_peaks):
    height = heights[i]
    if height < 3.5:
      break
    peak_count_35 += 1
    if height >= 4.0 and len(positive) < 20:
      frac = sites[i]
      cart = list(uc.orthogonalize(frac))
      near = _find_nearest_residue(
        cart, sites_cart, atom_labels
      )
      positive.append({
        "height": round(float(height), 1),
        "xyz": [round(c, 1) for c in cart],
        "near_residue": near["residue"],
        "near_chain": near["chain"],
        "distance": near["distance"],
      })

  # --- Process negative peaks ---
  # Invert the map to find holes as peaks.
  # deep_copy ensures we don't modify the original.
  neg_map = map_data.deep_copy() * (-1.0)
  neg_tags = flex.int(neg_map.accessor(), 0)
  neg_peaks_raw = maptbx.peak_list(
    data=neg_map,
    tags=neg_tags,
    peak_search_level=1,
    max_peaks=50,
    interpolate=True,
  )

  negative = []
  neg_heights = neg_peaks_raw.heights()
  neg_sites = neg_peaks_raw.sites()
  n_neg = len(neg_heights)
  for i in range(n_neg):
    depth = neg_heights[i]
    if depth < 4.0:
      break
    if len(negative) >= 10:
      break
    frac = neg_sites[i]
    cart = list(uc.orthogonalize(frac))
    near = _find_nearest_residue(
      cart, sites_cart, atom_labels
    )
    negative.append({
      "height": round(-float(depth), 1),
      "xyz": [round(c, 1) for c in cart],
      "near_residue": near["residue"],
      "near_chain": near["chain"],
      "distance": near["distance"],
    })

  return {
    "positive": positive,
    "negative": negative,
    "peak_count": peak_count_35,
  }


def _get_atom_lookup(model, map_coeffs_path,
                     crystal_symmetry):
  """Build or retrieve cached atom position lookup.

  Cache key: (path, mtime, size, space_group,
  unit_cell). Rebuilds if any component changes.

  Returns dict with:
    sites_cart: flex.vec3_double of atom positions
    labels: list of {chain, resname, resid} per atom
  """
  global _atom_lookup_cache, _atom_lookup_cache_key

  # Build cache key
  sg_str = str(
    crystal_symmetry.space_group_info()
  ) if crystal_symmetry else ""
  uc_params = tuple(
    crystal_symmetry.unit_cell().parameters()
  ) if crystal_symmetry else ()
  try:
    stat = os.stat(map_coeffs_path)
    file_key = (
      map_coeffs_path, stat.st_mtime,
      stat.st_size,
    )
  except OSError:
    file_key = (map_coeffs_path, 0, 0)
  cache_key = file_key + (sg_str, uc_params)

  if (
    _atom_lookup_cache_key == cache_key
    and _atom_lookup_cache
  ):
    return _atom_lookup_cache

  # Build fresh lookup
  hierarchy = model.get_hierarchy()
  sites_cart = model.get_sites_cart()
  labels = []
  for atom in hierarchy.atoms():
    parent = atom.parent()
    rg = parent.parent() if parent else None
    chain = rg.parent() if rg else None
    labels.append({
      "chain": (
        chain.id.strip() if chain else "?"
      ),
      "resname": (
        parent.resname.strip() if parent
        else "?"
      ),
      "resid": (
        rg.resseq_as_int() if rg else 0
      ),
    })

  result = {
    "sites_cart": sites_cart,
    "labels": labels,
  }
  _atom_lookup_cache = result
  _atom_lookup_cache_key = cache_key
  return result


def _find_nearest_residue(peak_xyz, sites_cart,
                          atom_labels):
  """Find the nearest residue to a peak position.

  Uses vectorized flex arithmetic for efficient
  distance computation across all atoms at once.

  Args:
    peak_xyz: list of [x, y, z].
    sites_cart: flex.vec3_double of atom positions.
    atom_labels: list of {chain, resname, resid}.

  Returns dict:
    residue: str, e.g. "A/His47"
    chain: str
    distance: float (Angstroms)
  """
  import math
  from scitbx.array_family import flex

  if sites_cart.size() == 0:
    return {
      "residue": "unknown",
      "chain": "?",
      "distance": 999.9,
    }

  # Decompose into x, y, z flex.double arrays
  # (.parts() is stable across all cctbx versions)
  xs, ys, zs = sites_cart.parts()
  px = float(peak_xyz[0])
  py = float(peak_xyz[1])
  pz = float(peak_xyz[2])
  dx = xs - px
  dy = ys - py
  dz = zs - pz
  d_sq = dx * dx + dy * dy + dz * dz
  min_idx = flex.min_index(d_sq)
  dist = math.sqrt(d_sq[min_idx])

  if min_idx < len(atom_labels):
    label = atom_labels[min_idx]
    residue = "%s/%s%d" % (
      label["chain"],
      label["resname"],
      label["resid"],
    )
    chain = label["chain"]
  else:
    residue = "unknown"
    chain = "?"

  return {
    "residue": residue,
    "chain": chain,
    "distance": round(dist, 1),
  }
