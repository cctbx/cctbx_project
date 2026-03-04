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

  # Data-dependent validation (ligand RSCC, etc.)
  if data_path and experiment_type == "xray":
    result["data_model"] = _validate_data_model(
      model, data_path
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


def _validate_data_model(model, data_path):
  """Data-dependent validation metrics.

  Currently returns a stub. Future work (follow-up F1)
  will compute per-ligand real-space CC using
  mmtbx.real_space_correlation or phenix.get_cc_mtz_pdb.

  R-factors are NOT computed here — they come from
  log_analysis which is more authoritative (reads them
  directly from the refinement log).

  Returns dict:
    ligand_rscc: list of {name, chain, resid, rscc,
      b_mean} — currently empty, filled in F1.
  """
  return {"ligand_rscc": []}


def _find_diff_peaks(model, map_coeffs_path):
  """Find significant peaks in the difference map.

  Currently returns a stub. Future work (follow-up F2)
  will use mmtbx.find_peaks_holes to find peaks >4 sigma
  in the mFo-DFc map and assign each to the nearest
  residue.

  Returns dict:
    positive: list of {height, near_residue}
      — currently empty, filled in F2.
    negative: list of {height, near_residue}
      — currently empty, filled in F2.
  """
  return {"positive": [], "negative": []}
