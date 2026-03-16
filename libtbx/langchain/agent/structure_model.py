"""
Structure Model for the Goal-Directed Agent.

Maintains a running understanding of the specific structure
being solved — not just metrics, but structural knowledge
that persists across cycles and influences decisions.

Updated after each cycle from validation results, log
analysis, and user input. Persisted to session JSON for
resume support.

NOT derived from LLM reasoning — all data comes from
ground truth (validation, log parsing, file inspection).

Entry point: StructureModel class.

2-space indentation, 80-char line width.
"""

from __future__ import absolute_import, division, print_function

import copy
import hashlib
import json
import logging
import time

logger = logging.getLogger(__name__)


# ── Hypothesis (Phase 4 data structure) ─────────────────
#
# Defined here because the StructureModel owns the
# hypothesis list and enforces the single-active budget.
# Phase 4 adds generation and evaluation logic; the data
# structure is stable from Phase 1 onward.

class Hypothesis(object):
  """A testable claim about the structure.

  Lifecycle: proposed → testing → pending → confirmed
                                           | refuted
                                           | abandoned

  The StructureModel enforces at most ONE active hypothesis
  (status in {"testing", "pending"}) at any time.
  """

  VALID_STATUSES = frozenset([
    "proposed", "testing", "pending",
    "confirmed", "refuted", "abandoned",
  ])
  ACTIVE_STATUSES = frozenset(["testing", "pending"])

  def __init__(self, id, statement, test_program="",
               test_parameters=None, confirm_if="",
               refute_if="", status="proposed",
               proposed_at_cycle=0,
               resolved_at_cycle=None,
               test_cycles_remaining=1,
               revalidation_reason=""):
    self.id = str(id)
    self.statement = str(statement)
    self.test_program = str(test_program)
    self.test_parameters = (
      dict(test_parameters) if test_parameters else {}
    )
    self.confirm_if = str(confirm_if)
    self.refute_if = str(refute_if)
    self.status = (
      status if status in self.VALID_STATUSES
      else "proposed"
    )
    self.proposed_at_cycle = int(proposed_at_cycle)
    self.resolved_at_cycle = (
      int(resolved_at_cycle)
      if resolved_at_cycle is not None else None
    )
    self.test_cycles_remaining = max(
      0, int(test_cycles_remaining)
    )
    # Set when a confirmed hypothesis is re-evaluated
    # (Phase 4 re-validation). Empty string if not
    # a re-proposal.
    self.revalidation_reason = str(revalidation_reason)

  @property
  def is_active(self):
    """True if this hypothesis occupies the active slot."""
    return self.status in self.ACTIVE_STATUSES

  @property
  def is_resolved(self):
    """True if confirmed, refuted, or abandoned."""
    return self.status in (
      "confirmed", "refuted", "abandoned"
    )

  def to_dict(self):
    """Serialize to JSON-safe dict."""
    return {
      "id": self.id,
      "statement": self.statement,
      "test_program": self.test_program,
      "test_parameters": dict(self.test_parameters),
      "confirm_if": self.confirm_if,
      "refute_if": self.refute_if,
      "status": self.status,
      "proposed_at_cycle": self.proposed_at_cycle,
      "resolved_at_cycle": self.resolved_at_cycle,
      "test_cycles_remaining": self.test_cycles_remaining,
      "revalidation_reason": self.revalidation_reason,
    }

  @classmethod
  def from_dict(cls, d):
    """Reconstruct from dict. Tolerant of missing keys."""
    if not isinstance(d, dict):
      return cls(id="unknown", statement="unknown")
    return cls(
      id=d.get("id", "unknown"),
      statement=d.get("statement", "unknown"),
      test_program=d.get("test_program", ""),
      test_parameters=d.get("test_parameters"),
      confirm_if=d.get("confirm_if", ""),
      refute_if=d.get("refute_if", ""),
      status=d.get("status", "proposed"),
      proposed_at_cycle=d.get("proposed_at_cycle", 0),
      resolved_at_cycle=d.get("resolved_at_cycle"),
      test_cycles_remaining=d.get(
        "test_cycles_remaining", 1
      ),
      revalidation_reason=d.get(
        "revalidation_reason", ""
      ),
    )

  def __repr__(self):
    return (
      "Hypothesis(id=%r, status=%r, statement=%r)"
      % (self.id, self.status, self.statement[:50])
    )


# ── StructureModel ──────────────────────────────────────

class StructureModel(object):
  """Running understanding of the current structure.

  Updated after each cycle from validation results,
  log analysis, and user input. Persisted to session
  JSON for resume support.

  NOT derived from LLM reasoning — all data comes from
  ground truth (validation, log parsing, file inspection).

  Sections:
    data_characteristics — set early, rarely change
    model_state — updated every cycle from validation
    progress — cumulative per-cycle metric snapshots
    strategy_blacklist — failed strategies (anti-oscillation)
    hypotheses — testable claims (Phase 4)
  """

  def __init__(self):
    # --- Data characteristics (set early, rarely change) ---
    # Populated by update_from_xtriage and update_from_phaser.
    self.data_characteristics = {
      "resolution": None,         # float, Angstroms
      "completeness": None,       # float, 0-1
      "redundancy": None,         # float
      "space_group": None,        # str, e.g. "P 21 21 21"
      "unit_cell": None,          # tuple of 6 floats
      "twinning": {
        "is_twinned": None,       # bool or None
        "twin_law": None,         # str
        "twin_fraction": None,    # float
      },
      "anomalous": {
        "has_anomalous_data": None,   # bool or None
        "anomalous_d_min": None,      # float
        "anomalous_measurability": None,  # float
      },
      "data_quality": "",         # "good", "twinned", etc.
      "experiment_type": None,    # "xray" or "cryoem"
      "i_over_sigma": None,       # float, overall <I/σ(I)>
    }

    # --- Model state (updated every cycle) ---
    # Populated by update_from_validation.
    self.model_state = {
      "chains": [],
      # list of {chain_id, residues_built,
      #          residues_expected, completeness,
      #          gaps, avg_b_factor, disordered_fraction}
      "ligands": [],
      # list of {name, chain, resid, n_atoms,
      #          rscc, z_rscc, b_mean, occupancy}
      "waters": 0,
      "ions": [],
      # list of {name, chain, resid}
      "geometry": {
        "rama_favored": None,
        "rama_outliers": None,
        "rama_outlier_list": [],
        "rotamer_outliers": None,
        "rotamer_outlier_list": [],
        "clashscore": None,
        "bonds_rmsd": None,
        "angles_rmsd": None,
      },
      "diff_peaks": {
        "positive": [],
        # list of {height, xyz, near_residue,
        #          near_chain, distance}
        "negative": [],
        "peak_count": 0,
      },
      "r_work": None,
      "r_free": None,
      "model_map_cc": None,   # cryo-EM overall CC
      "problems": [],
      # list of {region, type, severity, description}
      # severity: "high", "medium", "low"
    }

    # --- Progress tracking (cumulative) ---
    # One entry per cycle that updated the model.
    # Each entry tagged with space_group for
    # enantiomorph-aware trajectory tracking.
    self.progress = []
    # list of {cycle, program, r_work, r_free,
    #          model_map_cc, space_group, annotation,
    #          timestamp}

    # --- Strategy blacklist (anti-oscillation) ---
    self.strategy_blacklist = []
    # list of {strategy_id, reason, metrics_at_retreat,
    #          cycle, timestamp}

    # --- Hypotheses (Phase 4) ---
    self.hypotheses = []
    # list of Hypothesis objects

    # --- Runtime cache (NOT serialized) ---
    # KD-tree for nearest-residue lookups in
    # _find_diff_peaks (Step 1.1). Cached here so
    # it persists across calls within the same cycle.
    # Invalidated by cache key:
    #   (model_path, mtime, size, space_group, unit_cell)
    self._kdtree_cache = None
    self._kdtree_cache_key = None

  # ── Update methods ──────────────────────────────────

  def update_from_validation(self, validation_result,
                             log_metrics, cycle_number,
                             program_name):
    """Update model state from validation results.

    Args:
      validation_result: dict from run_validation()
        (keys: model_contents, geometry, data_model,
         diff_peaks). May be None.
      log_metrics: dict with r_work, r_free, etc.
        from log analysis. May be None or empty.
      cycle_number: int.
      program_name: str, e.g. "phenix.refine".

    Never raises.
    """
    try:
      self._update_from_validation_inner(
        validation_result, log_metrics,
        cycle_number, program_name,
      )
    except Exception:
      logger.debug(
        "StructureModel.update_from_validation failed",
        exc_info=True,
      )

  def _update_from_validation_inner(
    self, validation_result, log_metrics,
    cycle_number, program_name,
  ):
    """Inner update logic. May raise."""
    metrics = log_metrics if log_metrics else {}
    vr = validation_result if validation_result else {}

    # --- R-factors from log_metrics (authoritative) ---
    r_work = _safe_float(metrics.get("r_work"))
    r_free = _safe_float(metrics.get("r_free"))
    if r_work is not None:
      self.model_state["r_work"] = r_work
    if r_free is not None:
      self.model_state["r_free"] = r_free

    # Model-map CC (cryo-EM)
    cc = _safe_float(metrics.get("model_map_cc"))
    if cc is not None:
      self.model_state["model_map_cc"] = cc

    # --- Model contents from validation ---
    mc = vr.get("model_contents")
    if mc and isinstance(mc, dict):
      # Chains
      chains = mc.get("chains", [])
      if chains:
        self._update_chains(chains, mc)
      # Ligands
      ligands = mc.get("ligands", [])
      if isinstance(ligands, list):
        self._update_ligands(ligands, vr)
      # Waters
      wc = mc.get("waters")
      if wc is not None:
        self.model_state["waters"] = int(wc)
      # Ions
      ions = mc.get("ions", [])
      if isinstance(ions, list):
        self.model_state["ions"] = [
          dict(ion) for ion in ions
        ]

    # --- Geometry from validation ---
    geom = vr.get("geometry")
    if geom and isinstance(geom, dict):
      g = self.model_state["geometry"]
      for key in (
        "rama_favored", "rama_outliers",
        "rama_outlier_list", "rotamer_outliers",
        "rotamer_outlier_list", "clashscore",
        "bonds_rmsd", "angles_rmsd",
      ):
        val = geom.get(key)
        if val is not None:
          g[key] = val

    # --- Difference density from validation ---
    dp = vr.get("diff_peaks")
    if dp and isinstance(dp, dict):
      pos = dp.get("positive", [])
      neg = dp.get("negative", [])
      self.model_state["diff_peaks"] = {
        "positive": list(pos),
        "negative": list(neg),
        "peak_count": dp.get(
          "peak_count", len(pos)
        ),
      }

    # --- Data-model metrics (ligand RSCC) ---
    dm = vr.get("data_model")
    if dm and isinstance(dm, dict):
      rscc_list = dm.get("ligand_rscc", [])
      if rscc_list:
        self._merge_ligand_rscc(rscc_list)

    # --- Chain completeness (Step 1.3) ---
    cc = vr.get("chain_completeness")
    if cc and isinstance(cc, list):
      self.update_chain_completeness(cc)

    # --- Detect problems ---
    self._detect_problems(metrics)

    # --- Record progress entry ---
    sg = self.data_characteristics.get("space_group")
    entry = {
      "cycle": cycle_number,
      "program": program_name,
      "r_work": r_work,
      "r_free": r_free,
      "model_map_cc": cc,
      "space_group": sg,
      "annotation": "",
      "timestamp": time.time(),
    }
    self.progress.append(entry)
    # Cap to 50 entries
    if len(self.progress) > 50:
      self.progress = self.progress[-50:]

  def _update_chains(self, chain_ids, model_contents):
    """Update chain inventory.

    Preserves completeness data from _chain_completeness
    (Step 1.3) if present, otherwise records chain IDs
    with minimal info. Existing per-chain data is kept
    for chains that are still present; removed for chains
    that vanish (e.g. after rebuilding).
    """
    # Build lookup of existing chain data
    existing = {
      c["chain_id"]: c
      for c in self.model_state["chains"]
      if isinstance(c, dict) and "chain_id" in c
    }
    updated = []
    for cid in chain_ids:
      if cid in existing:
        updated.append(existing[cid])
      else:
        updated.append({
          "chain_id": cid,
          "residues_built": None,
          "residues_expected": None,
          "completeness": None,
          "gaps": [],
          "avg_b_factor": None,
          "disordered_fraction": None,
        })
    self.model_state["chains"] = updated

  def _update_ligands(self, ligands_from_validation, vr):
    """Update ligand inventory from validation contents.

    Merges with existing RSCC data if available. New
    ligands are added; ligands no longer in the model
    are removed.
    """
    # Build lookup of existing RSCC data by (name, chain, resid)
    existing_rscc = {}
    for lig in self.model_state.get("ligands", []):
      if isinstance(lig, dict):
        key = (
          lig.get("name"),
          lig.get("chain"),
          lig.get("resid"),
        )
        existing_rscc[key] = {
          "rscc": lig.get("rscc"),
          "z_rscc": lig.get("z_rscc"),
          "b_mean": lig.get("b_mean"),
          "occupancy": lig.get("occupancy"),
        }

    updated = []
    for lig in ligands_from_validation:
      if not isinstance(lig, dict):
        continue
      name = lig.get("name", "UNK")
      chain = lig.get("chain", "")
      resid = lig.get("resid", 0)
      key = (name, chain, resid)
      prev = existing_rscc.get(key, {})
      updated.append({
        "name": name,
        "chain": chain,
        "resid": resid,
        "n_atoms": lig.get("n_atoms", 0),
        "rscc": prev.get("rscc"),
        "z_rscc": prev.get("z_rscc"),
        "b_mean": prev.get("b_mean"),
        "occupancy": prev.get("occupancy"),
      })
    self.model_state["ligands"] = updated

  def _merge_ligand_rscc(self, rscc_list):
    """Merge per-ligand RSCC from data_model validation.

    Matches by (name, chain, resid) against existing
    ligand inventory.
    """
    # Build lookup by (name, chain, resid)
    rscc_lookup = {}
    for entry in rscc_list:
      if not isinstance(entry, dict):
        continue
      key = (
        entry.get("name"),
        entry.get("chain"),
        entry.get("resid"),
      )
      rscc_lookup[key] = entry

    for lig in self.model_state.get("ligands", []):
      key = (
        lig.get("name"),
        lig.get("chain"),
        lig.get("resid"),
      )
      if key in rscc_lookup:
        entry = rscc_lookup[key]
        if entry.get("rscc") is not None:
          lig["rscc"] = _safe_float(entry["rscc"])
        if entry.get("z_rscc") is not None:
          lig["z_rscc"] = _safe_float(entry["z_rscc"])
        if entry.get("b_mean") is not None:
          lig["b_mean"] = _safe_float(entry["b_mean"])
        if entry.get("occupancy") is not None:
          lig["occupancy"] = _safe_float(
            entry["occupancy"]
          )

  def _detect_problems(self, metrics):
    """Identify current structural problems.

    Rebuilds the problems list from scratch each cycle
    based on current model_state. Problems are transient
    observations, not accumulated history — a problem
    that disappears after rebuilding should no longer
    appear.
    """
    problems = []
    geom = self.model_state.get("geometry", {})

    # High R-free gap (overfitting)
    r_work = _safe_float(self.model_state.get("r_work"))
    r_free = _safe_float(self.model_state.get("r_free"))
    if r_work is not None and r_free is not None:
      gap = r_free - r_work
      if gap > 0.10:
        problems.append({
          "region": "global",
          "type": "r_free_gap",
          "severity": "high",
          "description": (
            "R-free gap %.3f (R-work=%.3f, R-free=%.3f)"
            " suggests overfitting"
            % (gap, r_work, r_free)
          ),
        })

    # R-free stalled (3+ cycles without improvement)
    # Note: _detect_problems runs before the current
    # progress entry is appended. So we build the full
    # series from progress + current model_state r_free.
    rfree_vals = [
      _safe_float(p.get("r_free"))
      for p in self.progress
      if p.get("r_free") is not None
    ]
    rfree_vals = [v for v in rfree_vals if v is not None]
    current_rf = _safe_float(
      self.model_state.get("r_free"))
    if current_rf is not None:
      rfree_vals.append(current_rf)
    if len(rfree_vals) >= 3:
      recent = rfree_vals[-3:]
      if all(
        abs(r - recent[0]) < 0.002
        for r in recent[1:]
      ):
        problems.append({
          "region": "global",
          "type": "r_free_stalled",
          "severity": "high",
          "description": (
            "R-free stalled at %.3f for 3+ cycles"
            % recent[-1]
          ),
        })

    # R-free regression (worsened from last cycle)
    if len(rfree_vals) >= 2:
      if rfree_vals[-1] > rfree_vals[-2] + 0.01:
        problems.append({
          "region": "global",
          "type": "r_free_regression",
          "severity": "medium",
          "description": (
            "R-free worsened (%.3f -> %.3f)"
            % (rfree_vals[-2], rfree_vals[-1])
          ),
        })

    # Twinning (critical — affects strategy)
    tw = self.data_characteristics.get(
      "twinning", {}
    )
    if tw.get("is_twinned"):
      tf = tw.get("twin_fraction")
      desc = "Data is twinned"
      if tf is not None:
        desc += " (fraction=%.2f)" % tf
      problems.append({
        "region": "global",
        "type": "twinning",
        "severity": "high",
        "description": desc,
      })

    # High clashscore
    cs = geom.get("clashscore")
    if cs is not None and cs > 10:
      sev = "high" if cs > 20 else "medium"
      problems.append({
        "region": "global",
        "type": "clashscore",
        "severity": sev,
        "description": (
          "Clashscore %.1f (target: <5)" % cs
        ),
      })

    # Ramachandran outliers
    rama_out = geom.get("rama_outlier_list", [])
    if rama_out:
      problems.append({
        "region": ", ".join(rama_out[:5]),
        "type": "rama_outlier",
        "severity": (
          "medium" if len(rama_out) <= 3
          else "high"
        ),
        "description": (
          "%d Ramachandran outlier(s)" % len(rama_out)
        ),
      })

    # Poor ligand fit
    for lig in self.model_state.get("ligands", []):
      rscc = lig.get("rscc")
      if rscc is not None and rscc < 0.7:
        problems.append({
          "region": "%s %s/%d" % (
            lig.get("name", "UNK"),
            lig.get("chain", "?"),
            lig.get("resid", 0),
          ),
          "type": "poor_ligand_fit",
          "severity": (
            "high" if rscc < 0.5 else "medium"
          ),
          "description": (
            "Ligand %s RSCC=%.2f (target: >0.7)"
            % (lig.get("name", "UNK"), rscc)
          ),
        })

    # Significant unmodeled density (positive peaks)
    dp = self.model_state.get("diff_peaks", {})
    pos_peaks = dp.get("positive", [])
    high_peaks = [
      p for p in pos_peaks
      if isinstance(p, dict)
      and _safe_float(p.get("height", 0)) is not None
      and _safe_float(p.get("height", 0)) > 5.0
    ]
    if high_peaks:
      # Describe the strongest peak
      best = max(
        high_peaks,
        key=lambda p: _safe_float(
          p.get("height", 0)
        ),
      )
      near = best.get("near_residue", "unknown")
      ht = _safe_float(best.get("height", 0))
      problems.append({
        "region": "near %s" % near,
        "type": "unmodeled_density",
        "severity": (
          "high" if ht > 8.0 else "medium"
        ),
        "description": (
          "%d peak(s) >5σ in difference density"
          " (strongest: %.1fσ near %s)"
          % (len(high_peaks), ht, near)
        ),
      })

    # Wrongly-placed atoms (negative peaks)
    neg_peaks = dp.get("negative", [])
    deep_holes = [
      p for p in neg_peaks
      if isinstance(p, dict)
      and _safe_float(p.get("height", 0)) is not None
      and abs(
        _safe_float(p.get("height", 0))
      ) > 4.0
    ]
    if deep_holes:
      worst = min(
        deep_holes,
        key=lambda p: _safe_float(
          p.get("height", 0)
        ),
      )
      near = worst.get("near_residue", "unknown")
      ht = _safe_float(worst.get("height", 0))
      problems.append({
        "region": "near %s" % near,
        "type": "misplaced_atoms",
        "severity": (
          "high" if abs(ht) > 6.0 else "medium"
        ),
        "description": (
          "%d negative peak(s) <-4σ "
          "(deepest: %.1fσ near %s) — "
          "atoms may be wrongly placed"
          % (len(deep_holes), ht, near)
        ),
      })

    # Disordered chains
    for ch in self.model_state.get("chains", []):
      if not isinstance(ch, dict):
        continue
      dis = ch.get("disordered_fraction")
      if dis is not None and dis > 0.3:
        problems.append({
          "region": "chain %s" % ch.get(
            "chain_id", "?"
          ),
          "type": "disorder",
          "severity": (
            "high" if dis > 0.5 else "medium"
          ),
          "description": (
            "Chain %s %.0f%% disordered"
            % (ch.get("chain_id", "?"), dis * 100)
          ),
        })

    # Sort by severity (high first)
    severity_order = {"high": 0, "medium": 1, "low": 2}
    problems.sort(
      key=lambda p: severity_order.get(
        p.get("severity", "low"), 2
      )
    )
    self.model_state["problems"] = problems

  def update_from_xtriage(self, xtriage_results):
    """Update data characteristics from xtriage.

    Args:
      xtriage_results: dict extracted from xtriage log.
        Expected keys (all optional):
          resolution, completeness, redundancy,
          space_group, unit_cell,
          is_twinned, twin_law, twin_fraction,
          has_anomalous_data, anomalous_d_min,
          anomalous_measurability,
          i_over_sigma, data_quality.

    Never raises.
    """
    try:
      self._update_from_xtriage_inner(xtriage_results)
    except Exception:
      logger.debug(
        "StructureModel.update_from_xtriage failed",
        exc_info=True,
      )

  def _update_from_xtriage_inner(self, xr):
    """Inner xtriage update. May raise."""
    if not isinstance(xr, dict):
      return
    dc = self.data_characteristics

    # Scalar data characteristics
    for key in (
      "resolution", "completeness", "redundancy",
      "i_over_sigma",
    ):
      val = _safe_float(xr.get(key))
      if val is not None:
        dc[key] = val

    # Space group and unit cell
    sg = xr.get("space_group")
    if sg:
      dc["space_group"] = str(sg)
    uc = xr.get("unit_cell")
    if uc and isinstance(uc, (list, tuple)):
      dc["unit_cell"] = tuple(
        float(x) for x in uc
      )

    # Twinning
    tw = xr.get("is_twinned")
    if tw is not None:
      dc["twinning"]["is_twinned"] = bool(tw)
    tl = xr.get("twin_law")
    if tl is not None:
      dc["twinning"]["twin_law"] = str(tl)
    tf = _safe_float(xr.get("twin_fraction"))
    if tf is not None:
      dc["twinning"]["twin_fraction"] = tf
    # Infer is_twinned from twin_fraction if not
    # explicitly provided. The _extract_xtriage helper
    # provides twin_fraction but not is_twinned, so
    # without this inference, twinned data from history
    # won't be flagged.
    if dc["twinning"]["is_twinned"] is None and tf is not None:
      dc["twinning"]["is_twinned"] = tf > 0.05

    # Anomalous signal — gates hypothesis engine
    anom = xr.get("has_anomalous_data")
    if anom is not None:
      dc["anomalous"]["has_anomalous_data"] = bool(anom)
    ad = _safe_float(xr.get("anomalous_d_min"))
    if ad is not None:
      dc["anomalous"]["anomalous_d_min"] = ad
    am = _safe_float(xr.get("anomalous_measurability"))
    if am is not None:
      dc["anomalous"]["anomalous_measurability"] = am
    # Infer has_anomalous_data from measurability if
    # not explicitly provided. Same rationale as the
    # twinning inference above.
    if (dc["anomalous"]["has_anomalous_data"] is None
        and am is not None):
      dc["anomalous"]["has_anomalous_data"] = (
        am > 0.05
      )

    # Overall quality assessment
    dq = xr.get("data_quality")
    if dq:
      dc["data_quality"] = str(dq)

    # Experiment type
    et = xr.get("experiment_type")
    if et:
      dc["experiment_type"] = str(et)

  def update_from_phaser(self, phaser_results):
    """Update with MR solution quality.

    Args:
      phaser_results: dict with:
        tfz: float (translation function Z-score)
        llg: float (log-likelihood gain)
        space_group: str (confirmed by Phaser)
        unit_cell: tuple of 6 floats
        n_copies: int (copies in ASU)
        Any key may be absent.

    Never raises.
    """
    try:
      self._update_from_phaser_inner(phaser_results)
    except Exception:
      logger.debug(
        "StructureModel.update_from_phaser failed",
        exc_info=True,
      )

  def _update_from_phaser_inner(self, pr):
    """Inner phaser update. May raise."""
    if not isinstance(pr, dict):
      return
    dc = self.data_characteristics

    # Phaser may confirm or change the space group
    sg = pr.get("space_group")
    if sg:
      dc["space_group"] = str(sg)
    uc = pr.get("unit_cell")
    if uc and isinstance(uc, (list, tuple)):
      dc["unit_cell"] = tuple(
        float(x) for x in uc
      )

    # Store MR quality in data_characteristics
    # (these are session-level facts, not per-cycle)
    tfz = _safe_float(pr.get("tfz"))
    if tfz is not None:
      dc["mr_tfz"] = tfz
    llg = _safe_float(pr.get("llg"))
    if llg is not None:
      dc["mr_llg"] = llg
    nc = pr.get("n_copies")
    if nc is not None:
      dc["n_copies_asu"] = int(nc)

  def update_from_autosol(self, autosol_results):
    """Update with experimental phasing results.

    Args:
      autosol_results: dict with:
        fom: float (figure of merit)
        bayes_cc: float (Bayesian CC)
        sites_found: int (number of heavy atom sites)
        Any key may be absent.

    Never raises.
    """
    try:
      self._update_from_autosol_inner(
        autosol_results)
    except Exception:
      logger.debug(
        "StructureModel.update_from_autosol failed",
        exc_info=True,
      )

  def _update_from_autosol_inner(self, ar):
    """Inner autosol update. May raise."""
    if not isinstance(ar, dict):
      return
    dc = self.data_characteristics

    fom = _safe_float(ar.get("fom"))
    if fom is not None:
      dc["phasing_fom"] = fom
    bcc = _safe_float(ar.get("bayes_cc"))
    if bcc is not None:
      dc["phasing_bayes_cc"] = bcc
    sites = ar.get("sites_found")
    if sites is not None:
      try:
        dc["sites_found"] = int(sites)
      except (ValueError, TypeError):
        pass
    # Mark phasing method
    dc["phasing_method"] = "SAD"

  def update_chain_completeness(self, chain_data):
    """Update per-chain completeness from Step 1.3.

    Args:
      chain_data: list of dicts from
        _chain_completeness(), each with:
          chain_id, residues_built,
          residues_expected, completeness_fraction,
          gaps, avg_b_factor, disordered_fraction.

    Never raises.
    """
    try:
      if not isinstance(chain_data, list):
        return
      for cd in chain_data:
        if not isinstance(cd, dict):
          continue
        cid = cd.get("chain_id")
        if not cid:
          continue
        # Find matching chain in model_state
        for ch in self.model_state["chains"]:
          if ch.get("chain_id") == cid:
            for key in (
              "residues_built", "residues_expected",
              "gaps", "avg_b_factor",
              "disordered_fraction",
            ):
              val = cd.get(key)
              if val is not None:
                ch[key] = val
            comp = cd.get("completeness_fraction")
            if comp is not None:
              ch["completeness"] = comp
            break
        else:
          # Chain not yet in inventory — add it
          self.model_state["chains"].append({
            "chain_id": cid,
            "residues_built": cd.get(
              "residues_built"
            ),
            "residues_expected": cd.get(
              "residues_expected"
            ),
            "completeness": cd.get(
              "completeness_fraction"
            ),
            "gaps": cd.get("gaps", []),
            "avg_b_factor": cd.get("avg_b_factor"),
            "disordered_fraction": cd.get(
              "disordered_fraction"
            ),
          })
    except Exception:
      logger.debug(
        "StructureModel.update_chain_completeness"
        " failed",
        exc_info=True,
      )

  def annotate_last_progress(self, annotation):
    """Add annotation to the most recent progress entry.

    Used to record what changed in a cycle, e.g.
    "added 187 waters", "rebuilt chain B loop 145-160".

    Never raises.
    """
    try:
      if self.progress and annotation:
        self.progress[-1]["annotation"] = str(
          annotation
        )[:200]
    except Exception:
      pass

  # ── Query methods ───────────────────────────────────

  def get_summary(self, detail_level="normal"):
    """Produce text summary for prompts or display.

    Args:
      detail_level:
        "brief": 2-3 lines for cycle-level display
        "normal": ~500 chars for THINK prompt
        "detailed": full report for session end

    Returns:
      str. Empty string if no data yet.

    Never raises — returns "" on error.
    """
    try:
      if detail_level == "brief":
        return self._summary_brief()
      elif detail_level == "detailed":
        return self._summary_detailed()
      else:
        return self._summary_normal()
    except Exception:
      logger.debug(
        "StructureModel.get_summary failed",
        exc_info=True,
      )
      return ""

  def _summary_brief(self):
    """2-3 line summary for cycle-level display."""
    parts = []

    # Chains
    chains = self.model_state.get("chains", [])
    if chains:
      chain_strs = []
      for ch in chains:
        cid = ch.get("chain_id", "?")
        comp = ch.get("completeness")
        if comp is not None:
          chain_strs.append(
            "%s: %.0f%%" % (cid, comp * 100)
          )
        else:
          chain_strs.append(cid)
      parts.append(
        "%d chain(s) (%s)"
        % (len(chains), ", ".join(chain_strs))
      )

    # Ligands, ions, waters on one line
    contents = []
    ligs = self.model_state.get("ligands", [])
    if ligs:
      lig_names = []
      for lig in ligs[:3]:
        name = lig.get("name", "UNK")
        rscc = lig.get("rscc")
        if rscc is not None:
          lig_names.append(
            "%s (RSCC=%.2f)" % (name, rscc)
          )
        else:
          lig_names.append(name)
      contents.append(
        "Ligands: %s" % ", ".join(lig_names)
      )
    ions = self.model_state.get("ions", [])
    if ions:
      ion_names = [
        i.get("name", "?") for i in ions[:3]
      ]
      contents.append(
        "%d ion(s) (%s)"
        % (len(ions), ", ".join(ion_names))
      )
    wc = self.model_state.get("waters", 0)
    if wc:
      contents.append("%d waters" % wc)
    if contents:
      parts.append("; ".join(contents))

    # Problems (one line)
    problems = self.model_state.get("problems", [])
    if problems:
      p_strs = [
        p.get("description", "")
        for p in problems[:2]
      ]
      parts.append(
        "Problems: %s" % "; ".join(p_strs)
      )

    return "\n".join(parts) if parts else ""

  def _summary_normal(self):
    """~500 char summary for THINK prompt."""
    parts = []

    # Data line
    dc = self.data_characteristics
    data_parts = []
    res = dc.get("resolution")
    if res is not None:
      data_parts.append("%.1fÅ" % res)
    sg = dc.get("space_group")
    if sg:
      data_parts.append(sg)
    tw = dc.get("twinning", {}).get("is_twinned")
    if tw:
      data_parts.append("twinned")
    if data_parts:
      parts.append("Data: %s" % ", ".join(data_parts))

    # R-factors
    r_work = _safe_float(self.model_state.get("r_work"))
    r_free = _safe_float(self.model_state.get("r_free"))
    if r_work is not None and r_free is not None:
      parts.append(
        "R-work=%.3f R-free=%.3f" % (r_work, r_free)
      )
    elif r_free is not None:
      parts.append("R-free=%.3f" % r_free)

    # Model-map CC (cryo-EM)
    cc = _safe_float(self.model_state.get("model_map_cc"))
    if cc is not None:
      parts.append("Model-map CC=%.3f" % cc)

    # Chains
    chains = self.model_state.get("chains", [])
    if chains:
      chain_strs = []
      for ch in chains:
        cid = ch.get("chain_id", "?")
        comp = ch.get("completeness")
        if comp is not None:
          chain_strs.append(
            "%s: %.0f%% complete" % (cid, comp * 100)
          )
        else:
          chain_strs.append(cid)
      parts.append(
        "Chains: %s" % ", ".join(chain_strs)
      )

    # Ligands
    ligs = self.model_state.get("ligands", [])
    if ligs:
      lig_strs = []
      for lig in ligs[:4]:
        name = lig.get("name", "UNK")
        chain = lig.get("chain", "")
        resid = lig.get("resid", 0)
        rscc = lig.get("rscc")
        if rscc is not None:
          lig_strs.append(
            "%s (%s/%d, RSCC=%.2f)"
            % (name, chain, resid, rscc)
          )
        else:
          lig_strs.append(
            "%s (%s/%d)" % (name, chain, resid)
          )
      parts.append(
        "Ligands: %s" % ", ".join(lig_strs)
      )

    # Ions
    ions = self.model_state.get("ions", [])
    if ions:
      ion_strs = [
        "%s (%s/%d)" % (
          i.get("name", "?"),
          i.get("chain", "?"),
          i.get("resid", 0),
        )
        for i in ions[:4]
      ]
      parts.append("Ions: %s" % ", ".join(ion_strs))

    # Waters
    wc = self.model_state.get("waters", 0)
    if wc:
      parts.append("Waters: %d" % wc)

    # Geometry (compact)
    geom = self.model_state.get("geometry", {})
    geom_parts = []
    rama_fav = geom.get("rama_favored")
    if rama_fav is not None:
      geom_parts.append(
        "Rama=%.1f%%" % (rama_fav * 100)
      )
    cs = geom.get("clashscore")
    if cs is not None:
      geom_parts.append("Clash=%.1f" % cs)
    if geom_parts:
      parts.append(
        "Geometry: %s" % ", ".join(geom_parts)
      )

    # Diff peaks (compact)
    dp = self.model_state.get("diff_peaks", {})
    pc = dp.get("peak_count", 0)
    if pc > 0:
      parts.append(
        "Diff density: %d peak(s) >3.5σ" % pc
      )

    # Problems
    problems = self.model_state.get("problems", [])
    if problems:
      p_strs = [
        p.get("description", "")
        for p in problems[:3]
      ]
      parts.append(
        "Problems: %s" % "; ".join(p_strs)
      )

    # R-free trend (last 5)
    rfree_trend = self.get_rfree_trend()
    if len(rfree_trend) > 1:
      trend_str = " → ".join(
        "%.3f" % v for _, v in rfree_trend[-5:]
      )
      parts.append("R-free trend: %s" % trend_str)

    return "\n".join(parts) if parts else ""

  def _summary_detailed(self):
    """Full report for session end."""
    sections = []

    # Data characteristics
    dc = self.data_characteristics
    dc_parts = []
    res = dc.get("resolution")
    if res is not None:
      dc_parts.append("Resolution: %.2fÅ" % res)
    sg = dc.get("space_group")
    if sg:
      uc = dc.get("unit_cell")
      if uc:
        uc_str = ", ".join("%.1f" % x for x in uc)
        dc_parts.append(
          "Space group: %s (%s)" % (sg, uc_str)
        )
      else:
        dc_parts.append("Space group: %s" % sg)
    comp = dc.get("completeness")
    if comp is not None:
      dc_parts.append(
        "Completeness: %.1f%%" % (comp * 100)
      )
    tw = dc.get("twinning", {})
    if tw.get("is_twinned"):
      tl = tw.get("twin_law", "unknown")
      tf = tw.get("twin_fraction")
      if tf is not None:
        dc_parts.append(
          "Twinned: %s (fraction=%.2f)"
          % (tl, tf)
        )
      else:
        dc_parts.append("Twinned: %s" % tl)
    elif tw.get("is_twinned") is False:
      dc_parts.append("Not twinned")
    anom = dc.get("anomalous", {})
    if anom.get("has_anomalous_data"):
      ad = anom.get("anomalous_d_min")
      if ad is not None:
        dc_parts.append(
          "Anomalous signal to %.1fÅ" % ad
        )
      else:
        dc_parts.append("Anomalous signal present")
    elif anom.get("has_anomalous_data") is False:
      dc_parts.append("No anomalous signal")
    if dc_parts:
      sections.append(
        "Data:\n  %s" % "\n  ".join(dc_parts)
      )

    # R-factors
    r_work = _safe_float(self.model_state.get("r_work"))
    r_free = _safe_float(self.model_state.get("r_free"))
    if r_work is not None and r_free is not None:
      sections.append(
        "R-factors: R-work=%.3f R-free=%.3f"
        % (r_work, r_free)
      )

    # Model contents
    mc_parts = []
    chains = self.model_state.get("chains", [])
    for ch in chains:
      cid = ch.get("chain_id", "?")
      built = ch.get("residues_built")
      exp = ch.get("residues_expected")
      comp = ch.get("completeness")
      gaps = ch.get("gaps", [])
      parts_ch = ["Chain %s:" % cid]
      if built is not None and exp is not None:
        parts_ch.append(
          "residues %d/%d" % (built, exp)
        )
      if comp is not None:
        parts_ch.append("(%.0f%% complete)" % (
          comp * 100
        ))
      if gaps:
        gap_strs = [
          "%d-%d" % (g.get("start", 0),
                     g.get("end", 0))
          for g in gaps[:3]
        ]
        parts_ch.append(
          "gaps: %s" % ", ".join(gap_strs)
        )
      mc_parts.append(" ".join(parts_ch))
    ligs = self.model_state.get("ligands", [])
    for lig in ligs:
      name = lig.get("name", "UNK")
      chain = lig.get("chain", "")
      resid = lig.get("resid", 0)
      rscc = lig.get("rscc")
      bm = lig.get("b_mean")
      lig_parts = [
        "%s (%s/%d)" % (name, chain, resid)
      ]
      if rscc is not None:
        lig_parts.append("RSCC=%.2f" % rscc)
      if bm is not None:
        lig_parts.append("B=%.1f" % bm)
      mc_parts.append(
        "Ligand: %s" % ", ".join(lig_parts)
      )
    ions = self.model_state.get("ions", [])
    for ion in ions:
      mc_parts.append(
        "Ion: %s (%s/%d)" % (
          ion.get("name", "?"),
          ion.get("chain", "?"),
          ion.get("resid", 0),
        )
      )
    wc = self.model_state.get("waters", 0)
    mc_parts.append("Waters: %d" % wc)
    if mc_parts:
      sections.append(
        "Model contents:\n  %s"
        % "\n  ".join(mc_parts)
      )

    # Geometry
    geom = self.model_state.get("geometry", {})
    geom_parts = []
    rama_fav = geom.get("rama_favored")
    if rama_fav is not None:
      geom_parts.append(
        "Ramachandran favored: %.1f%%"
        % (rama_fav * 100)
      )
    rama_out = geom.get("rama_outliers")
    if rama_out is not None:
      geom_parts.append(
        "Ramachandran outliers: %.1f%%"
        % (rama_out * 100)
      )
    rot_out = geom.get("rotamer_outliers")
    if rot_out is not None:
      geom_parts.append(
        "Rotamer outliers: %.1f%%"
        % (rot_out * 100)
      )
    cs = geom.get("clashscore")
    if cs is not None:
      geom_parts.append("Clashscore: %.1f" % cs)
    bonds = geom.get("bonds_rmsd")
    if bonds is not None:
      geom_parts.append(
        "Bonds RMSD: %.4f" % bonds
      )
    angles = geom.get("angles_rmsd")
    if angles is not None:
      geom_parts.append(
        "Angles RMSD: %.2f" % angles
      )
    if geom_parts:
      sections.append(
        "Geometry:\n  %s"
        % "\n  ".join(geom_parts)
      )

    # Diff peaks
    dp = self.model_state.get("diff_peaks", {})
    pos = dp.get("positive", [])
    neg = dp.get("negative", [])
    if pos or neg:
      dp_parts = []
      if pos:
        dp_parts.append(
          "%d positive peak(s)" % len(pos)
        )
        for p in pos[:3]:
          ht = p.get("height", 0)
          near = p.get("near_residue", "?")
          dp_parts.append(
            "  %.1fσ near %s" % (ht, near)
          )
      if neg:
        dp_parts.append(
          "%d negative peak(s)" % len(neg)
        )
      sections.append(
        "Difference density:\n  %s"
        % "\n  ".join(dp_parts)
      )

    # Problems
    problems = self.model_state.get("problems", [])
    if problems:
      p_strs = [
        "[%s] %s" % (
          p.get("severity", "?"),
          p.get("description", ""),
        )
        for p in problems
      ]
      sections.append(
        "Problems:\n  %s" % "\n  ".join(p_strs)
      )

    # Progress history
    if self.progress:
      prog_strs = []
      for entry in self.progress:
        cycle = entry.get("cycle", "?")
        prog = entry.get("program", "?")
        rf = entry.get("r_free")
        ann = entry.get("annotation", "")
        line = "Cycle %s (%s)" % (cycle, prog)
        if rf is not None:
          line += ": R-free=%.3f" % rf
        if ann:
          line += " [%s]" % ann
        prog_strs.append(line)
      sections.append(
        "Progress:\n  %s"
        % "\n  ".join(prog_strs)
      )

    # Strategy blacklist
    if self.strategy_blacklist:
      bl_strs = [
        "%s — %s" % (
          b.get("strategy_id", "?"),
          b.get("reason", ""),
        )
        for b in self.strategy_blacklist
      ]
      sections.append(
        "Blacklisted strategies:\n  %s"
        % "\n  ".join(bl_strs)
      )

    # Hypotheses
    hyps = [h for h in self.hypotheses if h.is_resolved]
    if hyps:
      h_strs = []
      for h in hyps:
        symbol = (
          "✓" if h.status == "confirmed"
          else "✗" if h.status == "refuted"
          else "—"
        )
        h_strs.append(
          '%s "%s" (%s, cycle %s)'
          % (symbol, h.statement, h.status,
             h.resolved_at_cycle)
        )
      sections.append(
        "Hypotheses tested:\n  %s"
        % "\n  ".join(h_strs)
      )

    return "\n\n".join(sections) if sections else ""

  def get_current_problems(self):
    """List unresolved problems, ordered by severity.

    Returns:
      list of {problem, severity, suggested_action}.
      Empty list if no problems.

    Never raises.
    """
    try:
      problems = self.model_state.get("problems", [])
      result = []
      for p in problems:
        suggested = _suggest_action(p)
        result.append({
          "problem": p.get("description", ""),
          "severity": p.get("severity", "low"),
          "suggested_action": suggested,
        })
      return result
    except Exception:
      logger.debug(
        "StructureModel.get_current_problems failed",
        exc_info=True,
      )
      return []

  def get_rfree_trend(self, space_group=None):
    """Extract R-free time series from progress.

    Args:
      space_group: str or None. If provided, return
        only entries matching this space group. Used
        for enantiomorph-aware trajectory comparison.

    Returns:
      list of (cycle_number, r_free) tuples.
    """
    try:
      result = []
      for entry in self.progress:
        rf = entry.get("r_free")
        if rf is None:
          continue
        if space_group is not None:
          if entry.get("space_group") != space_group:
            continue
        result.append(
          (entry.get("cycle", 0), rf)
        )
      return result
    except Exception:
      return []

  def get_metric(self, metric_name):
    """Get the current value of a named metric.

    Used by the gate evaluator (Phase 3) to check
    success criteria like {"r_free": "<0.35"}.

    Supports:
      Direct: r_free, r_work, model_map_cc
      Geometry: clashscore, rama_favored, rama_outliers,
        rotamer_outliers, bonds_rmsd, angles_rmsd
      Data: resolution, tfz, llg
      Counts: waters, water_count, ligand_count,
        ion_count, positive_peak_count,
        negative_peak_count
      Derived: r_free_gap, ligand_cc, ligand_rscc,
        ligand_z_rscc

    Returns:
      float, int, or None.
    """
    try:
      ms = self.model_state
      dc = self.data_characteristics

      # Direct model_state fields (float)
      if metric_name in ("r_free", "r_work",
                         "model_map_cc"):
        val = ms.get(metric_name)
        return _safe_float(val) if val is not None \
          else None

      # Waters — return as int
      if metric_name in ("waters", "water_count"):
        val = ms.get("waters", 0)
        return int(val) if val is not None else 0

      # Geometry fields
      geom = ms.get("geometry", {})
      if metric_name in geom:
        return _safe_float(geom.get(metric_name))

      # Data characteristics
      if metric_name == "resolution":
        return _safe_float(dc.get("resolution"))
      if metric_name == "tfz":
        return _safe_float(dc.get("mr_tfz"))
      if metric_name == "llg":
        return _safe_float(dc.get("mr_llg"))

      # Derived: R-free gap (overfitting indicator)
      if metric_name == "r_free_gap":
        rw = _safe_float(ms.get("r_work"))
        rf = _safe_float(ms.get("r_free"))
        if rw is not None and rf is not None:
          return rf - rw
        return None

      # Derived: worst (minimum) ligand RSCC
      # Plan uses "ligand_rscc"; alias "ligand_cc"
      if metric_name in ("ligand_cc", "ligand_rscc"):
        ligs = ms.get("ligands", [])
        rsccs = [
          _safe_float(l.get("rscc"))
          for l in ligs
          if l.get("rscc") is not None
        ]
        return min(rsccs) if rsccs else None

      # Derived: worst (minimum) ligand Z-score
      if metric_name == "ligand_z_rscc":
        ligs = ms.get("ligands", [])
        zs = [
          _safe_float(l.get("z_rscc"))
          for l in ligs
          if l.get("z_rscc") is not None
        ]
        return min(zs) if zs else None

      # Count metrics
      if metric_name == "ligand_count":
        return len(ms.get("ligands", []))
      if metric_name == "ion_count":
        return len(ms.get("ions", []))
      if metric_name == "positive_peak_count":
        dp = ms.get("diff_peaks", {})
        return len(dp.get("positive", []))
      if metric_name == "negative_peak_count":
        dp = ms.get("diff_peaks", {})
        return len(dp.get("negative", []))

      # Anomalous data characteristics (v115.05)
      anom = dc.get("anomalous", {})
      if metric_name == "anomalous_measurability":
        return _safe_float(
          anom.get("anomalous_measurability"))
      if metric_name == "has_anomalous":
        return anom.get("has_anomalous_data")

      return None
    except Exception:
      return None

  # ── Strategy blacklist ──────────────────────────────

  def blacklist_strategy(self, strategy_id, reason,
                         metrics_at_retreat=None):
    """Record a failed strategy to prevent oscillation.

    Args:
      strategy_id: str, unique identifier for the
        strategy (e.g. "MR_with_beta.pdb",
        "refine_P3121").
      reason: str, why it failed.
      metrics_at_retreat: dict, e.g.
        {"r_free": 0.482, "cycle": 5}.

    Never raises.
    """
    try:
      self.strategy_blacklist.append({
        "strategy_id": str(strategy_id),
        "reason": str(reason)[:200],
        "metrics_at_retreat": (
          dict(metrics_at_retreat)
          if metrics_at_retreat else {}
        ),
        "cycle": (
          self.progress[-1].get("cycle", 0)
          if self.progress else 0
        ),
        "timestamp": time.time(),
      })
      # Cap to 20 entries
      if len(self.strategy_blacklist) > 20:
        self.strategy_blacklist = (
          self.strategy_blacklist[-20:]
        )
    except Exception:
      logger.debug(
        "StructureModel.blacklist_strategy failed",
        exc_info=True,
      )

  def is_blacklisted(self, strategy_id):
    """Check if a strategy was previously tried and failed.

    Returns:
      bool
    """
    try:
      sid = str(strategy_id)
      return any(
        b.get("strategy_id") == sid
        for b in self.strategy_blacklist
      )
    except Exception:
      return False

  def get_blacklist_reason(self, strategy_id):
    """Get the reason a strategy was blacklisted.

    Returns:
      str or None.
    """
    try:
      sid = str(strategy_id)
      for b in self.strategy_blacklist:
        if b.get("strategy_id") == sid:
          return b.get("reason", "")
      return None
    except Exception:
      return None

  # ── Hypothesis management ───────────────────────────

  def get_active_hypothesis(self):
    """Get the currently active hypothesis, if any.

    Returns:
      Hypothesis or None.
    """
    try:
      for h in self.hypotheses:
        if h.is_active:
          return h
      return None
    except Exception:
      return None

  def has_active_hypothesis(self):
    """Check if there's an active hypothesis.

    Returns:
      bool.
    """
    return self.get_active_hypothesis() is not None

  def add_hypothesis(self, hypothesis):
    """Add a hypothesis, enforcing single-active budget.

    Args:
      hypothesis: Hypothesis instance.

    Returns:
      bool — True if added, False if rejected (because
      an active hypothesis already exists and the new
      one would also be active).

    Never raises.
    """
    try:
      if not isinstance(hypothesis, Hypothesis):
        return False
      # Enforce single-active budget
      if hypothesis.is_active and \
          self.has_active_hypothesis():
        logger.debug(
          "Rejected hypothesis %r: active slot "
          "occupied", hypothesis.id,
        )
        return False
      self.hypotheses.append(hypothesis)
      # Cap total to 30 (including resolved)
      if len(self.hypotheses) > 30:
        # Remove oldest resolved hypotheses first
        resolved = [
          (i, h) for i, h in enumerate(self.hypotheses)
          if h.is_resolved
        ]
        if resolved:
          self.hypotheses.pop(resolved[0][0])
      return True
    except Exception:
      logger.debug(
        "StructureModel.add_hypothesis failed",
        exc_info=True,
      )
      return False

  def get_hypothesis(self, hypothesis_id):
    """Look up a hypothesis by ID.

    Returns:
      Hypothesis or None.
    """
    try:
      for h in self.hypotheses:
        if h.id == hypothesis_id:
          return h
      return None
    except Exception:
      return None

  # ── Serialization ───────────────────────────────────

  def to_dict(self):
    """Serialize to JSON-safe dict for session persistence.

    Returns:
      dict. All values are JSON-serializable.
    """
    try:
      return {
        "data_characteristics": copy.deepcopy(
          self.data_characteristics
        ),
        "model_state": copy.deepcopy(
          self.model_state
        ),
        "progress": [
          dict(p) for p in self.progress
        ],
        "strategy_blacklist": [
          dict(b) for b in self.strategy_blacklist
        ],
        "hypotheses": [
          h.to_dict() for h in self.hypotheses
        ],
        "_version": 1,
      }
    except Exception:
      logger.debug(
        "StructureModel.to_dict failed",
        exc_info=True,
      )
      return {}

  @classmethod
  def from_dict(cls, d):
    """Deserialize from session data.

    Tolerant of missing/extra keys — the session may
    be from an earlier version of the code.

    Args:
      d: dict from to_dict().

    Returns:
      StructureModel instance.
    """
    model = cls()
    if not isinstance(d, dict):
      return model
    try:
      # Data characteristics — merge to preserve
      # defaults for any new keys added later
      dc = d.get("data_characteristics")
      if isinstance(dc, dict):
        _deep_merge(
          model.data_characteristics, dc
        )

      # Model state — merge similarly
      ms = d.get("model_state")
      if isinstance(ms, dict):
        _deep_merge(model.model_state, ms)

      # Coerce all numeric fields after merging raw
      # JSON values.  _deep_merge copies values as-is
      # from the deserialized dict, so strings like
      # "0.385" can end up in numeric fields.
      # _coerce_numerics converts them to float/int
      # or None, preventing TypeError in arithmetic
      # (e.g. _detect_problems gap = r_free - r_work).
      _coerce_numerics(model.model_state)

      # Progress
      prog = d.get("progress")
      if isinstance(prog, list):
        model.progress = [
          dict(p) for p in prog
          if isinstance(p, dict)
        ]
        # Coerce progress entry numerics too
        for p in model.progress:
          for k in ("r_work", "r_free", "model_map_cc"):
            if k in p and p[k] is not None:
              p[k] = _safe_float(p[k])

      # Strategy blacklist
      bl = d.get("strategy_blacklist")
      if isinstance(bl, list):
        model.strategy_blacklist = [
          dict(b) for b in bl
          if isinstance(b, dict)
        ]

      # Hypotheses
      hyps = d.get("hypotheses")
      if isinstance(hyps, list):
        model.hypotheses = [
          Hypothesis.from_dict(h)
          for h in hyps
          if isinstance(h, dict)
        ]
    except Exception:
      logger.debug(
        "StructureModel.from_dict partial failure",
        exc_info=True,
      )
    return model

  def compute_fingerprint(self):
    """Content hash of key model state.

    Used by the plan generator (Phase 2) to detect
    when the structure model has changed enough to
    warrant plan revision. Covers R-free, ligands
    placed, space group, and problem count.

    Returns:
      str (hex digest).
    """
    try:
      data = {
        "r_free": self.model_state.get("r_free"),
        "r_work": self.model_state.get("r_work"),
        "space_group": (
          self.data_characteristics.get("space_group")
        ),
        "n_ligands": len(
          self.model_state.get("ligands", [])
        ),
        "n_problems": len(
          self.model_state.get("problems", [])
        ),
        "n_blacklisted": len(self.strategy_blacklist),
      }
      raw = json.dumps(
        data, sort_keys=True, default=str
      )
      return hashlib.md5(
        raw.encode("utf-8")
      ).hexdigest()[:12]
    except Exception:
      return ""


# ── Helpers (module-private) ────────────────────────────

def _safe_float(val):
  """Convert to float or return None."""
  if val is None:
    return None
  try:
    return float(val)
  except (ValueError, TypeError):
    return None


def _deep_merge(target, source):
  """Recursively merge source dict into target.

  Overwrites leaf values in target with source values.
  For dicts, recurses. For lists and scalars, replaces.
  """
  for key, val in source.items():
    if (
      key in target
      and isinstance(target[key], dict)
      and isinstance(val, dict)
    ):
      _deep_merge(target[key], val)
    else:
      target[key] = val


# Fields in model_state that must be float (or None).
_MODEL_STATE_FLOAT_FIELDS = (
  "r_work", "r_free", "model_map_cc",
)

# Fields in model_state.geometry that must be float (or None).
_GEOMETRY_FLOAT_FIELDS = (
  "rama_favored", "rama_outliers", "rotamer_outliers",
  "clashscore", "bonds_rmsd", "angles_rmsd",
)


def _coerce_numerics(model_state):
  """Coerce all known numeric fields in model_state.

  After JSON deserialization via _deep_merge, values may
  be strings (e.g. "0.385" instead of 0.385).  This
  causes TypeError when arithmetic is performed
  (e.g. r_free - r_work in _detect_problems).

  Converts to float via _safe_float; invalid values
  become None.  Runs once after from_dict().
  """
  for key in _MODEL_STATE_FLOAT_FIELDS:
    if key in model_state and model_state[key] is not None:
      model_state[key] = _safe_float(model_state[key])

  geom = model_state.get("geometry")
  if isinstance(geom, dict):
    for key in _GEOMETRY_FLOAT_FIELDS:
      if key in geom and geom[key] is not None:
        geom[key] = _safe_float(geom[key])

  # diff_peaks.peak_count → int
  dp = model_state.get("diff_peaks")
  if isinstance(dp, dict) and "peak_count" in dp:
    try:
      dp["peak_count"] = int(dp["peak_count"])
    except (ValueError, TypeError):
      dp["peak_count"] = 0

  # waters → int
  if "waters" in model_state:
    try:
      model_state["waters"] = int(
        model_state["waters"])
    except (ValueError, TypeError):
      model_state["waters"] = 0


def _suggest_action(problem):
  """Suggest an action for a detected problem.

  Returns str. Used by get_current_problems().
  """
  ptype = problem.get("type", "")
  if ptype == "r_free_gap":
    return (
      "Consider fewer refinement parameters "
      "(group B-factors, remove TLS)"
    )
  if ptype == "r_free_stalled":
    return (
      "Try a different strategy: rebuild model, "
      "add/remove features, or check space group"
    )
  if ptype == "r_free_regression":
    return (
      "Check last cycle for errors; consider "
      "reverting to the previous model"
    )
  if ptype == "twinning":
    return (
      "Use twin refinement in phenix.refine "
      "(twin_law from xtriage)"
    )
  if ptype == "clashscore":
    return "Run refinement with geometry optimization"
  if ptype == "rama_outlier":
    return (
      "Inspect outliers in Coot; rebuild if in "
      "poorly-fit regions"
    )
  if ptype == "poor_ligand_fit":
    return (
      "Re-examine ligand placement; consider "
      "alternative ligand identity"
    )
  if ptype == "unmodeled_density":
    return (
      "Inspect peaks in Coot; consider ligand, "
      "ion, or buffer molecule"
    )
  if ptype == "misplaced_atoms":
    return (
      "Inspect negative peaks in Coot; atoms near "
      "these holes may need rebuilding or removal"
    )
  if ptype == "disorder":
    return (
      "Expect higher R-free; do not over-refine "
      "disordered regions"
    )
  return ""
