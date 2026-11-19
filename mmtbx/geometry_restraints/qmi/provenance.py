"""Provenance comment block for the qmi ``geometry_restraints.edits`` header
(plan Sec 9.3).

Three kinds of lines, all emitted as ``# ...`` comments by
:func:`mmtbx.geometry_restraints.qmi.restraints.as_edits_phil`: a GLOBAL block
(tool / versions / model / protonation, once per file), a SETTINGS block (the
effective region/freeze partition and per-engine optimisation knobs, once per
file), and a per-SITE block (committed state, source, confidence, oxidation/spin
reasons, advisory BVS, restraint origin, flags). No JSON sidecar; every emitted
restraint, and the settings that produced it, are auditable from the ``.edits``
alone.
"""

import libtbx.version
from libtbx import easy_run

from mmtbx.geometry_restraints import xtb_manager


def xtb_version():
  """The local xtb version (parsed from ``xtb --version``), or ``"unknown"``."""
  exe = xtb_manager.get_exe()
  if not exe:
    return "unknown"
  try:
    out = easy_run.fully_buffered(f"{exe} --version")
    for line in out.stdout_lines + out.stderr_lines:
      toks = line.split()
      if "version" in toks:                       # "* xtb version 6.7.1 (...)"
        i = toks.index("version")
        if i + 1 < len(toks):
          return toks[i + 1]
  except Exception:
    pass
  return "unknown"


def cctbx_version():
  """The cctbx/Phenix version string, or ``"unknown"``."""
  try:
    return libtbx.version.get_version()
  except Exception:
    return "unknown"


def global_block(method, solvent, temperature, model_name, protonation,
                 resolution=None, xtb_ver=None, cctbx_ver=None):
  """The global provenance lines (once per file).

  *protonation* is a short description of the protonation source/opt-out state.
  *resolution* is included only when not None (qmi does not require it). Versions
  are queried lazily unless *xtb_ver* / *cctbx_ver* are supplied (for tests)."""
  xtb_ver = xtb_version() if xtb_ver is None else xtb_ver
  cctbx_ver = cctbx_version() if cctbx_ver is None else cctbx_ver
  model_line = f"model: {model_name}"
  if resolution is not None:
    model_line += f"; resolution {resolution:.2f} A"
  return [
    f"qmi metal restraints  ({method} / {solvent or 'gas'} / {temperature:.1f} K)",
    f"versions: xtb {xtb_ver}; cctbx {cctbx_ver}",
    model_line,
    f"protonation: {protonation}",
  ]


def settings_block(engine, radius, free_radius, freeze_waters,
                   freeze_far_carboxylate_o, optimise_all,
                   optimise_hydrogens=None, optimisation=None,
                   hessian=None, basis_set=None, nproc=1):
  """The effective computational settings (once per file), as header lines.

  Records the region/freeze partition and the per-engine optimisation knobs that
  produced the targets, so the ``.edits`` documents how it was generated, not
  only the chemistry outcome. *optimise_hydrogens* / *optimisation* / *hessian* /
  *basis_set* are ORCA-only (xTB frees every non-cap hydrogen). *optimisation* is
  the collapsed convergence knob ('looser' or an ORCA convergence tier)."""
  def onoff(flag):
    return "on" if flag else "off"
  region = (f"region: radius={radius} "
            f"free_radius={free_radius}; freeze_waters={freeze_waters} "
            f"freeze_far_carboxylate_o={onoff(freeze_far_carboxylate_o)} "
            f"optimise_all={onoff(optimise_all)}")
  if engine == "orca":
    eng = f"engine: orca H={optimise_hydrogens} opt={optimisation} hessian={hessian}"
    if basis_set:
      eng += f" basis={basis_set}"
    eng += f"; nproc={nproc}"
  else:
    eng = f"engine: xtb (GFN2) all non-cap H free; nproc={nproc}"
  return [region, eng]


def site_block(metal_element, site_tag, oxidation, multiplicity, charge,
               state_source, confidence, oxidation_reason, spin_reason,
               bvs_est, sigma_origin, n_skipped_symmetry, n_perturb,
               sanity=None, donor_summary=None):
  """The per-site provenance lines (committed state + how it was decided).

  *sanity* is the :func:`mmtbx.geometry_restraints.qmi.sanity.check_site` result
  (or None); any atypical / implausible emitted distances are flagged here
  (Sec 9.4, flag-not-fix), with a per-bond detail line.

  *donor_summary* is the :func:`mmtbx.geometry_restraints.qmi.donors.summarize`
  result (or None); donors the residue-name classifier could not place are
  reported so an empty donor summary is never silent."""
  src = state_source + (f"/{confidence}" if confidence else "")
  lines = [f"site {site_tag}: {metal_element} ox={oxidation:+d} "
           f"mult={multiplicity} charge={charge:+d}  [source={src}]"]
  reasons = []
  if oxidation_reason:
    reasons.append(f"oxidation ({oxidation_reason})")
  if spin_reason:
    reasons.append(f"spin ({spin_reason})")
  if reasons:
    lines.append("  " + "; ".join(reasons))
  if bvs_est is not None:
    if getattr(bvs_est, "coverage_asymmetric", False):
      # Candidate states were summed over different bond counts, so the
      # comparison is meaningless; report the refusal rather than a biased
      # winner (e.g. Mn(3+) has no vendored S parameters).
      lines.append("  BVS: no verdict, parameter coverage differs between "
                   "candidate states (advisory)")
    else:
      margin = getattr(bvs_est, "margin", None)
      note = ""
      if margin is not None:
        note = f", margin {margin:.2f}"
        if getattr(bvs_est, "close_call", False):
          note += " (close call: not stable against re-refinement [84])"
      lines.append(f"  BVS most-likely ox={bvs_est.best:+d}{note} (advisory)")
      # The advisory and the committed state can disagree, and until now the
      # header printed both without ever saying so: a reader had to notice the
      # contradiction unaided. Say it. This does NOT promote BVS to a decider
      # (it is a bond-length function, and the restraint qmi emits moves the
      # very lengths it reads [84]); it makes the disagreement actionable.
      # No redox-lability test is needed: a BVS verdict exists only where two or
      # more oxidation states are accessible, which is that condition.
      if bvs_est.best != oxidation and state_source != "bvs":
        lines.append(f"  BVS DISAGREES with the committed ox={oxidation:+d} "
                     f"(advisory, from observed bond lengths)")
        if bvs_est.close_call:
          lines.append("    the margin is a close call, so the disagreement is "
                       "not itself evidence against the commit [84]")
        if state_source != "user-supplied":
          lines.append(f"    to commit ox={bvs_est.best:+d} instead, add a "
                       f"metal_state block; states=reasonable emits both")
  # Donors the residue-name classifier could not place. Silence here is what let
  # a bis-citrate Fe commit a ferrous default off an empty summary.
  n_unc = getattr(donor_summary, "n_unclassified", 0) if donor_summary else 0
  if n_unc:
    which = ", ".join(getattr(donor_summary, "unclassified", []) or [])
    n_bg = getattr(donor_summary, "n_anionic_unplaced", 0)
    if n_bg:
      # The residue name told us nothing, but the bond graph did: report the
      # charge that actually informed the state.
      tail = (f"; {n_bg} anionic O group(s) recovered from the bond graph "
              f"(q_hard {donor_summary.q_hard:+d})")
    elif not getattr(donor_summary, "motif_evidence", True):
      tail = "; NO donor was classified, the state is a default"
    else:
      tail = ""
    lines.append(f"  unclassified donors: {n_unc} ({which}){tail}")
  flags = [f"sigmas: {sigma_origin}"]
  if n_perturb:
    outcome = "recovered" if sigma_origin == "hessian" else "unresolved"
    flags.append(f"{outcome} after {n_perturb} perturbation(s)")
  if n_skipped_symmetry:
    flags.append(f"{n_skipped_symmetry} symmetry donor(s) skipped")
  n_mode_changed = getattr(sanity, "n_mode_changed", 0) if sanity else 0
  if sanity is not None and (sanity.n_atypical or sanity.n_implausible
                             or n_mode_changed):
    counts = []
    if sanity.n_atypical:
      counts.append(f"{sanity.n_atypical} atypical")
    if sanity.n_implausible:
      counts.append(f"{sanity.n_implausible} implausible")
    if n_mode_changed:
      counts.append(f"{n_mode_changed} mode-changed")
    flags.append("sanity " + " + ".join(counts))
  lines.append("  " + "; ".join(flags))
  if sanity is not None and sanity.details:
    lines.append("  sanity: " + "; ".join(sanity.details))
  return lines
