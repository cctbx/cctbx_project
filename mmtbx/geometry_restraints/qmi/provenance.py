"""Provenance comment block for the qmi ``geometry_restraints.edits`` header
(plan Sec 9.3).

Two blocks, both emitted as ``# ...`` comment lines by
:func:`mmtbx.geometry_restraints.qmi.restraints.as_edits_phil`: a GLOBAL block
(tool / versions / model / protonation, once per file) and a per-SITE block
(committed state, source, confidence, oxidation/spin reasons, advisory BVS,
restraint origin, flags). No JSON sidecar; every emitted restraint is auditable
from the ``.edits`` alone.
"""

from __future__ import absolute_import, division, print_function

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


def site_block(metal_element, site_tag, oxidation, multiplicity, charge,
               state_source, confidence, oxidation_reason, spin_reason,
               bvs_est, sigma_origin, n_skipped_symmetry, n_perturb,
               sanity=None):
  """The per-site provenance lines (committed state + how it was decided).

  *sanity* is the :func:`mmtbx.geometry_restraints.qmi.sanity.check_site` result
  (or None); any atypical / implausible emitted distances are flagged here
  (Sec 9.4, flag-not-fix), with a per-bond detail line."""
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
    lines.append(f"  BVS most-likely ox={bvs_est.best:+d} (advisory)")
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
