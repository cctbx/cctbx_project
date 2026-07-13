"""In-pipeline protonation for the qmi metal-restraint Program.

The QM cluster needs hydrogens, and metal-coordinating sidechain donors need
the *right* protonation (reduce2 protonates every sidechain to its neutral
biological-pH form, which is wrong at a metal site). This module ports the
protonation prep from the calibration ``prepare_qm_calcs.py``:

- ``add_hydrogens``: run ``mmtbx.reduce2`` on the model to add H.
- ``deprotonate_metal_donors``: fix the protonation of donors coordinating any
  metal. Cys/Tyr/Asp/Glu carboxylate donors are stripped to their anionic
  form; His is handled tautomer-aware (it stays NEUTRAL, not an imidazolate):
  the coordinating ring N is made H-free while the partner ring N keeps a
  proton.

The Program's validate() guard uses ``mmtbx.model.manager.has_hd()`` to decide
whether the model is already protonated. Water protonation is intentionally
left out for now.
"""

from __future__ import absolute_import, division, print_function

import os
import shutil
import tempfile
from contextlib import contextmanager

from libtbx import easy_run
from libtbx.utils import null_out, Sorry
from scitbx import matrix

from mmtbx.geometry_restraints.qmi import metals as qmi_metals


# Distance defining a metal-coordinating donor (Harding / CheckMyMetal).
DEFAULT_COORD_CUTOFF_A = 3.0

# Anionic donors: (resname, donor heavy-atom) -> H names to strip so charge
# accounting reads them as anionic, not neutral.
_COORD_DEPROT_RULES = {
  ("CYS", "SG"):  {"HG", "DG"},
  ("TYR", "OH"):  {"HH", "DH"},
  ("ASP", "OD1"): {"HD2", "DD2"},
  ("ASP", "OD2"): {"HD2", "DD2"},
  ("GLU", "OE1"): {"HE2", "DE2"},
  ("GLU", "OE2"): {"HE2", "DE2"},
}

# Histidine is tautomer-aware (neutral, not imidazolate): the coordinating ring
# N donates through its lone pair and must be H-free; the partner ring N keeps
# the proton.
#   ring N -> (its two heavy ring neighbours, used to place a relocated N-H)
_HIS_RING_NEIGHBOURS = {
  "ND1": ("CG", "CE1"),
  "NE2": ("CD2", "CE1"),
}
#   ring N -> the H/D names that sit on it
_HIS_RING_H = {
  "ND1": ("HD1", "DD1"),
  "NE2": ("HE2", "DE2"),
}
_NH_BOND_A = 1.01  # placed N-H length; the optimiser frees all H, so it relaxes


@contextmanager
def _in_directory(path):
  saved = os.getcwd()
  os.chdir(path)
  try:
    yield
  finally:
    os.chdir(saved)


def add_hydrogens(model, use_neutron_distances=True, log=None):
  """Add hydrogens to *model* by running ``mmtbx.reduce2`` and return a fresh
  ``mmtbx.model.manager`` for the H-completed structure.

  Neutron X-H distances by default (closer to the QM-optimal positions the
  downstream optimisation relaxes to). Runs in a scratch directory (removed on
  return) so reduce2's side output does not litter the cwd.
  """
  import iotbx.pdb
  import mmtbx.model
  if log is None:
    log = null_out()
  tmp = tempfile.mkdtemp(prefix="qmi_reduce_")
  try:
    in_path = os.path.join(tmp, "in.cif")
    out_path = os.path.join(tmp, "in_H.cif")
    with open(in_path, "w") as fh:
      fh.write(model.get_hierarchy().as_mmcif_string(
        crystal_symmetry=model.crystal_symmetry()))
    cmd = (f"mmtbx.reduce2 {in_path} output.filename={out_path} "
           f"use_neutron_distances={use_neutron_distances} --overwrite")
    print(f"  adding hydrogens: {cmd}", file=log)
    with _in_directory(tmp):
      result = easy_run.fully_buffered(cmd)
    if result.return_code != 0 or not os.path.exists(out_path):
      tail = "\n".join((result.stdout_lines + result.stderr_lines)[-15:])
      raise Sorry(f"mmtbx.reduce2 failed to protonate the model "
                  f"(rc={result.return_code}):\n{tail}")
    # Carry the original cell forward explicitly (an explicit crystal_symmetry=
    # wins over model_input): symmetry-aware endo_exo cannot silently degrade if
    # reduce2's mmCIF ever drops CRYST1. iotbx.pdb.input reads the file fully, so
    # the scratch dir is safe to remove in the finally.
    return mmtbx.model.manager(
      model_input=iotbx.pdb.input(file_name=out_path),
      crystal_symmetry=model.crystal_symmetry())
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------------------
# His tautomer helpers (ported from prepare_qm_calcs.py)
# ---------------------------------------------------------------------------
def _imidazole_nh_xyz(atoms_by_name, n_name, nh_len=_NH_BOND_A):
  """In-plane bisector position for an N-H on the imidazole ring nitrogen
  *n_name*, or None if a ring atom needed for the construction is missing."""
  nbr1, nbr2 = _HIS_RING_NEIGHBOURS[n_name]
  n = atoms_by_name.get(n_name)
  a1 = atoms_by_name.get(nbr1)
  a2 = atoms_by_name.get(nbr2)
  if n is None or a1 is None or a2 is None:
    return None
  npos = matrix.col(n.xyz)
  v1 = (matrix.col(a1.xyz) - npos).normalize()
  v2 = (matrix.col(a2.xyz) - npos).normalize()
  bisector = v1 + v2
  if bisector.length() < 1e-6:
    return None
  return tuple(npos - bisector.normalize() * nh_len)  # outward, in ring plane


def _relocate_or_remove_his_h(ag, coord_n, coord_h, metal_d):
  """Fix a metal-coordinating His whose coordinating ring N (*coord_n*) still
  carries a proton (*coord_h*). Returns a short action string for logging.

  - both ring N protonated (HIP, +1): remove the coordinating-N proton, giving
    a neutral His.
  - only the coordinating N protonated (wrong tautomer): relocate that proton
    to the partner N, which also renames it (HD1<->HE2 / DD1<->DE2 encode the
    tautomer site).
  """
  atoms_by_name = {a.name.strip().upper(): a for a in ag.atoms()}
  other_n = "NE2" if coord_n == "ND1" else "ND1"
  other_has_h = any(h in atoms_by_name for h in _HIS_RING_H[other_n])
  if other_has_h:
    old = coord_h.name.strip()
    ag.remove_atom(coord_h)
    return f"removed {old} from {coord_n} (HIP to neutral, d={metal_d:.2f} A)"
  new_xyz = _imidazole_nh_xyz(atoms_by_name, other_n)
  if new_xyz is None:
    old = coord_h.name.strip()
    ag.remove_atom(coord_h)
    return (f"removed {old} from {coord_n} (no ring geom to relocate, "
            f"d={metal_d:.2f} A)")
  old_raw = coord_h.name
  old = old_raw.strip()
  isotope = "D" if old.upper().startswith("D") else "H"
  new_name = isotope + ("E2" if other_n == "NE2" else "D1")
  # preserve the 4-char field justification (old and new are both 3 chars)
  lead = len(old_raw) - len(old_raw.lstrip())
  coord_h.name = " " * lead + new_name + " " * (len(old_raw) - lead - len(new_name))
  coord_h.set_xyz(new_xyz)
  return (f"relocated {old}({coord_n}) to {new_name}({other_n}) "
          f"wrong tautomer, d={metal_d:.2f} A")


# ---------------------------------------------------------------------------
# metal-donor deprotonation
# ---------------------------------------------------------------------------
def deprotonate_metal_donors(hier, metals=None,
                             coord_cutoff_A=DEFAULT_COORD_CUTOFF_A,
                             log=None):
  """Fix, IN PLACE, the protonation of sidechain donors in *hier* coordinating
  any metal atom (within *coord_cutoff_A*).

  - Cys/Tyr/Asp/Glu (:data:`_COORD_DEPROT_RULES`): strip the donor H so charge
    accounting reads them as anionic.
  - His: tautomer-aware (neutral, not imidazolate) via
    :func:`_relocate_or_remove_his_h`.

  *metals* is the element set the donors are scored against (default: all
  METALS); the in-pipeline caller passes the subset qmi actually models, so
  donors near an out-of-scope metal (e.g. an alkali or alkaline-earth ion) are
  left untouched.

  Canonicalises atom order at the end (a relocated His proton is renamed in
  place). Returns ``(n_donors_corrected, n_his_corrected)``. A no-op (0, 0)
  when the model has no in-scope metals.
  """
  if log is None:
    log = null_out()
  metal_xyzs = [matrix.col(a.xyz)
                for a in qmi_metals.metal_atoms(hier.atoms(), metals=metals)]
  if not metal_xyzs:
    return 0, 0

  # Cys / Tyr / Asp / Glu: strip the donor proton.
  to_remove, seen = [], set()
  reports = []
  for atom in hier.atoms():
    ag = atom.parent()
    key = (ag.resname.strip().upper(), atom.name.strip().upper())
    if key not in _COORD_DEPROT_RULES:
      continue
    donor_xyz = matrix.col(atom.xyz)
    d_min = min((donor_xyz - mxyz).length() for mxyz in metal_xyzs)
    if d_min > coord_cutoff_A:
      continue
    rg = ag.parent()
    ch = rg.parent()
    h_names = _COORD_DEPROT_RULES[key]
    for h_atom in ag.atoms():
      # A bidentate carboxylate matches under both O keys (OD1 and OD2 map to
      # the same H), so dedup by atom identity or the H is removed twice.
      if (h_atom.name.strip().upper() in h_names
          and h_atom.memory_id() not in seen):
        seen.add(h_atom.memory_id())
        to_remove.append((ag, h_atom))
        reports.append((ch.id.strip(), rg.resseq.strip(),
                        key[0], key[1], h_atom.name.strip(), d_min))
  for ag, h_atom in to_remove:
    ag.remove_atom(h_atom)

  # Histidine: tautomer-aware. The coordinating ring N is the one nearest a
  # metal; the partner ring N is ~4 A across the ring, so the nearest-N pick
  # does the right thing even for a rare both-within-cutoff bridging His.
  his_reports = []
  for atom_group in hier.atom_groups():
    if atom_group.resname.strip().upper() != "HIS":
      continue
    abn = {a.name.strip().upper(): a for a in atom_group.atoms()}
    coord_n, d_coord = None, None
    for nn in ("ND1", "NE2"):
      na = abn.get(nn)
      if na is None:
        continue
      nxyz = matrix.col(na.xyz)
      dmin = min((nxyz - mxyz).length() for mxyz in metal_xyzs)
      if d_coord is None or dmin < d_coord:
        coord_n, d_coord = nn, dmin
    if coord_n is None or d_coord > coord_cutoff_A:
      continue
    coord_h = next((abn[h] for h in _HIS_RING_H[coord_n] if h in abn), None)
    if coord_h is None:
      continue  # coordinating N already H-free: correct neutral tautomer
    rg = atom_group.parent()
    ch = rg.parent()
    action = _relocate_or_remove_his_h(atom_group, coord_n, coord_h, d_coord)
    his_reports.append((ch.id.strip(), rg.resseq.strip(), action))

  # Canonicalise atom order (a relocated His proton is renamed in place).
  hier.sort_atoms_in_place()

  if reports:
    print(f"  deprotonated {len(reports)} metal-coordinating donor(s):",
          file=log)
    for ch, rs, rn, donor, hname, d in reports:
      print(f"    {rn} {ch}{rs} {donor} (d={d:.2f} A): removed {hname}",
            file=log)
  if his_reports:
    print(f"  corrected {len(his_reports)} metal-coordinating His:", file=log)
    for ch, rs, action in his_reports:
      print(f"    HIS {ch}{rs}: {action}", file=log)
  return len(reports), len(his_reports)
