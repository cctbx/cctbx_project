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
  proton. A coordinating Cu/Ni N-terminus is deprotonated to its ATCUN form: a
  backbone amide becomes an amidate and the free alpha-amino terminus becomes
  -NH2.
- ``protonate_waters``: place H on bare waters (H-bond and cation aware) via
  the naiad placer, so a coordinating water is neutral H2O with its H kept off
  the metal.
- ``prepare_model``: the single entry point the Program calls (reduce2 +
  water H + metal-donor deprotonation), returning a fresh model.

The Program's validate() guard uses ``mmtbx.model.manager.has_hd()`` to decide
whether the model is already protonated.

INTERIM: water H come from the qmi-local ``naiad`` copy of
``place_water_hydrogens``; that placer is being reworked into reduce2. Once
reduce2 places good water H, request it inside ``add_hydrogens`` and delete the
marked ``protonate_waters`` call in ``prepare_model`` (and the ``qmi/naiad``
copy). No other qmi code, and no caller, changes.
"""

from __future__ import absolute_import, division, print_function

import os
import shutil
import tempfile

from libtbx import easy_run
from libtbx.utils import null_out, Sorry
from scitbx import matrix
from iotbx.pdb import common_residue_names_get_class

from mmtbx.geometry_restraints.qmi import metals as qmi_metals
from mmtbx.geometry_restraints.qmi import util


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

# Cys-thiolate and Tyr-phenolate formation is metal-INDUCED: only a soft or
# borderline, Lewis-acidic transition-metal / Zn cation supplies the covalency
# that stabilises the anion. Hard, low-charge ions (Na, K, Ca, Mg) do NOT, and
# rarely bind S/O-phenol at all, so a Cys/Tyr near one must be left neutral
# (HSAB: Pearson 1963 [37]; Dudev & Lim 2003 [38]). Carboxylate (Asp/Glu) is
# anionic by pH regardless of the metal, so it is NOT gated. These two donors
# fire only when the coordinating metal is in _ANION_INDUCING_METALS.
_SOFT_DONOR_KEYS = {("CYS", "SG"), ("TYR", "OH")}

# The d-block transition metals + group 12 (Zn/Cd/Hg): acidic/soft enough to
# induce thiolate/phenolate. Excludes hard alkali/alkaline-earth (Na, K, Mg,
# Ca, ...) and the p-block metals (conservative; a Pb/Sn thiolate would need a
# user override). Cited [37, 38].
_ANION_INDUCING_METALS = frozenset({
  "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
  "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
  "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
})

# Metal-induced N-terminal deprotonation at Cu(II)/Ni(II) (the ATCUN motif [71]).
# Binding collapses the pKa of the free alpha-amino terminus and of the first
# backbone amides, so they deprotonate and the N coordinate. reduce2 protonates
# them as for a free peptide, so at one of these ions strip:
#   - the single peptide amide H (_BACKBONE_AMIDE_H) on a coordinating backbone N;
#   - one of the three N-terminal amine H (_NTERM_AMINE_H) on a coordinating
#     alpha-amino N (-NH3+ -> -NH2), freeing its lone pair for the metal.
# Kept to the ions that do this biologically: a broader set would risk
# deprotonating an N that merely sits near a metal.
_ATCUN_METALS = frozenset({"Cu", "Ni"})
# Backbone amide H (single peptide H on N) vs the N-terminal amine H (three, on
# the alpha-amino N); the two donors are handled separately below.
_BACKBONE_AMIDE_H = {"H", "HN", "D", "DN"}
_NTERM_AMINE_H = {"H1", "H2", "H3", "D1", "D2", "D3"}

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
    with util.in_directory(tmp):
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
# water protonation (naiad)
# ---------------------------------------------------------------------------
def _count_water_h(hier):
  n = 0
  for ag in hier.atom_groups():
    if common_residue_names_get_class(ag.resname) != "common_water":
      continue
    for a in ag.atoms():
      if a.element_is_hydrogen():
        n += 1
  return n


def protonate_waters(hier, n_refine=5, refine_tol=1, n_basin=0,
                     lone_pair_directed=False, joint=False,
                     report_stats=False, log=None):
  """Place the two H on every bare water IN PLACE, H-bond-aware and cation-aware
  (the H are kept off coordinating metals), via the naiad placer.

  reduce2 leaves waters bare, so a metal-coordinating water would otherwise be
  read as an oxide/hydroxide by the charge accounting; protonating it to a
  neutral H2O with its H directed away from the metal fixes the cluster charge.
  Idempotent (waters already carrying H are untouched). Returns the number of
  water H atoms added.

  *n_refine* / *refine_tol* / *n_basin* / *lone_pair_directed* / *joint* are the
  placement-quality knobs forwarded to
  :func:`...naiad.water_protonation.place_water_hydrogens` (defaults match it).
  With *report_stats*, the per-sweep water-H clash table is streamed and the
  residual inter-water H-H contacts are listed afterwards (grouped by
  1.5/1.8/2.0 A), mirroring the naiad water_protonation Program.
  """
  if log is None:
    log = null_out()
  from mmtbx.geometry_restraints.qmi.naiad import water_protonation
  before = _count_water_h(hier)

  on_state = None
  if report_stats:
    header = {"printed": False}
    def on_state(label, stats):
      if not header["printed"]:
        print(f"  water H clash report ({stats[0]} placed; H-H between "
              f"different waters):", file=log)
        header["printed"] = True
      water_protonation._clash_row(label, stats, log)

  water_protonation.place_water_hydrogens(
    hier, n_refine=n_refine, refine_tol=refine_tol, n_basin=n_basin,
    lone_pair_directed=lone_pair_directed, joint=joint, on_state=on_state)
  added = _count_water_h(hier) - before
  print(f"  placed {added} water hydrogen(s)", file=log)
  if report_stats:
    contacts = water_protonation._worst_water_clashes(hier)
    if not contacts:
      print("  no water-H contacts < 2.0 A", file=log)
    else:
      print("  residual water-H contacts:", file=log)
      for label, lo, hi in (("< 1.5 A", 0.0, 1.5),
                            ("1.5-1.8 A", 1.5, 1.8),
                            ("1.8-2.0 A", 1.8, 2.0)):
        band = [c for c in contacts if lo <= c[0] < hi]
        print(f"    {label} ({len(band)}):", file=log)
        for d, id_a, id_b in band:
          print(f"      {d:.2f} A  {id_a}  <->  {id_b}", file=log)
  return added


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
    # No ring geometry to place the partner-N proton; keep the coordinating-N
    # proton so the ring stays neutral (removing the only remaining ring proton
    # would leave both ring N bare, an imidazolate).
    old = coord_h.name.strip()
    return (f"kept {old} on {coord_n} (no ring geom to relocate, "
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
                             coord_cutoff=None,
                             log=None):
  """Fix, IN PLACE, the protonation of donors in *hier* coordinating any metal
  atom (within ``coord_cutoff(metal_element, donor_element)``).

  - Asp/Glu carboxylate (:data:`_COORD_DEPROT_RULES`): strip the acid H for any
    coordinating metal (anionic by pH, metal-independent).
  - Cys thiolate / Tyr phenolate: strip the donor H ONLY when coordinated by an
    anion-inducing metal (:data:`_ANION_INDUCING_METALS`, transition metals +
    Zn); a Cys/Tyr near a hard Na/K/Ca/Mg is left neutral (Pearson 1963 [37];
    Dudev & Lim 2003 [38]).
  - ATCUN N-terminal deprotonation at a coordinating Cu(II)/Ni(II)
    (:data:`_ATCUN_METALS`, ATCUN-type [71]): a standard-AA backbone N is a
    deprotonated peptide amidate (strip its amide H); a free alpha-amino terminus
    (-NH3+, all three H) binds as -NH2 (drop the H nearest the metal).
  - His: tautomer-aware (neutral, not imidazolate) via
    :func:`_relocate_or_remove_his_h`, for any coordinating metal.

  *metals* is the element set the donors are scored against (default: all
  METALS); the in-pipeline caller passes the metals qmi may model, so donors
  near an entirely out-of-scope metal are left untouched. The Cys/Tyr gate is
  applied on top of that, against the anion-inducing subset.

  Canonicalises atom order at the end (a relocated His proton is renamed in
  place). Returns ``(n_donors_corrected, n_his_corrected)``. A no-op (0, 0)
  when the model has no in-scope metals.
  """
  if log is None:
    log = null_out()
  if coord_cutoff is None:
    coord_cutoff = util.make_coordination_cutoff()
  metal_atoms = qmi_metals.metal_atoms(hier.atoms(), metals=metals)
  if not metal_atoms:
    return 0, 0
  # Cys/Tyr score only against anion-inducing metals (see _SOFT_DONOR_KEYS).
  anion_atoms = [a for a in metal_atoms
                 if qmi_metals.normalize_element(a.element) in _ANION_INDUCING_METALS]

  # Asp/Glu (any metal) and Cys/Tyr (anion-inducing metals only): strip the H.
  to_remove, seen = [], set()
  reports = []
  for atom in hier.atoms():
    ag = atom.parent()
    key = (ag.resname.strip().upper(), atom.name.strip().upper())
    if key not in _COORD_DEPROT_RULES:
      continue
    scored_metals = anion_atoms if key in _SOFT_DONOR_KEYS else metal_atoms
    if not scored_metals:
      continue
    donor_xyz = matrix.col(atom.xyz)
    nearest, d_min = min(
      ((m, (donor_xyz - matrix.col(m.xyz)).length()) for m in scored_metals),
      key=lambda md: md[1])
    if d_min > coord_cutoff(nearest.element, atom.element):
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

  # ATCUN (Cu/Ni) N-terminal deprotonation. (a) backbone amidate: a standard-AA
  # backbone N coordinating an ATCUN metal is a deprotonated peptide N, so strip
  # the amide H reduce2 placed. A backbone N carrying only H1/H2/H3 is an
  # N-terminal amine (handled in (b)) and is skipped here (no _BACKBONE_AMIDE_H).
  atcun_metals = [a for a in metal_atoms if qmi_metals.normalize_element(
                  a.element) in _ATCUN_METALS]
  for atom in (hier.atoms() if atcun_metals else ()):
    if atom.name.strip().upper() != "N":
      continue
    ag = atom.parent()
    if common_residue_names_get_class(ag.resname) not in util.AA_RESIDUE_CLASSES:
      continue
    h_atom = next((a for a in ag.atoms()
                   if a.name.strip().upper() in _BACKBONE_AMIDE_H
                   and a.memory_id() not in seen), None)
    if h_atom is None:
      continue
    n_xyz = matrix.col(atom.xyz)
    nearest, d_min = min(
      ((m, (n_xyz - matrix.col(m.xyz)).length()) for m in atcun_metals),
      key=lambda md: md[1])
    if d_min > coord_cutoff(nearest.element, atom.element):
      continue
    rg = ag.parent()
    ch = rg.parent()
    seen.add(h_atom.memory_id())
    to_remove.append((ag, h_atom))
    reports.append((ch.id.strip(), rg.resseq.strip(),
                    ag.resname.strip(), "N", h_atom.name.strip(), d_min))

  # (b) N-terminal amine: a coordinating alpha-amino N still carrying all three
  # -NH3+ protons binds as -NH2; drop the H nearest the metal so the vacated site
  # becomes the lone pair pointing at it. Requiring all three keeps this off a
  # backbone amide (single H, done in (a)) or an already-neutral amine.
  for atom in (hier.atoms() if atcun_metals else ()):
    if atom.name.strip().upper() != "N":
      continue
    ag = atom.parent()
    if common_residue_names_get_class(ag.resname) not in util.AA_RESIDUE_CLASSES:
      continue
    amine_h = [a for a in ag.atoms()
               if a.name.strip().upper() in _NTERM_AMINE_H
               and a.memory_id() not in seen]
    if len(amine_h) < 3:
      continue
    n_xyz = matrix.col(atom.xyz)
    nearest, d_min = min(
      ((m, (n_xyz - matrix.col(m.xyz)).length()) for m in atcun_metals),
      key=lambda md: md[1])
    if d_min > coord_cutoff(nearest.element, atom.element):
      continue
    m_xyz = matrix.col(nearest.xyz)
    h_atom = min(amine_h, key=lambda h: (matrix.col(h.xyz) - m_xyz).length())
    rg = ag.parent()
    ch = rg.parent()
    seen.add(h_atom.memory_id())
    to_remove.append((ag, h_atom))
    reports.append((ch.id.strip(), rg.resseq.strip(),
                    ag.resname.strip(), "N", h_atom.name.strip(), d_min))

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
    coord_n, d_coord, coord_metal = None, None, None
    for nn in ("ND1", "NE2"):
      na = abn.get(nn)
      if na is None:
        continue
      nxyz = matrix.col(na.xyz)
      m, dmin = min(
        ((mm, (nxyz - matrix.col(mm.xyz)).length()) for mm in metal_atoms),
        key=lambda md: md[1])
      if d_coord is None or dmin < d_coord:
        coord_n, d_coord, coord_metal = nn, dmin, m
    if coord_n is None or d_coord > coord_cutoff(coord_metal.element,
                                                 abn[coord_n].element):
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


def prepare_model(model, metals=None, water_kwargs=None, coord_cutoff=None,
                  log=None):
  """Protonate a raw model for qmi and fix metal-donor protonation.

  Runs reduce2 (protein H), places water H, and deprotonates metal-coordinating
  donors, returning a fresh model ready for cluster extraction. This is the one
  entry point the Program calls, so the reduce2-vs-naiad split is hidden from
  the caller: switching water H to reduce2 later touches only this function.

  *metals* scopes the deprotonation (see :func:`deprotonate_metal_donors`);
  *water_kwargs* is a dict forwarded to :func:`protonate_waters` (None uses its
  defaults); *coord_cutoff* is the shared metal-donor cutoff callable (None uses
  the fixed default).
  """
  import mmtbx.model
  if log is None:
    log = null_out()
  model = add_hydrogens(model, log=log)
  hier = model.get_hierarchy()
  # --- INTERIM: qmi-local naiad water protonation. When reduce2 places good
  #     water H, request it in add_hydrogens (so it returns a model with waters
  #     already protonated) and delete the next call; nothing else changes. ---
  protonate_waters(hier, log=log, **(water_kwargs or {}))
  # -------------------------------------------------------------------------
  deprotonate_metal_donors(hier, metals=metals, coord_cutoff=coord_cutoff,
                           log=log)
  return mmtbx.model.manager(
    model_input=None, pdb_hierarchy=hier,
    crystal_symmetry=model.crystal_symmetry())
