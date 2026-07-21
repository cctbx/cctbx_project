"""Biological donor classification for the committing electronic predictor
(plan Sec 6).

Classifies each metal-coordinating atom by residue/element identity into the
refinement-ROBUST properties the predictor commits from: donor element,
qualitative ligand-field strength, formal charge, and structural role
(porphyrin pyrrole N, strong axial candidate). It reads residue names, atom
names, elements, and the parent-residue grouping ONLY; it NEVER uses metal-
ligand bond lengths (refinement is what corrects those), QM energy, or
diffraction data. Residue->donor assignments follow the biological metal-site
donor surveys [9, 38]; the ferric tyrosinate diagnostic is [50].

Two calls:

- :func:`classify` -- one coordinating atom -> its donor descriptor.
- :func:`summarize` -- a list of descriptors (each carrying a ``res_id``) ->
  the coordination summary the predictor reads: heme (FeN4 topology) and its
  axial partition, plus the donor-hardness motif counts for the oxidation call.
"""

from __future__ import absolute_import, division, print_function

from libtbx import group_args
from iotbx.pdb import common_residue_names_get_class


WEAK, INTERMEDIATE, STRONG = "weak", "intermediate", "strong"

# Side-chain donor residues in biological metal sites [9, 38].
_CARBOXYLATE = frozenset(("ASP", "GLU"))     # OD1/OD2, OE1/OE2: anionic hard O
_TYROSINATE = frozenset(("TYR",))            # OH: anionic hard O, ferric diag [50]
_THIOLATE = frozenset(("CYS",))              # SG: anionic soft S
_THIOETHER = frozenset(("MET",))             # SD: neutral soft S
_IMIDAZOLE = frozenset(("HIS",))             # ND1/NE2: neutral N
_HYDROXYL = frozenset(("SER", "THR"))        # OG/OG1: neutral O
_SIDE_AMIDE = frozenset(("ASN", "GLN"))      # OD1/OE1 carbonyl O: neutral
_AA_CLASSES = frozenset(("common_amino_acid", "modified_amino_acid"))


def classify(element, resname, atom_name):
  """Classify one coordinating atom into a donor descriptor from residue/element
  identity alone (no bond lengths).

  Returns ``group_args(element, kind, field, charge, strong_axial)`` where
  ``charge`` is the formal charge the donor's RESIDUE contributes (a carboxylate,
  tyrosinate or thiolate is -1 per residue; :func:`summarize` dedups by residue
  so a bidentate carboxylate counts once), and ``strong_axial`` flags a donor
  able to complete a low-spin heme axial pair (imidazole N, thioether S) [45].
  """
  el = element.strip().upper()
  rn = resname.strip().upper()
  cls = common_residue_names_get_class(rn)

  if cls == "common_water" and el == "O":
    return group_args(element="O", kind="water_o", field=WEAK, charge=0,
                      strong_axial=False)
  if rn in _CARBOXYLATE and el == "O":
    return group_args(element="O", kind="carboxylate_o", field=WEAK, charge=-1,
                      strong_axial=False)
  if rn in _TYROSINATE and el == "O":
    return group_args(element="O", kind="tyrosinate_o", field=WEAK, charge=-1,
                      strong_axial=False)
  if rn in _THIOLATE and el == "S":
    return group_args(element="S", kind="cys_s", field=WEAK, charge=-1,
                      strong_axial=False)
  if rn in _THIOETHER and el == "S":
    return group_args(element="S", kind="met_s", field=INTERMEDIATE, charge=0,
                      strong_axial=True)
  if rn in _IMIDAZOLE and el == "N":
    return group_args(element="N", kind="his_n", field=INTERMEDIATE, charge=0,
                      strong_axial=True)
  if rn in _HYDROXYL and el == "O":
    return group_args(element="O", kind="hydroxyl_o", field=WEAK, charge=0,
                      strong_axial=False)
  if rn in _SIDE_AMIDE and el == "O":
    return group_args(element="O", kind="amide_o", field=WEAK, charge=0,
                      strong_axial=False)
  # Main-chain donors of a standard residue: carbonyl O / amide N, both neutral.
  if cls in _AA_CLASSES and el == "O":
    return group_args(element="O", kind="backbone_o", field=WEAK, charge=0,
                      strong_axial=False)
  if cls in _AA_CLASSES and el == "N":
    return group_args(element="N", kind="backbone_n", field=INTERMEDIATE,
                      charge=0, strong_axial=True)
  # Heteroatom / non-standard residue (cofactor, ligand, ion): keep the element,
  # unknown charge. An N here may be a porphyrin pyrrole N -- resolved by the
  # FeN4 topology test in summarize(), not by a residue-name whitelist.
  field = {"N": INTERMEDIATE}.get(el, WEAK)
  return group_args(element=el, kind="het_%s" % el.lower(), field=field,
                    charge=0, strong_axial=(el in ("N", "S")))


def summarize(donors):
  """Aggregate classified donors into the coordination summary the predictor
  commits from. Each donor must carry a hashable ``res_id`` (its parent residue)
  in addition to the :func:`classify` fields.

  Heme is detected by FeN4 topology -- a SINGLE residue contributing four or more
  N donors (a porphyrin/chlorin/corrin macrocycle), which distinguishes a heme
  (4 N from one residue) from a 4-histidine site (4 N from four residues) without
  a residue-name whitelist. The non-porphyrin donors are the axial set. Motif
  counts (per residue, so a bidentate carboxylate counts once) drive the Fe
  oxidation call.

  Returns group_args(heme, n_axial, axial_all_strong, n_carboxylate,
  n_tyrosinate, n_cys_s, n_his, q_hard).
  """
  donors = list(donors)
  n_by_res = {}
  for d in donors:
    if d.element == "N":
      n_by_res.setdefault(d.res_id, []).append(d)
  heme_res = next((r for r, ns in n_by_res.items() if len(ns) >= 4), None)
  heme = heme_res is not None
  axial = [d for d in donors if d.res_id != heme_res] if heme else []

  # Motif counts, deduped per residue (one carboxylate/tyrosinate/His per residue).
  kinds_by_res = {}
  for d in donors:
    kinds_by_res.setdefault(d.res_id, set()).add(d.kind)
  def _count(kind):
    return sum(1 for kinds in kinds_by_res.values() if kind in kinds)
  n_carboxylate = _count("carboxylate_o")
  n_tyrosinate = _count("tyrosinate_o")
  n_cys_s = _count("cys_s")
  n_his = _count("his_n")

  return group_args(
    heme=heme,
    n_axial=len(axial),
    axial_all_strong=(len(axial) >= 2 and all(d.strong_axial for d in axial)),
    n_carboxylate=n_carboxylate,
    n_tyrosinate=n_tyrosinate,
    n_cys_s=n_cys_s,
    n_his=n_his,
    # Hard anionic O charge (carboxylate + tyrosinate), one -1 per group.
    q_hard=-(n_carboxylate + n_tyrosinate))
