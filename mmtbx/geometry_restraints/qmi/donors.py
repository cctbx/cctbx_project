"""Biological donor classification for the committing electronic predictor
(plan Sec 6).

Classifies each metal-coordinating atom by residue/element identity into the
refinement-ROBUST properties the predictor commits from: donor element, formal
charge, and structural role (porphyrin pyrrole N, strong axial candidate). It
reads residue names, elements, and the parent-residue grouping ONLY; it NEVER
uses metal-ligand bond lengths (refinement is what corrects those), QM energy,
or diffraction data. Residue->donor assignments follow the biological metal-site
donor surveys [9, 38]; the ferric tyrosinate diagnostic is [50].

Two calls:

- :func:`classify`: one coordinating atom -> its donor descriptor.
- :func:`summarize`: a list of descriptors (each carrying a ``res_id``) ->
  the coordination summary the predictor reads: heme (FeN4 topology) and its
  axial partition, plus the donor-hardness motif counts for the oxidation call.
"""

from libtbx import group_args
from iotbx.pdb import common_residue_names_get_class

from mmtbx.geometry_restraints.qmi import metals, util


# Side-chain donor residues in biological metal sites [9, 38].
_CARBOXYLATE = frozenset(("ASP", "GLU"))     # OD1/OD2, OE1/OE2: anionic hard O
_TYROSINATE = frozenset(("TYR",))  # OH: anionic hard O, ferric diag [50]
_THIOLATE = frozenset(("CYS",))              # SG: anionic soft S
_THIOETHER = frozenset(("MET",))             # SD: neutral soft S
_IMIDAZOLE = frozenset(("HIS",))             # ND1/NE2: neutral N
_HYDROXYL = frozenset(("SER", "THR"))        # OG/OG1: neutral O
_SIDE_AMIDE = frozenset(("ASN", "GLN"))      # OD1/OE1 carbonyl O: neutral
_AA_CLASSES = util.AA_RESIDUE_CLASSES


def classify(element, resname):
  """Classify one coordinating atom into a donor descriptor from residue/element
  identity alone (no bond lengths).

  Returns ``group_args(element, kind, charge, strong_axial)`` where ``charge``
  is the formal charge the donor's RESIDUE contributes (a carboxylate,
  tyrosinate or thiolate is -1 per residue; :func:`summarize` dedups by residue
  so a bidentate carboxylate counts once), and ``strong_axial`` flags a donor
  able to complete a low-spin heme axial pair (imidazole N, thioether S) [45].
  """
  el = metals.normalize_element(element)   # canonical Title-case symbol
  rn = resname.strip().upper()
  cls = common_residue_names_get_class(rn)

  if cls == "common_water" and el == "O":
    return group_args(element="O", kind="water_o", charge=0, strong_axial=False)
  if rn in _CARBOXYLATE and el == "O":
    return group_args(element="O", kind="carboxylate_o", charge=-1,
                      strong_axial=False)
  if rn in _TYROSINATE and el == "O":
    return group_args(element="O", kind="tyrosinate_o", charge=-1,
                      strong_axial=False)
  if rn in _THIOLATE and el == "S":
    return group_args(element="S", kind="cys_s", charge=-1, strong_axial=False)
  if rn in _THIOETHER and el == "S":
    return group_args(element="S", kind="met_s", charge=0, strong_axial=True)
  if rn in _IMIDAZOLE and el == "N":
    return group_args(element="N", kind="his_n", charge=0, strong_axial=True)
  if rn in _HYDROXYL and el == "O":
    return group_args(element="O", kind="hydroxyl_o", charge=0,
                      strong_axial=False)
  if rn in _SIDE_AMIDE and el == "O":
    return group_args(element="O", kind="amide_o", charge=0, strong_axial=False)
  # Main-chain donors of a standard residue: carbonyl O / amide N, both neutral.
  if cls in _AA_CLASSES and el == "O":
    return group_args(element="O", kind="backbone_o", charge=0,
                      strong_axial=False)
  if cls in _AA_CLASSES and el == "N":
    return group_args(element="N", kind="backbone_n", charge=0,
                      strong_axial=True)
  # Heteroatom / non-standard residue (cofactor, ligand, ion): keep the element,
  # unknown charge. An N here may be a porphyrin pyrrole N, resolved by the FeN4
  # topology test in summarize(), not by a residue-name whitelist.
  return group_args(element=el, kind=f"het_{el.lower()}", charge=0,
                    strong_axial=(el in ("N", "S")))


def summarize(donors):
  """Aggregate classified donors into the coordination summary the predictor
  commits from. Each donor must carry a hashable ``res_id`` (its parent residue)
  in addition to the :func:`classify` fields.

  Heme is detected by FeN4 topology: a SINGLE residue contributing four or more
  N donors (a porphyrin/chlorin/corrin macrocycle), which distinguishes a heme
  (4 N from one residue) from a 4-histidine site (4 N from four residues)
  without a residue-name whitelist. The non-porphyrin donors are the axial set.
  Motif counts (per residue, so a bidentate carboxylate counts once) drive the
  Fe oxidation call.

  UNCLASSIFIED DONORS. :func:`classify` recognises donors by residue name, so a
  donor on a cofactor, buffer or substrate ligand falls through to ``het_*``
  with charge 0 and contributes to no motif count. That is silent: the predictor
  then reads zero carboxylates and ``q_hard = 0`` for a site that may be
  strongly anionic. 6G7Q (an Fe chelated by two citrates) produced a summary
  bit-for-bit identical to passing no donors at all, and committed the ferrous
  fallthrough. So the count is reported: ``n_unclassified`` excludes the
  porphyrin N that the FeN4 test above legitimately consumes, and
  ``motif_evidence`` is False when NO donor contributed any motif count, i.e.
  the summary carries no information. See
  docs/donor_classification_generalisation.md.

  BOND-GRAPH CHARGE. A donor may additionally carry ``group_charge`` (the formal
  charge on its resonance group) and ``group_id`` (that group's identity), read
  from the bond graph by :mod:`mmtbx.geometry_restraints.qmi.donor_perception`.
  It is used ONLY for donors the residue-name pass could not place, and only for
  hard O, so the motif counts, both HIGH-confidence rules, and every state the
  whitelist already gets right are bit-for-bit unaffected; the charge fills the
  hole the whitelist leaves rather than replacing what it knows.

  Each anionic group contributes exactly -1 however negative it is, so a
  carbonate at -2 does not count double. Note the unit is the resonance GROUP,
  not the residue: for the twenty standard residues the two coincide (one
  carboxylate per Asp/Glu), so nothing the whitelist counts is rescaled, but a
  multi-arm het ligand contributes one -1 PER COORDINATING ARM. That is
  deliberate, and it is what makes a mono-citrate reach the ferric threshold on
  its own; the consequence to be aware of is that a single het residue can
  therefore reach ``q_hard = -2``, which no single whitelisted residue can.
  Restricting to O keeps the rule's HSAB reading [37, 44] honest and, as a side
  effect, excludes the porphyrin macrocycle (whose pyrrole N carry the ring
  charge) from a heme's ``q_hard``.

  Returns group_args(heme, n_axial, axial_all_strong, n_carboxylate,
  n_tyrosinate, n_cys_s, n_his, n_water_o, q_hard, n_unclassified,
  unclassified, motif_evidence, n_anionic_unplaced).
  """
  donors = list(donors)
  # Indices, not object identity: a caller may legitimately pass the same
  # descriptor object more than once (e.g. [desc] * 4), and id() would then
  # collapse them into one.
  n_idx_by_res = {}
  for i, d in enumerate(donors):
    if d.element == "N":
      n_idx_by_res.setdefault(d.res_id, []).append(i)
  heme_res = next((r for r, ix in n_idx_by_res.items() if len(ix) >= 4), None)
  heme = heme_res is not None
  axial = [d for d in donors if d.res_id != heme_res] if heme else []

  # Motif counts, deduped per residue (one carboxylate/tyrosinate/His per
  # residue).
  kinds_by_res = {}
  for d in donors:
    kinds_by_res.setdefault(d.res_id, set()).add(d.kind)
  def _count(kind):
    return sum(1 for kinds in kinds_by_res.values() if kind in kinds)
  n_carboxylate = _count("carboxylate_o")
  n_tyrosinate = _count("tyrosinate_o")
  n_cys_s = _count("cys_s")
  n_his = _count("his_n")
  n_water_o = _count("water_o")

  # Donors classify() could not place. A porphyrin N is het_n but IS consumed,
  # by the FeN4 test above, so exclude the macrocycle from the count.
  consumed = set()
  if heme:
    # Only the four N the FeN4 test consumed are explained; a fifth N on the
    # same residue (an intramolecular axial, as in cobalamin) or a macrocycle
    # propionate O is NOT, and must still be reported.
    consumed = set(n_idx_by_res[heme_res][:4])
  unplaced = [d for i, d in enumerate(donors)
              if str(d.kind).startswith("het_") and i not in consumed]
  # `heme` is deliberately NOT evidence here: it is derived from het_n donors
  # the residue-name classifier failed to place, so counting it would report a
  # 100%-unclassified site (a chlorophyll Mg, an Fe tetraaza chelate) as
  # informative. Only donors classify() actually placed count.
  motif_evidence = any((n_carboxylate, n_tyrosinate, n_cys_s, n_his, n_water_o))

  # Anionic hard-O groups among the donors the residue-name pass could not place
  # (see BOND-GRAPH CHARGE above). Deduped by resonance group, so the two O of
  # one carboxylate arm count once and the two arms of a citrate count twice.
  anionic_groups = set()
  for d in unplaced:
    q = getattr(d, "group_charge", None)
    if d.element == "O" and q is not None and q < 0:
      gid = getattr(d, "group_id", None)
      # Fall back to object identity only when there is genuinely no group id;
      # an EMPTY group id is still an id, and `or` would silently drop two
      # donors of one group into two separate counts.
      anionic_groups.add(gid if gid is not None else id(d))
  n_anionic_unplaced = len(anionic_groups)

  return group_args(
    heme=heme,
    n_axial=len(axial),
    axial_all_strong=(len(axial) >= 2 and all(d.strong_axial for d in axial)),
    n_carboxylate=n_carboxylate,
    n_tyrosinate=n_tyrosinate,
    n_cys_s=n_cys_s,
    n_his=n_his,
    # Donors the residue-name classifier could not place, and whether ANY motif
    # evidence survived. motif_evidence False means the predictor is deciding
    # with no donor information at all.
    n_unclassified=len(unplaced),
    unclassified=sorted({
      f"{d.res_id[-1] if isinstance(d.res_id, tuple) else d.res_id}"
      f".{d.element}" for d in unplaced}),
    motif_evidence=motif_evidence,
    # Anionic hard-O groups recovered from the bond graph for donors the
    # residue-name pass could not place (0 when no bond graph was available).
    n_anionic_unplaced=n_anionic_unplaced,
    # Coordinating water O (aquo/hydroxo; both modelled as HOH). For a heme
    # every non-porphyrin donor is axial, so this is the axial-water count
    # there.
    n_water_o=n_water_o,
    # Hard anionic O charge, one -1 per group: the carboxylates and tyrosinates
    # the residue names identify, plus the anionic O groups only the bond graph
    # can see (a citrate, a phosphate, a carbonate).
    q_hard=-(n_carboxylate + n_tyrosinate + n_anionic_unplaced))
