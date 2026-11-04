"""
Complete the conformers of residues with partial alternative conformations,
so that every conformer label of a residue has all of the residue's atoms.
Implementation of iotbx.pdb.hierarchy.root.complete_altlocs.

Conformer labels of a residue group are the altloc labels with at least one
non-hydrogen atom; hydrogen-only labels (e.g. neutron H/D exchange sites)
are not conformers and neither borrow nor lend atoms. Conformer X of a
residue group starts as its blank atoms plus its atoms labelled X. Every
atom name the residue group has under another label, but neither blank nor
X, is a gap. A gap is filled by a candidate, the atom of that name under
another label P of an atom group with the same residue name as conformer
X (a conformer never borrows from a different residue of a
microheterogeneity site), when the candidate is bonded to an atom already
in conformer X: closer than the sum of covalent radii plus BOND_TOLERANCE.
A non-hydrogen candidate must be bonded to a non-hydrogen atom of
conformer X. A hydrogen candidate is accepted only if its parent, the
nearest non-hydrogen atom bonded to it among the atoms of label P and the
blank atoms, is already in conformer X: a truncated conformer with its caps
lends nothing and borrows nothing through them.
An accepted atom joins conformer X and can anchor further candidates;
filling stops when no candidate is accepted. When more than one candidate
of the same name is bonded, the candidate must have occupancy to give:
its occupancy must be at least the occupancy of label X plus that of label
P plus that of every other label already borrowing it, less
OCCUPANCY_TOLERANCE. The occupancy of a label in a residue group is the
smallest occupancy among its non-hydrogen atoms. If not exactly one
candidate passes, the case is ambiguous and Sorry is raised. Gaps that stay
open are reported as truncations. Borrowing never crosses residue groups.

An atom borrowed by every other conformer label of the residue group is
relabelled blank (shared) with its occupancy unchanged. Otherwise the atom
is copied into each borrowing conformer at the same coordinates and
B-factor, and its occupancy is divided among its own label and the
borrowers in proportion to their label occupancies; a residue group whose
label occupancies do not add up to the source atom's occupancy (within
OCCUPANCY_TOLERANCE) is reported as inconsistent, and one whose label
occupancies are all zero raises Sorry. Every atom site keeps its total
occupancy, so structure factors are unchanged.
"""
import math
from libtbx import group_args
from libtbx.utils import null_out, Sorry

BOND_TOLERANCE = 0.4
OCCUPANCY_TOLERANCE = 0.05

def _distance(a, b):
  return math.sqrt(sum([(x - y)**2 for x, y in zip(a.xyz, b.xyz)]))

def bonded(atom_a, atom_b):
  """
  Distance-based covalent bond test.

  Parameters
  ----------
  atom_a, atom_b : iotbx.pdb.hierarchy.atom

  Returns
  -------
  bool
    True if the atoms are closer than the sum of their covalent radii plus
    BOND_TOLERANCE.
  """
  from iotbx.pdb.altloc_reconciliation import _element, _covalent_radius
  limit = _covalent_radius(_element(atom_a)) + \
    _covalent_radius(_element(atom_b)) + BOND_TOLERANCE
  return _distance(atom_a, atom_b) < limit

def conformer_labels(residue_group):
  """
  Altloc labels of a residue group that have at least one non-hydrogen atom.

  Parameters
  ----------
  residue_group : iotbx.pdb.hierarchy.residue_group

  Returns
  -------
  list of str
    Sorted labels.
  """
  result = set()
  for ag in residue_group.atom_groups():
    if (ag.altloc.strip() == ""): continue
    for atom in ag.atoms():
      if (not atom.element_is_hydrogen()):
        result.add(ag.altloc)
        break
  return sorted(result)

def label_occupancy(residue_group, label):
  """
  Occupancy of a conformer label in a residue group.

  Parameters
  ----------
  residue_group : iotbx.pdb.hierarchy.residue_group
  label : str
    Altloc label.

  Returns
  -------
  float or None
    The smallest occupancy among the non-hydrogen atoms with that label,
    None if there are none.
  """
  result = None
  for ag in residue_group.atom_groups():
    if (ag.altloc != label): continue
    for atom in ag.atoms():
      if (atom.element_is_hydrogen()): continue
      if (result is None or atom.occ < result): result = atom.occ
  return result

def _residue_id(atom_group):
  """Residue id with resname and a blank altloc position, e.g. ' SER A 205 '"""
  return " " + atom_group.id_str()[1:]

def _nearest(atom, atoms):
  """(distance, atom) of the atom in atoms nearest to atom, None if empty"""
  result = None
  for other in atoms:
    d = _distance(atom, other)
    if (result is None or d < result[0]): result = (d, other)
  return result

def _parent_in(atom, rg, atoms):
  """True if the parent of hydrogen atom, the nearest non-hydrogen atom
  bonded to it in its own atom group or a blank one of rg, is in atoms"""
  heavy = [a for ag in rg.atom_groups()
    if ag.altloc == "" or ag.memory_id() == atom.parent().memory_id()
    for a in ag.atoms() if not a.element_is_hydrogen()]
  parent = _nearest(atom, [a for a in heavy if bonded(atom, a)])
  return parent is not None and parent[1] in atoms

class _borrow(object):
  """An atom of label source_label borrowed by one or more labels"""
  def __init__(self, atom, source_label):
    self.atom = atom
    self.source_label = source_label
    self.residue_id = _residue_id(atom.parent())
    self.borrowers = [] # (label, atom_group)

def _complete_conformer(rg, labels, occ_r, atom_group, borrows, truncated):
  """Fill the gaps of the conformer of atom_group; record accepted atoms in
  borrows (dict keyed by (label, resname, name)) and open gaps in truncated"""
  x = atom_group.altloc
  rx = atom_group.resname
  in_x = [a for ag in rg.atom_groups() if ag.altloc == ""
    for a in ag.atoms()]
  in_x.extend(atom_group.atoms())
  names_in_x = set([a.name for a in in_x])
  sources = [ag for ag in rg.atom_groups()
    if ag.altloc in labels and ag.altloc != x and ag.resname == rx]
  gaps = set()
  for ag in sources:
    for atom in ag.atoms():
      if (atom.name not in names_in_x): gaps.add(atom.name)
  accepted = True
  while (accepted):
    accepted = False
    for name in sorted(gaps):
      passing = [] # (label, atom, distance, nearest bonded anchor)
      for ag in sources:
        for atom in ag.atoms():
          if (atom.name != name): continue
          if (atom.element_is_hydrogen() and
              not _parent_in(atom, rg, in_x)): continue
          # a hydrogen never anchors a heavy atom
          bonded_anchor = _nearest(atom, [b for b in in_x
            if bonded(atom, b) and (atom.element_is_hydrogen() or
              not b.element_is_hydrogen())])
          if (bonded_anchor is None): continue
          d, anchor = bonded_anchor
          passing.append((ag.altloc, atom, d, anchor))
      if (len(passing) == 0): continue
      if (len(passing) > 1):
        budget_passing = []
        for p, atom, d, anchor in passing:
          key = (p, rx, name)
          need = occ_r[x] + occ_r[p] - OCCUPANCY_TOLERANCE
          if (key in borrows):
            for y, ag_y in borrows[key].borrowers:
              need += occ_r[y]
          if (atom.occ >= need):
            budget_passing.append((p, atom, d, anchor))
        if (len(budget_passing) != 1):
          raise Sorry("Cannot complete conformer %s of residue \"%s\": "
            "atom %s could come from %s; fix the model by hand" % (
            x, _residue_id(atom_group), name.strip(),
            ", ".join(["altloc %s (%.2f A from %s, occ %.2f)" % (
              p, d, anchor.name.strip(), atom.occ)
              for p, atom, d, anchor in passing])))
        passing = budget_passing
      p, atom, d, anchor = passing[0]
      key = (p, rx, name)
      if (key not in borrows): borrows[key] = _borrow(atom, p)
      borrows[key].borrowers.append((x, atom_group))
      in_x.append(atom)
      gaps.remove(name)
      accepted = True
  if (len(gaps) > 0):
    names = sorted(gaps)
    nearest = None # (distance, label, anchor)
    for ag in sources:
      for atom in ag.atoms():
        if (atom.name != names[0]): continue
        d, anchor = _nearest(atom, in_x)
        if (nearest is None or d < nearest[0]):
          nearest = (d, ag.altloc, anchor)
    d, p, anchor = nearest
    detail = (names[0].strip(), p, d, anchor.name.strip(), x)
    reason = "%s: %s copy is %.1f A from %s %s" % detail
    truncated.append(
      (_residue_id(atom_group), x, [n.strip() for n in names], reason, detail))

def _blank_atom_group(rg, resname):
  from iotbx.pdb import hierarchy
  for ag in rg.atom_groups():
    if (ag.altloc == "" and ag.resname == resname): return ag
  ag = hierarchy.atom_group(altloc="", resname=resname)
  rg.insert_atom_group(0, ag)
  return ag

def run(hierarchy, check_only=False, log=None):
  """
  Complete partial alternative conformations of every residue group.

  Parameters
  ----------
  hierarchy : iotbx.pdb.hierarchy.root
    Edited in place unless check_only.
  check_only : bool
    Report what would change without editing the hierarchy.
  log : file-like or None
    Receives report_lines(result), one per line.

  Returns
  -------
  libtbx.group_args
    Residue ids are atom group id_str with a blank altloc position, e.g.
    " SER A 205 ".
    shared: list of (residue id, atom name, source label);
    shared_needed_by: per shared entry, the sorted labels that borrowed it;
    copied: list of (residue id, atom name, source label, target label,
    occupancy of the copy); truncated: list of (residue id, label,
    atom names, reason); truncated_details: per truncated entry, (atom
    name, source label, distance, anchor atom name, label);
    inconsistent: list of residue id; changed: bool.

  Raises
  ------
  Sorry
    If a gap can be filled from more than one label and occupancies do not
    tell which one, or if the occupancy of an atom must be split between
    labels whose occupancies are all zero.
  """
  if (log is None): log = null_out()
  plans = [] # (rg, labels, occ_r, borrows)
  truncated = []
  for rg in hierarchy.residue_groups():
    labels = conformer_labels(rg)
    if (len(labels) < 2): continue
    occ_r = dict([(label, label_occupancy(rg, label)) for label in labels])
    borrows = {}
    for ag in rg.atom_groups():
      if (ag.altloc not in labels): continue
      _complete_conformer(rg=rg, labels=labels, occ_r=occ_r, atom_group=ag,
        borrows=borrows, truncated=truncated)
    if (len(borrows) > 0): plans.append((rg, labels, occ_r, borrows))
  for rg, labels, occ_r, borrows in plans:
    for b in borrows.values():
      borrower_labels = [y for y, ag_y in b.borrowers]
      if (set(borrower_labels) == set(labels) - set([b.source_label])):
        continue # shared: occupancy kept, nothing to split
      if (sum([occ_r[y] for y in [b.source_label] + borrower_labels]) == 0):
        raise Sorry("all conformers of residue \"%s\" have zero occupancy; "
          "cannot split" % b.residue_id)
  shared, shared_needed_by, copied, inconsistent = [], [], [], []
  for rg, labels, occ_r, borrows in plans:
    for b in borrows.values():
      atom, p, rg_id = b.atom, b.source_label, b.residue_id
      name = atom.name.strip()
      borrower_labels = [y for y, ag_y in b.borrowers]
      if (set(borrower_labels) == set(labels) - set([p])):
        shared.append((rg_id, name, p))
        shared_needed_by.append(sorted(borrower_labels))
        if (check_only): continue
        ag_p = atom.parent()
        blank = _blank_atom_group(rg, ag_p.resname)
        moved = atom.detached_copy()
        ag_p.remove_atom(atom)
        blank.append_atom(moved)
        continue
      weights = [occ_r[p]] + [occ_r[y] for y in borrower_labels]
      total = atom.occ
      if (abs(sum(weights) - total) > OCCUPANCY_TOLERANCE and
          rg_id not in inconsistent):
        inconsistent.append(rg_id)
      for (y, ag_y), w in zip(b.borrowers, weights[1:]):
        occ = total * w / sum(weights)
        copied.append((rg_id, name, p, y, occ))
        if (check_only): continue
        copy = atom.detached_copy()
        copy.set_occ(occ)
        ag_y.append_atom(copy)
      if (not check_only): atom.set_occ(total * weights[0] / sum(weights))
    # only after every edit: a later copy may still go into such a group
    if (not check_only):
      for ag in rg.atom_groups():
        if (ag.atoms_size() == 0): rg.remove_atom_group(ag)
  changed = len(shared) + len(copied) > 0
  if (changed and not check_only):
    hierarchy.atoms_reset_serial()
    hierarchy.atoms().reset_i_seq()
  result = group_args(
    shared            = shared,
    shared_needed_by  = shared_needed_by,
    copied            = copied,
    truncated         = [t[:4] for t in truncated],
    truncated_details = [t[4] for t in truncated],
    inconsistent      = inconsistent,
    changed           = changed)
  for line in report_lines(result):
    print(line, file=log)
  return result

def report_lines(result):
  """
  Log lines for the result of run.

  Parameters
  ----------
  result : libtbx.group_args
    As returned by run.

  Returns
  -------
  list of str
    One line per shared atom, per residue group and pair of labels copied,
    and per truncated conformer; empty if nothing happened.
  """
  lines = []
  for (rg_id, name, p), needed_by in zip(result.shared,
                                         result.shared_needed_by):
    lines.append("%s %s: altloc %s -> shared (needed by %s)" % (
      rg_id, name, p, ", ".join(needed_by)))
  pairs = [] # (rg_id, source label, target label)
  names = {}
  for rg_id, name, p, y, occ in result.copied:
    key = (rg_id, p, y)
    if (key not in names):
      pairs.append(key)
      names[key] = []
    names[key].append(name)
  for rg_id, p, y in pairs:
    lines.append("%s %s: copied %s -> %s" % (
      rg_id, ", ".join(names[(rg_id, p, y)]), p, y))
  for rg_id, (name, p, d, anchor, x) in zip(
      [t[0] for t in result.truncated], result.truncated_details):
    lines.append("%s %s: %s copy is %.1f A from %s, not borrowed; "
      "truncated in %s" % (rg_id, name, p, d, "%s %s" % (anchor, x), x))
  return lines
