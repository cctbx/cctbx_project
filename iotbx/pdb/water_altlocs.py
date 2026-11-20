"""Altloc labels for unlabelled partial-occupancy alternates around waters.

Deposited models often hold two waters, or a water and a ligand, at partial
occupancy on overlapping sites without altloc labels: they are alternates
(their occupancies share a site), but every tool that reads altlocs sees one
conformer with a clash. run() gives such atoms altloc labels and changes
nothing else: coordinates, occupancies and B factors stay, so the scattering
of the model is unchanged.

Three kinds of overlap (pairs within distance_cutoff, symmetry mates
included), where "unlabelled partial" means a non-H atom with occupancy
below 1 and a blank altloc:
1. water - water, both unlabelled partial: both are labelled;
2. water (unlabelled partial) - non-polymer residue whose atoms all have one
   partial occupancy (a ligand at 0.6 next to a water at 0.4): both are
   labelled, the residue as a whole;
3. water (unlabelled partial) - an atom of a residue in conformers: the water
   takes the one conformer label the residue does not have at that atom
   (a water next to OG B of a Ser in A and B becomes A); a residue with
   three or more labels is reported for manual work.
Metal ions are never alternates of a water: a water 2.1 A from a partial Mg
coordinates it. Overlaps of a water with a blank partial atom of anything else (a protein
atom at 0.5 without labels) are reported, not labelled. The overlap graph of
1 and 2 is coloured per connected component with the fewest labels (A, B,
C, ...), ligands first, then by occupancy, keeping labels already present.
"""
from __future__ import absolute_import, division, print_function
from cctbx.eltbx import chemical_elements
from libtbx import group_args
from iotbx.pdb import common_residue_names_get_class

LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
# Two waters, or a water and an atom of an alternative conformation, closer
# than this cannot be present at the same time (shortest real contacts of
# water with N, O are about 2.4 A). Water building in
# mmtbx.solvent.ordered_solvent uses the same number so that its output
# needs no labelling here.
DISTANCE_CUTOFF = 2.2
METALS = chemical_elements.metals_upper_set()

def is_metal(atom):
  """
  True if atom is a metal: water can be as close as about 2 A to it.
  """
  return atom.element.strip().upper() in METALS
POLYMER = ("common_amino_acid", "modified_amino_acid", "d_amino_acid",
           "common_rna_dna", "modified_rna_dna", "ccp4_mon_lib_rna_dna")

def _residue_class(atom):
  return common_residue_names_get_class(name=atom.parent().resname)

def _is_water(atom):
  return _residue_class(atom) == "common_water"

def _unlabelled_partial(atom):
  return (atom.occ < 1.0 and not atom.parent().altloc.strip()
          and not atom.element_is_hydrogen())

def _uniform_nonpolymer(atom):
  if _residue_class(atom) in POLYMER + ("common_water",): return False
  return len(set(round(a.occ, 2) for a in atom.parent().atoms())) == 1

def _label(atom):
  ag = atom.parent()
  rg = ag.parent()
  return "%s %s %s%s" % (ag.resname.strip(), rg.parent().id.strip(),
                         rg.resseq.strip(), rg.icode.strip())

def run(hierarchy, crystal_symmetry=None, distance_cutoff=DISTANCE_CUTOFF):
  """
  Label the alternates of hierarchy in place (module docstring).
  crystal_symmetry: for the symmetry mates; None: a box around the model.
  Returns group_args(n_labelled, labelled: [(residue label, altloc)] in
  model order, manual: [text] of the overlaps left to the user).
  """
  from cctbx import crystal
  atoms = hierarchy.atoms()
  if crystal_symmetry is None or crystal_symmetry.unit_cell() is None:
    crystal_symmetry = crystal.non_crystallographic_symmetry(
      sites_cart=atoms.extract_xyz(), buffer_layer=5)
  xs = hierarchy.extract_xray_structure(crystal_symmetry=crystal_symmetry)
  pst = xs.pair_asu_table(
    distance_cutoff=distance_cutoff).extract_pair_sym_table()
  def node(atom):
    if _is_water(atom): return ("w", atom.i_seq)
    return ("r", atom.parent().memory_id())
  edges, owner, forced, manual = {}, {}, {}, []
  def report(line):
    if line not in manual: manual.append(line)
  def connect(a, b):
    na, nb = node(a), node(b)
    edges.setdefault(na, set()).add(nb)
    edges.setdefault(nb, set()).add(na)
    owner[na], owner[nb] = a, b
  for i, pairs in enumerate(pst):
    ai = atoms[i]
    for j in pairs.keys():
      if j == i: continue # own symmetry copy (water on special position)
      aj = atoms[j]
      if ai.element_is_hydrogen() or aj.element_is_hydrogen(): continue
      w, x = sorted([ai, aj], key=lambda a: not _is_water(a))
      if not _is_water(w) or not _unlabelled_partial(w): continue
      if is_metal(x): continue
      if _is_water(x):
        if _unlabelled_partial(x): connect(w, x)
      elif _unlabelled_partial(x) and _uniform_nonpolymer(x):
        connect(w, x)
      elif x.parent().altloc.strip():
        own = sorted(set(ag.altloc for ag in x.parent().parent().atom_groups()
                         if ag.altloc.strip()))
        others = [l for l in own if l != x.parent().altloc]
        if len(others) == 1:
          forced.setdefault(w.i_seq, set()).add(others[0])
        else:
          report("%s overlaps %s %s (residue altlocs %s): choose its altloc"
            % (_label(w), x.name.strip(), _label(x), "".join(own)))
      else:
        report("%s (occ %.2f) overlaps blank partial %s %s (occ %.2f): label "
          "them" % (_label(w), w.occ, x.name.strip(), _label(x), x.occ))
  labelled = {}
  for i, letters in sorted(forced.items()):
    if len(letters) == 1:
      atoms[i].parent().altloc = list(letters)[0]
      labelled[i] = list(letters)[0]
    else:
      report("%s would need altlocs %s: choose one" % (
        _label(atoms[i]), "".join(sorted(letters))))
  seen = set()
  for start in sorted(edges):
    if start in seen: continue
    component, stack = [], [start]
    while stack:
      n = stack.pop()
      if n in seen: continue
      seen.add(n)
      component.append(n)
      stack.extend(edges[n])
    colour = {}
    for n in sorted(component, key=lambda n: (n[0] != "r", -owner[n].occ,
                                              owner[n].i_seq)):
      if owner[n].parent().altloc.strip():
        colour[n] = LETTERS.index(owner[n].parent().altloc)
        continue
      used = set(colour.get(m) for m in edges[n] if m in colour)
      c = 0
      while c in used: c += 1
      colour[n] = c
    for n, c in colour.items():
      a = owner[n]
      if a.parent().altloc.strip(): continue
      a.parent().altloc = LETTERS[c]
      labelled[a.i_seq] = LETTERS[c]
  out = [(_label(atoms[i]), labelled[i]) for i in sorted(labelled)]
  return group_args(n_labelled=len(out), labelled=out, manual=manual)
