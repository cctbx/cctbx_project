"""
Split a model with alternative conformations into complete conformers.

The polymer (protein, nucleic acid) is rewritten so that every atom carries
an altloc label and each of the N labels found in the input ("A", "B", ...)
makes a complete copy of it. All conformers have the same number of polymer
atoms, except where two residue types share one position
(microheterogeneity). Ligands and ions are not split: their atoms, copies
and occupancies stay as they are and only their labels follow the polymer.

What is done, in this order (see run()):
  1. H/D atoms and waters are removed. A model without alternative
     conformations is returned unchanged.
  2. Missing atoms of polymer residues that join atoms present in the model
     (CB between CA and CG, for example) are added from ideal geometry. Open
     ends are left.
  3. Every state of a polymer residue is made complete: atoms it lacks are
     shared with or copied from the other states
     (iotbx.pdb.altloc_completion, then the nearest state for what is
     left). The function stops (Sorry) if a state is too far from an atom
     it has to take from outside (a blank atom of the residue, or an atom
     that only other states have) to be bonded to it
     (bond_stretch_threshold).
  4. Every polymer residue is expanded to N atom groups, one per label:
       - a residue with N states starts with its states under their labels;
       - with fewer states, the N conformer slots are distributed evenly
         over the states (extra slots go to the states of highest
         occupancy) and the occupancy of a state is divided equally among
         its copies;
       - two or more labelled states of an atom adding up to more than 1
         are scaled down to 1; if they add up to less than 1 they are left
         as they are or, with fill_occupancy_to_one, the missing portion is
         added equally to each;
       - a blank atom, or an atom with a single label, of occupancy q
         becomes N copies of occupancy q/N.
  5. Untangling: copies of an atom exchange what they hold (position,
     B-factor, occupancy) between conformers wherever this lowers the
     covalent restraint target, so that each conformer has intact geometry.
  6. Groups of polymer atoms that are not tied to each other by restraints
     are arranged so that, in this order of priority, atoms of one conformer
     do not clash and occupancy rises from the first conformer to the last.
  7. Ligands and ions get the labels that fit the polymer as arranged: a
     labelled group moves to the label under which it does not clash and
     stays with the polymer state it was with in the input; a blank group
     of partial occupancy that clashes with some conformers gets the label
     of one where it does not (see arrangement.relabel_non_polymers).

Atoms of the input never move. Only copies at the same position (and the
atoms of step 2) are added, so structure factors are the same as those of
the input without H and water, except for the atoms of step 2 and the
polymer atoms whose occupancies were changed to add up to 1.

The function stops (Sorry) if a polymer atom without alternative
conformations and with partial occupancy clashes with another atom: such
alternates, encoded by occupancy only, would clash in every conformer.
"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb.hierarchy
import mmtbx.model
from cctbx.eltbx import chemical_elements
from libtbx import group_args
from libtbx.utils import Sorry, null_out

# Occupancy sums that differ from 1 by less than this are not reported.
OCCUPANCY_TOLERANCE = 0.005

def strip_h_and_water(model):
  """
  New model without H/D atoms and waters. The input model is not changed.
  """
  selection = model.selection("not (element H or element D or water)")
  return model.select(selection)

def merge_residue_groups(hierarchy):
  """
  Merge consecutive residue groups of a chain that have the same residue
  number and insertion code. Reading a file gives such groups when two
  residue types share one position (microheterogeneity) and some atoms have
  no altloc label. Returns the number of merges; the hierarchy is changed
  in place.
  """
  n_merged = 0
  for chain in hierarchy.chains():
    primary = None
    for rg in chain.residue_groups():
      if(primary is not None and primary.resid() == rg.resid()):
        chain.merge_residue_groups(primary=primary, secondary=rg)
        n_merged += 1
      else:
        primary = rg
  if(n_merged > 0): hierarchy.atoms().reset_serial()
  return n_merged

def conformer_labels(hierarchy):
  """
  Sorted unique non-blank altloc labels of the hierarchy.
  """
  result = set()
  for ag in hierarchy.atom_groups():
    if(ag.altloc.strip() != ""): result.add(ag.altloc)
  return sorted(result)

def slot_states(state_occupancies, n_slots):
  """
  Distribute n_slots conformer slots over the states of a residue.

  Returns a list of length n_slots with the index of the state each slot
  takes. Slots are distributed evenly, extra slots go to the states of
  highest occupancy, and states are taken in the order of ascending
  occupancy so that low occupancy lands in the first slots.
  """
  n_states = len(state_occupancies)
  assert 0 < n_states <= n_slots, (n_states, n_slots)
  order = sorted(range(n_states), key=lambda i: (state_occupancies[i], i))
  counts = [n_slots//n_states]*n_states
  for i in order[n_states-(n_slots % n_states):] if n_slots % n_states else []:
    counts[i] += 1
  result = []
  for i in order: result.extend([i]*counts[i])
  return result

def normalized_state_occupancies(occupancies, fill=True):
  """
  Occupancies of the labelled states of one atom, made to add up to 1.

  With two or more states: if the sum is more than 1 all are scaled down; if
  it is less than 1 the missing portion is added equally to each state, or,
  with fill=False, nothing is changed. A single state, or states that are
  all zero, are returned unchanged.
  """
  result = list(occupancies)
  n = len(result)
  total = sum(result)
  if(n < 2 or total <= 0): return result
  if(total < 1):
    if(not fill): return result
    return [o + (1-total)/n for o in result]
  return [o/total for o in result]

def _is_polymer(residue_group):
  """True if the residue is an amino acid or a nucleotide."""
  from iotbx.pdb import common_residue_names_get_class
  for resname in residue_group.unique_resnames():
    if(common_residue_names_get_class(resname) in ["common_amino_acid",
        "modified_amino_acid", "common_rna_dna", "modified_rna_dna"]):
      return True
  return False

class _state(object):
  """
  One state of a residue group: the atoms under one altloc label and residue
  name. atoms is a list of [name, atom, occupancy] in output order.
  """
  def __init__(self, altloc, resname):
    self.altloc = altloc
    self.resname = resname
    self.atoms = []
    self.own_names = set()

  def names(self):
    return [a[0] for a in self.atoms]

  def occupancy(self):
    values = [a[1].occ for a in self.atoms if a[0] in self.own_names]
    if(len(values) == 0): return 0
    return sum(values)/len(values)

def _distance(a, b):
  return ((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)**0.5

def _fill_gaps(states, filled):
  """
  Make all states of one residue name have the same atoms: an atom name that
  a state lacks is taken, at the same position, from the state whose atom is
  closest to the atoms the lacking state already has. The occupancy of the
  source atom is divided between its own state and the borrowers in
  proportion to the state occupancies (equally if these are all zero).
  """
  by_resname = {}
  for s in states: by_resname.setdefault(s.resname, []).append(s)
  for group in by_resname.values():
    if(len(group) < 2): continue
    names = []
    for s in group:
      for name in s.names():
        if(name not in names): names.append(name)
    state_occ = dict([(id(s), s.occupancy()) for s in group])
    for name in names:
      have = [s for s in group if name in s.names()]
      lack = [s for s in group if name not in s.names()]
      if(len(lack) == 0): continue
      borrowers = {} # id(source state) -> list of states
      for s in lack:
        best = None
        for p in have:
          atom = [a for a in p.atoms if a[0] == name][0][1]
          d = min([_distance(atom.xyz, a[1].xyz) for a in s.atoms])
          if(best is None or d < best[0]): best = (d, p, atom)
        borrowers.setdefault(id(best[1]), [best[1]]).append(s)
      for item in borrowers.values():
        source, targets = item[0], item[1:]
        entry = [a for a in source.atoms if a[0] == name][0]
        weights = [state_occ[id(s)] for s in [source]+targets]
        if(sum(weights) <= 0): weights = [1]*len(weights)
        occ = entry[2]
        entry[2] = occ*weights[0]/sum(weights)
        for s, w in zip(targets, weights[1:]):
          s.atoms.append([name, entry[1], occ*w/sum(weights), True])
          filled.append((entry[1].id_str(), source.altloc, s.altloc))

def _residue_group_states(rg):
  """
  States of a residue group and its blank atoms.

  Returns (states, blank): states is a list of _state for the labelled atom
  groups (same altloc and residue name merged); blank is a dict residue name
  -> list of blank atoms. Without labelled atom groups there is one state
  per residue name and blank is empty.
  """
  states, blank = [], {}
  for ag in rg.atom_groups():
    if(ag.altloc.strip() == ""):
      blank.setdefault(ag.resname, []).extend(list(ag.atoms()))
      continue
    s = [s for s in states if (s.altloc, s.resname) == (ag.altloc, ag.resname)]
    if(len(s) == 0):
      s = _state(altloc=ag.altloc, resname=ag.resname)
      states.append(s)
    else: s = s[0]
    for a in ag.atoms():
      s.atoms.append([a.name, a, a.occ])
      s.own_names.add(a.name)
  return states, blank

def _expand_residue_group(rg, new_rg, labels, changes, filled, origin,
                          is_filled, fill_occupancy_to_one):
  n = len(labels)
  states, blank = _residue_group_states(rg)
  def add(ag, atom, occ, label="", fill=False):
    c = atom.detached_copy()
    c.set_occ(occ)
    ag.append_atom(c)
    origin.append(label)
    is_filled.append(fill)
  if(len(states) == 0):
    # no alternative conformations: every residue name is copied to all slots
    for label in labels:
      for resname in blank:
        ag = iotbx.pdb.hierarchy.atom_group(altloc=label, resname=resname)
        new_rg.append_atom_group(ag)
        for a in blank[resname]: add(ag, a, a.occ/n)
    return
  if(len(states) > n):
    raise Sorry("Residue %s has more states (%s) than the model has "
      "alternative conformation labels (%s)." % (
      rg.atoms()[0].id_str()[10:-1], ", ".join(["%s %s" % (s.altloc,
      s.resname) for s in states]), ", ".join(labels)))
  # Occupancies of the labelled states add up to 1. An atom that all states
  # have is treated by itself; an atom that only some states have (another
  # residue type, or a truncated state) follows the occupancy of its state.
  if(len(states) > 1):
    raw = [s.occupancy() for s in states]
    target = normalized_state_occupancies(raw, fill=fill_occupancy_to_one)
    names = []
    for s in states:
      for name in s.names():
        if(name not in names): names.append(name)
    for name in names:
      entries = [(i, a) for i, s in enumerate(states)
        for a in s.atoms if a[0] == name]
      old = [e[2] for i, e in entries]
      if(len(entries) == len(states)):
        new = normalized_state_occupancies(old, fill=fill_occupancy_to_one)
      elif(sum(raw) < 1):
        new = [e[2] + target[i] - raw[i] for i, e in entries]
      else:
        new = [e[2]*target[i]/raw[i] if raw[i] > 0 else e[2]
          for i, e in entries]
      for (i, e), o in zip(entries, new): e[2] = o
      if(abs(sum(old)-sum(new)) > OCCUPANCY_TOLERANCE):
        changes.append((entries[0][1][1].id_str(), sum(old), sum(new)))
  _fill_gaps(states, filled)
  if(sorted([s.altloc for s in states]) == list(labels)):
    # every label has its state: start from the labels of the input
    slots = sorted(range(n), key=lambda i: (states[i].altloc, i))
  else:
    slots = slot_states([s.occupancy() for s in states], n)
  n_copies = [slots.count(i) for i in range(len(states))]
  # Blank atoms belong to every state that does not have the atom itself,
  # whatever residue name the blank atom group carries (microheterogeneity).
  blank_atoms, blank_copies = [], {}
  for resname in blank:
    for a in blank[resname]:
      if(a.name in blank_copies): continue
      blank_atoms.append(a)
      blank_copies[a.name] = len([i for i in slots
        if a.name not in states[i].names()])
  for label, i in zip(labels, slots):
    s = states[i]
    ag = iotbx.pdb.hierarchy.atom_group(altloc=label, resname=s.resname)
    new_rg.append_atom_group(ag)
    names = s.names()
    for a in blank_atoms:
      if(a.name not in names):
        add(ag, a, a.occ/blank_copies[a.name], fill=True)
    for e in s.atoms:
      add(ag, e[1], e[2]/n_copies[i], label=s.altloc, fill=len(e) > 3)

def expand(hierarchy, labels, completed=None, fill_occupancy_to_one=True):
  """
  Expand every polymer residue group to len(labels) atom groups, one per
  label. Other residue groups (ligands, ions) are copied as they are.

  completed, if given, is a hierarchy with the polymer residue groups only,
  in the same order, to take the polymer residues from (the result of
  altloc completion). fill_occupancy_to_one: see run().

  Returns group_args with: hierarchy, the new hierarchy, in which polymer
  atoms have no blank altlocs (microheterogeneity: the atom group of a slot
  has the residue name of the state in that slot); occupancy_changes, a list
  of (atom id_str, old sum, new sum) for atoms whose labelled states did not
  add up to 1; filled, a list of (atom id_str, source label, target label)
  for atoms copied into a state that lacked them; origin, for every atom of
  the new hierarchy the input label of the state it comes from ("" for a
  blank atom, None for an atom that is not of a polymer); is_filled, for
  every atom True if a state got it from outside: one of the copies listed
  in filled, or a blank atom of a residue that has states. Occupancies
  follow the rules in the module docstring. The input hierarchy is not
  changed.
  """
  changes, filled, origin, is_filled = [], [], [], []
  polymer = None
  if(completed is not None): polymer = iter(completed.residue_groups())
  root = iotbx.pdb.hierarchy.root()
  for model in hierarchy.models():
    new_model = iotbx.pdb.hierarchy.model(id=model.id)
    root.append_model(new_model)
    for chain in model.chains():
      new_chain = iotbx.pdb.hierarchy.chain(id=chain.id)
      new_model.append_chain(new_chain)
      for rg in chain.residue_groups():
        if(not _is_polymer(rg)):
          new_chain.append_residue_group(rg.detached_copy())
          origin.extend([None]*rg.atoms_size())
          is_filled.extend([False]*rg.atoms_size())
          continue
        if(polymer is not None):
          source = next(polymer)
          assert source.resid() == rg.resid()
          rg = source
        new_rg = iotbx.pdb.hierarchy.residue_group(
          resseq=rg.resseq, icode=rg.icode)
        new_chain.append_residue_group(new_rg)
        _expand_residue_group(rg, new_rg, labels, changes, filled, origin,
          is_filled, fill_occupancy_to_one)
  root.atoms().reset_serial()
  assert len(origin) == root.atoms_size()
  # making restraints sorts the atoms of a residue: number them, see
  # _in_model_order()
  for i, a in enumerate(root.atoms()): a.tmp = i
  return group_args(hierarchy=root, occupancy_changes=changes, filled=filled,
    origin=origin, is_filled=is_filled)

def _in_model_order(model, values):
  """
  Per-atom values given in the order of the atoms of expand() in the order
  of the atoms of the model made from its hierarchy.
  """
  numbers = model.get_hierarchy().atoms().extract_tmp_as_size_t()
  assert sorted(numbers) == list(range(len(values)))
  return [values[i] for i in numbers]

def _model_from_hierarchy(hierarchy, template):
  """
  New model manager for a hierarchy, with the crystal symmetry and restraint
  (CIF) objects of the template model.
  """
  return mmtbx.model.manager(
    model_input       = None,
    pdb_hierarchy     = hierarchy,
    crystal_symmetry  = template.crystal_symmetry(),
    restraint_objects = template.get_restraint_objects(),
    monomer_parameters= template.get_monomer_parameters(),
    log               = null_out())

def _process(model):
  """
  Make geometry restraints. Every residue must be known: a ligand that is
  not in the monomer library needs a restraints (CIF) file.
  """
  try:
    model.process(make_restraints=True)
  except Sorry as e:
    if("nonbonded energy type" not in str(e)): raise
    raise Sorry("%s\nRestraints are needed for every residue: provide a "
      "restraints (CIF) file for each ligand that is not in the monomer "
      "library." % str(e))

# ---------------------------------------------------------------- bridges

def _ideal_bonds(mon_lib_srv, resname, cache):
  """
  Heavy-atom names and bonds of a residue from the monomer library:
  (set of names, dict name -> set of bonded names, dict name -> element) or
  None if the residue is not in the library.
  """
  key = resname.strip().upper()
  if(key not in cache):
    comp = mon_lib_srv.get_comp_comp_id_direct(comp_id=key)
    if(comp is None):
      cache[key] = None
    else:
      element = {}
      for a in comp.atom_list:
        if(a.type_symbol.strip().upper() not in ["H","D"]):
          element[a.atom_id.strip()] = a.type_symbol.strip().upper()
      bonds = dict([(n, set()) for n in element])
      for b in comp.bond_list:
        n1, n2 = b.atom_id_1.strip(), b.atom_id_2.strip()
        if(n1 in element and n2 in element):
          bonds[n1].add(n2)
          bonds[n2].add(n1)
      cache[key] = (set(element.keys()), bonds, element)
  return cache[key]

# Atoms that join consecutive residues of a polymer: (name in this residue,
# name in the neighbour, +1 for the next residue or -1 for the previous one).
_POLYMER_LINKS = [("N","C",-1), ("C","N",1), ("P","O3'",-1), ("O3'","P",1)]

def _has_polymer_neighbour(rgs, i, name, other, step):
  j = i + step
  if(j < 0 or j >= len(rgs)): return False
  if(abs(rgs[j].resseq_as_int()-rgs[i].resseq_as_int()) > 1): return False
  return other in [a.name.strip() for a in rgs[j].atoms()]

def find_bridge_atoms(model):
  """
  Missing atoms that join atoms present in the model.

  A heavy atom of a residue known to the monomer library is a bridge if no
  atom group of the residue has it and, in the bond graph of the library
  definition, it lies on a path between two atoms that are present (an atom
  of the neighbouring residue of a polymer counts). Missing open-end atoms
  and missing residues are not bridges.

  Returns a list of group_args with: residue_group, resname, altloc (blank
  if all atoms the bridge is bonded to are blank, otherwise one entry per
  label among them), names (the bridge atoms of one connected gap), bonds
  and elements (library bond graph and elements of the residue).
  """
  mon_lib_srv = model.get_mon_lib_srv()
  cache = {}
  result = []
  for chain in model.get_hierarchy().only_model().chains():
    rgs = chain.residue_groups()
    for i_rg, rg in enumerate(rgs):
      if(not _is_polymer(rg)): continue # ligands stay as they are
      for resname in rg.unique_resnames():
        ideal = _ideal_bonds(mon_lib_srv, resname, cache)
        if(ideal is None): continue
        expected, bonds, elements = ideal
        ags = [ag for ag in rg.atom_groups()
          if ag.resname == resname or ag.altloc.strip() == ""]
        present = set()
        for ag in ags:
          for a in ag.atoms(): present.add(a.name.strip())
        missing = expected - present
        if(len(missing) == 0): continue
        # atoms held from outside the residue
        anchored = set()
        for name, other, step in _POLYMER_LINKS:
          if(name in expected and
             _has_polymer_neighbour(rgs, i_rg, name, other, step)):
            anchored.add(name)
        # remove missing atoms that lead nowhere until none is left
        bridges = set(missing)
        while True:
          leaves = []
          for name in bridges:
            n = len([o for o in bonds[name] if o in present or o in bridges])
            if(name in anchored): n += 1
            if(n < 2): leaves.append(name)
          if(len(leaves) == 0): break
          bridges -= set(leaves)
        # connected gaps
        while len(bridges) > 0:
          gap, todo = set(), [sorted(bridges)[0]]
          while len(todo) > 0:
            name = todo.pop()
            if(name in gap): continue
            gap.add(name)
            todo.extend([o for o in bonds[name] if o in bridges])
          bridges -= gap
          neighbours = set()
          for name in gap: neighbours |= (bonds[name] & present)
          # a gap held at one place only (a ring on one atom) is an open end
          if(len(neighbours) + len(gap & anchored) < 2): continue
          altlocs = set()
          for ag in ags:
            for a in ag.atoms():
              if(a.name.strip() in neighbours): altlocs.add(ag.altloc)
          labelled = sorted([l for l in altlocs if l.strip() != ""])
          if(len(labelled) == 0): labelled = [""]
          for altloc in labelled:
            result.append(group_args(
              residue_group = rg,
              resname       = resname,
              altloc        = altloc,
              names         = sorted(gap),
              bonds         = bonds,
              elements      = elements))
  return result

def _ideal_sites(resname, cache):
  """
  Ideal coordinates of an amino acid residue: dict atom name -> (xyz,
  four-character name), or None for other residues.
  """
  key = resname.strip().lower()
  if(key not in cache):
    from mmtbx.monomer_library import idealized_aa
    cache["dict"] = cache.get("dict") or idealized_aa.residue_dict()
    h = cache["dict"].get(key)
    if(h is None): cache[key] = None
    else:
      cache[key] = dict([(a.name.strip(), (a.xyz, a.name)) for a in h.atoms()])
  return cache[key]

def _place_bridge(bridge, ideal_cache):
  """
  Add the atoms of one gap to the hierarchy, at approximate positions:
  from the ideal residue superposed on the nearby present atoms if ideal
  coordinates are available, otherwise next to the atoms they are bonded to.
  B-factor is the mean and occupancy the smallest of the bonded atoms.
  Returns the new atoms.
  """
  from scitbx.array_family import flex
  from scitbx.math import superpose
  from scitbx import matrix
  rg = bridge.residue_group
  ag_target = [ag for ag in rg.atom_groups()
    if ag.altloc.strip() == bridge.altloc.strip()
    and (ag.resname == bridge.resname or ag.altloc.strip() == "")][0]
  # atoms of this conformer: blank and own label; own label wins
  conformer = {}
  for ag in rg.atom_groups():
    if(ag.altloc.strip() == "" or
       (ag.altloc == bridge.altloc and ag.resname == bridge.resname)):
      for a in ag.atoms():
        if(a.name.strip() not in conformer or ag.altloc == bridge.altloc):
          conformer[a.name.strip()] = a
  bonds = bridge.bonds
  # present atoms up to three bonds away from the gap
  near, shell = set(bridge.names), set(bridge.names)
  for i in range(3):
    shell = set([o for n in shell for o in bonds[n]]) - near
    near |= shell
  ideal = _ideal_sites(bridge.resname, ideal_cache)
  fit = None
  if(ideal is not None and all([n in ideal for n in bridge.names])):
    for names in [near, set(conformer.keys())]:
      names = sorted([n for n in names if n in conformer and n in ideal])
      if(len(names) >= 3):
        fit = superpose.least_squares_fit(
          reference_sites = flex.vec3_double([conformer[n].xyz for n in names]),
          other_sites     = flex.vec3_double([ideal[n][0] for n in names]))
        break
  new_atoms = []
  for i_name, name in enumerate(bridge.names):
    bonded = [conformer[o] for o in bonds[name] if o in conformer]
    if(len(bonded) == 0): bonded = list(conformer.values())
    if(fit is not None):
      xyz = (fit.r * matrix.col(ideal[name][0]) + fit.t).elems
      full_name = ideal[name][1]
    else:
      # next to the bonded atoms, off the line that joins them
      c = matrix.col((0,0,0))
      for a in bonded: c += matrix.col(a.xyz)
      c /= len(bonded)
      xyz = (c + matrix.col((0.5, 0.3*(i_name+1), 0.1))).elems
      element = bridge.elements[name]
      full_name = "%-4s" % name
      if(len(element) == 1 and len(name) < 4): full_name = " %-3s" % name
    atom = iotbx.pdb.hierarchy.atom()
    atom.set_name(full_name[:4])
    atom.set_element("%2s" % bridge.elements[name])
    atom.set_xyz(xyz)
    atom.set_b(sum([a.b for a in bonded])/len(bonded))
    atom.set_occ(min([a.occ for a in bonded]))
    atom.set_hetero(bonded[0].hetero)
    atom.set_segid(bonded[0].segid)
    ag_target.append_atom(atom)
    new_atoms.append(atom)
  return new_atoms

def add_bridge_atoms(model, log=None):
  """
  Add missing bridge atoms (see find_bridge_atoms) to a model.

  New atoms are placed from ideal geometry and then only they are minimized
  against the geometry restraints; all other atoms stay where they are.

  Returns (model, added): a new processed model and the id_str of the added
  atoms, or the input model and an empty list if there is nothing to add.
  """
  import scitbx.lbfgs
  from cctbx import geometry_restraints
  from mmtbx.refinement import geometry_minimization
  bridges = find_bridge_atoms(model)
  if(len(bridges) == 0): return model, []
  hierarchy = model.get_hierarchy()
  for a in hierarchy.atoms(): a.tmp = 0
  ideal_cache = {}
  n_placed = 0
  for bridge in bridges:
    for a in _place_bridge(bridge, ideal_cache):
      a.tmp = 1
      n_placed += 1
  hierarchy.atoms().reset_serial()
  result = _model_from_hierarchy(hierarchy, template=model)
  _process(result)
  atoms = result.get_hierarchy().atoms()
  selection = atoms.extract_tmp_as_size_t() == 1
  assert selection.count(True) == n_placed, (selection.count(True), n_placed)
  sites_cart = result.get_sites_cart()
  sites_start = sites_cart.deep_copy()
  for cycle in range(3):
    geometry_minimization.lbfgs(
      sites_cart                         = sites_cart,
      geometry_restraints_manager        =
        result.get_restraints_manager().geometry,
      geometry_restraints_flags          =
        geometry_restraints.flags.flags(default=True),
      lbfgs_termination_params           =
        scitbx.lbfgs.termination_parameters(max_iterations=500),
      correct_special_position_tolerance = 1.0,
      sites_cart_selection               = selection)
  # atoms of the input must not have moved
  sites_cart = sites_start.set_selected(selection, sites_cart.select(selection))
  result.set_sites_cart(sites_cart)
  added = [a.id_str() for a in atoms.select(selection)]
  for line in added: print("  bridge atom added:", line, file=log)
  return result, added

# -------------------------------------------------------------- relabelling

# Contacts closer than this (A) are used to choose between arrangements of
# equal clash count and occupancy order.
SOFT_CONTACT_DISTANCE = 3.0
# Copies of an atom further apart than this (A) are different positions.
SAME_POSITION_TOLERANCE = 1.e-4
# A move must lower the restraint target by more than this to be accepted.
# The target is a sum of (deviation/sigma)**2: a wrong pairing costs tens to
# thousands, rounding of coordinates in a file changes it by about 0.1.
TARGET_TOLERANCE = 1.0
# Restraints are made for the arrangement at hand: which atoms are linked,
# and some ideal values, depend on it. So after untangling restraints are
# made again and the new arrangement is kept only if its own target is
# lower. This is repeated until nothing improves, at most this many times.
MAX_CYCLES = 5
# A bond known from one conformer but absent from the restraints of another
# (its atoms are too far apart there to be linked) enters the target of the
# other conformer up to this deviation from ideal (A). Beyond it the term is
# constant, so a bond that no arrangement can make has no influence.
MISSING_BOND_CAP = 1.0
# A metal and the atoms it coordinates are closer than other non-bonded
# atoms and are not always linked in the restraints: a pair with a metal
# clashes only below this distance (A), or below clash_distance if smaller.
METAL_CLASH_DISTANCE = 1.5
METALS = chemical_elements.metals_upper_set()
# Bonds between consecutive residues of a polymer, expected whether or not
# the restraints have them: (atom name, atom name in the next residue,
# ideal distance, sigma).
POLYMER_BONDS = [(" C  ", " N  ", 1.33, 0.014), (" O3'", " P  ", 1.607, 0.012)]

def _proxy_i_seqs(geometry, sites_cart):
  """
  Every covalent restraint: generator of (kind, i_seqs, proxy). Bonds to
  symmetry copies have kind "bond_asu".
  """
  bond_proxies = geometry.pair_proxies(sites_cart=sites_cart).bond_proxies
  for p in bond_proxies.simple: yield "bond", p.i_seqs, p
  if(bond_proxies.asu is not None):
    for p in bond_proxies.asu: yield "bond_asu", (p.i_seq, p.j_seq), p
  for kind, proxies in [
      ("angle",       geometry.angle_proxies),
      ("dihedral",    geometry.dihedral_proxies),
      ("chirality",   geometry.chirality_proxies),
      ("planarity",   geometry.planarity_proxies),
      ("parallelity", geometry.parallelity_proxies)]:
    if(proxies is None): continue
    for p in proxies:
      if(kind == "parallelity"):
        yield kind, tuple(p.i_seqs)+tuple(p.j_seqs), p
      else:
        yield kind, p.i_seqs, p

class arrangement(object):
  """
  Which input position sits in which conformer slot.

  The split model has, for every atom (a key: residue group, residue name,
  atom name), one copy per conformer slot. Copies of a key can exchange
  what they hold (position, B-factor, occupancy); nothing else ever changes,
  so structure factors are the same for every arrangement.

  A key is split if its copies are not all at the same position. Split keys
  that occur together in a covalent restraint (bond, angle, dihedral,
  chirality, planarity) belong to one component: inside a component the
  restraint target depends on the arrangement, between components it does
  not.
  """

  def __init__(self, model, labels, origin=None):
    self.model = model
    # input label of the state each atom comes from, see expand()
    self.origin = origin
    self.labels = labels
    self.n_slots = len(labels)
    hierarchy = model.get_hierarchy()
    atoms = hierarchy.atoms()
    self.n_atoms = atoms.size()
    self.geometry = model.get_restraints_manager().geometry
    # what the atoms of the split model hold at the start
    self.sites_start = atoms.extract_xyz()
    self.occ_start = atoms.extract_occ()
    self.is_metal = [e.strip().upper() in METALS
      for e in atoms.extract_element()]
    # source[i]: index of the starting atom whose content atom i holds now
    self.source = list(range(self.n_atoms))
    self.sites = self.sites_start.deep_copy()
    # keys
    slot_of_label = dict([(l, i) for i, l in enumerate(labels)])
    self.key_of = [None]*self.n_atoms
    self.slot_of = [None]*self.n_atoms
    self.copies = []     # per key: atom index per slot, None if absent
    self.key_rg = []     # per key: index of its residue group
    self.key_polymer = [] # per key: True if of a polymer residue
    self.rg_of = [None]*self.n_atoms
    self.is_polymer_atom = [None]*self.n_atoms
    index = {}
    i_rg = -1
    self.next_rg = {}    # residue group index -> index of the next in chain
    previous = None
    for rg in hierarchy.residue_groups():
      i_rg += 1
      if(previous is not None and
         previous.parent().memory_id() == rg.parent().memory_id() and
         abs(rg.resseq_as_int()-previous.resseq_as_int()) <= 1 and
         _is_polymer(previous) and _is_polymer(rg)):
        self.next_rg[i_rg-1] = i_rg
      previous = rg
      polymer = _is_polymer(rg)
      for ag in rg.atom_groups():
        # blank (ligands and ions only): in every conformer, slot None
        slot = slot_of_label.get(ag.altloc)
        assert slot is not None or not polymer, ag.id_str()
        for a in ag.atoms():
          k = (i_rg, ag.resname, a.name)
          if(k not in index):
            index[k] = len(self.copies)
            self.copies.append([None]*self.n_slots)
            self.key_rg.append(i_rg)
            self.key_polymer.append(polymer)
          i_key = index[k]
          if(slot is not None):
            assert self.copies[i_key][slot] is None, a.id_str()
            self.copies[i_key][slot] = a.i_seq
          self.key_of[a.i_seq] = i_key
          self.slot_of[a.i_seq] = slot
          self.rg_of[a.i_seq] = i_rg
          self.is_polymer_atom[a.i_seq] = polymer
    self.n_keys = len(self.copies)
    self.key_index = index
    # only polymer atoms are exchanged between conformers
    self.is_split = [p and self._differ(c)
      for p, c in zip(self.key_polymer, self.copies)]
    # slot of every atom at the end: changes for ligands and ions only
    self.new_slot = list(self.slot_of)
    self._contacts_cache = None
    self._set_components()
    self.n_moves = {"cut":0, "residue":0, "piece":0, "atom":0}
    self.n_clash_permutations = 0

  def _differ(self, copies):
    sites = [self.sites_start[i] for i in copies if i is not None]
    for s in sites[1:]:
      if(_distance(s, sites[0]) > SAME_POSITION_TOLERANCE): return True
    return False

  def _set_components(self):
    """
    Components of split keys, key pairs related by a restraint, bonds
    between split keys, and the atoms each component needs for its target.
    """
    parent = list(range(self.n_keys))
    def find(i):
      while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
      return i
    self.related = set()   # pairs of keys in one covalent restraint
    self.key_bonds = {}    # split key -> set of bonded split keys
    split_proxies = []     # atom indices of restraints with a split key
    bond_slots = {}        # pair of keys -> slots in which they are bonded
    bond_params = {}       # pair of keys -> (ideal distance, weight)
    bond_proxies = self.geometry.pair_proxies(
      sites_cart=self.sites_start).bond_proxies
    self.asu_mappings = None
    if(bond_proxies.asu is not None and bond_proxies.asu.size() > 0):
      self.asu_mappings = bond_proxies.asu_mappings()
    # restraints with a split key: [kind, proxy, slot]; only these can change
    self.records = []
    self.records_of_key = {} # key -> indices in self.records
    # bonds of ligands and ions: to polymer atoms (ligand atom, polymer
    # atom) and between residues (residue group, residue group)
    self.np_links = []
    self.np_rg_links = []
    for kind, i_seqs, proxy in _proxy_i_seqs(self.geometry, self.sites_start):
      keys = sorted(set([self.key_of[i] for i in i_seqs]))
      # the conformer of a restraint: that of its first labelled atom
      slot = None
      for i in i_seqs:
        if(self.slot_of[i] is not None):
          slot = self.slot_of[i]
          break
      if(kind in ["bond", "bond_asu"] and len(keys) == 2):
        i, j = i_seqs
        if(self.is_polymer_atom[i] != self.is_polymer_atom[j]):
          if(self.is_polymer_atom[i]): i, j = j, i
          self.np_links.append((i, j))
        elif(not self.is_polymer_atom[i] and self.rg_of[i] != self.rg_of[j]):
          self.np_rg_links.append((self.rg_of[i], self.rg_of[j]))
        bond_slots.setdefault(tuple(keys), set()).add(slot)
        if(kind == "bond"):
          bond_params[tuple(keys)] = (proxy.distance_ideal, proxy.weight)
      for i, k1 in enumerate(keys):
        for k2 in keys[i+1:]: self.related.add((k1, k2))
      split = [k for k in keys if self.is_split[k]]
      if(len(split) == 0): continue
      split_proxies.append((split[0], i_seqs))
      self.records.append((kind, proxy, slot))
      for k in split:
        self.records_of_key.setdefault(k, []).append(len(self.records)-1)
      for k in split[1:]: parent[find(k)] = find(split[0])
      if(kind in ["bond", "bond_asu"] and len(split) == 2):
        self.key_bonds.setdefault(split[0], set()).add(split[1])
        self.key_bonds.setdefault(split[1], set()).add(split[0])
    # bonds between consecutive residues that the restraints do not have in
    # any conformer (atoms too far apart to be linked)
    by_rg_name = {}
    for (i_rg, resname, name), k in self.key_index.items():
      by_rg_name.setdefault((i_rg, name), []).append(k)
    for (i_rg, resname, name), k1 in sorted(self.key_index.items()):
      j_rg = self.next_rg.get(i_rg)
      if(j_rg is None): continue
      for name_1, name_2, ideal, sigma in POLYMER_BONDS:
        if(name != name_1): continue
        for k2 in sorted(by_rg_name.get((j_rg, name_2), [])):
          pair = (min(k1,k2), max(k1,k2))
          if(pair in bond_slots): continue
          bond_slots[pair] = set()
          bond_params[pair] = (ideal, 1./sigma**2)
          self.related.add(pair)
          if(self.is_split[k1] and self.is_split[k2]):
            parent[find(k1)] = find(k2)
            self.key_bonds.setdefault(k1, set()).add(k2)
            self.key_bonds.setdefault(k2, set()).add(k1)
    roots = {}
    self.components = []       # list of lists of keys
    self.component_of = [None]*self.n_keys
    for k in range(self.n_keys):
      if(not self.is_split[k]): continue
      r = find(k)
      if(r not in roots):
        roots[r] = len(self.components)
        self.components.append([])
      self.component_of[k] = roots[r]
      self.components[roots[r]].append(k)
    self.component_atoms = [set() for c in self.components]
    for c, keys in enumerate(self.components):
      for k in keys:
        self.component_atoms[c].update(
          [i for i in self.copies[k] if i is not None])
    for k, i_seqs in split_proxies:
      self.component_atoms[self.component_of[k]].update(i_seqs)
    # bonds that some conformers have and others lack
    self.missing_bonds = [] # (key, key, slot, ideal distance, weight)
    self.component_missing_bonds = [[] for c in self.components]
    for pair in sorted(bond_params.keys()):
      k1, k2 = pair
      for slot in range(self.n_slots):
        if(slot in bond_slots[pair]): continue
        if(self.copies[k1][slot] is None or self.copies[k2][slot] is None):
          continue
        item = (k1, k2, slot) + bond_params[pair]
        self.missing_bonds.append(item)
        if(self.is_split[k1] or self.is_split[k2]):
          self.records.append(("missing", item, slot))
          for k in pair:
            if(self.is_split[k]):
              self.records_of_key.setdefault(k, []).append(
                len(self.records)-1)
        for k in pair:
          c = self.component_of[k]
          if(c is None): continue
          self.component_missing_bonds[c].append(item)
          for kk in pair:
            self.component_atoms[c].update(
              [i for i in self.copies[kk] if i is not None])
          break

  def broken_bonds(self):
    """
    Bonds that the restraints have in some conformers but not in others:
    list of (atom index, atom index, distance) for the conformers without.
    """
    result = []
    for k1, k2, slot, ideal, weight in self.missing_bonds:
      i, j = self.copies[k1][slot], self.copies[k2][slot]
      result.append((i, j, _distance(self.sites[i], self.sites[j])))
    return result

  def changed(self):
    """
    True if any atom holds something else than at the start or has another
    label.
    """
    for i, s in enumerate(self.source):
      if(i != s): return True
    return self.new_slot != self.slot_of

  def swap(self, key, p, q):
    """Exchange what the copies of a key in slots p and q hold."""
    i, j = self.copies[key][p], self.copies[key][q]
    self.source[i], self.source[j] = self.source[j], self.source[i]
    s = self.sites[i]
    self.sites[i] = self.sites[j]
    self.sites[j] = s

  def swappable(self, key, p, q):
    """True if slots p and q of a key exist and hold different positions."""
    i, j = self.copies[key][p], self.copies[key][q]
    if(i is None or j is None): return False
    return _distance(self.sites[i], self.sites[j]) > SAME_POSITION_TOLERANCE

  # ------------------------------------------------------------- untangle

  def _moves(self, keys):
    """
    Sets of keys of one component to try exchanging between two slots:
    all keys on one side of a bond ("cut": a chain tail or a side-chain
    branch), the keys of one residue, a part not bonded to the rest
    ("piece") and single keys. Largest first.
    """
    keys = set(keys)
    bonds = dict([(k, self.key_bonds.get(k, set()) & keys) for k in keys])
    seen = set()
    moves = []
    def add(kind, members):
      members = frozenset(members)
      if(len(members) == 0 or len(members) == len(keys)): return
      complement = frozenset(keys - members)
      if(members in seen or complement in seen): return
      seen.add(members)
      moves.append((kind, sorted(members)))
    def reach(start, blocked=None):
      result, todo = set(), [start]
      while len(todo) > 0:
        k = todo.pop()
        if(k in result): continue
        result.add(k)
        for o in bonds[k]:
          if(blocked is not None and (k, o) in blocked): continue
          todo.append(o)
      return result
    for k in sorted(keys):
      for o in sorted(bonds[k]):
        if(o < k): continue
        side = reach(k, blocked=set([(k, o), (o, k)]))
        if(o in side): continue # the bond is part of a ring
        if(len(side) > len(keys)-len(side)): side = keys - side
        add("cut", side)
    by_rg = {}
    for k in keys: by_rg.setdefault(self.key_rg[k], []).append(k)
    for members in by_rg.values(): add("residue", members)
    todo = set(keys)
    while len(todo) > 0:
      piece = reach(sorted(todo)[0])
      todo -= piece
      add("piece", piece)
    for k in sorted(keys): add("atom", [k])
    moves.sort(key=lambda m: -len(m[1]))
    return moves

  def _missing_bond_term(self, item):
    k1, k2, slot, ideal, weight = item
    i, j = self.copies[k1][slot], self.copies[k2][slot]
    delta = min(abs(_distance(self.sites[i], self.sites[j])-ideal),
      MISSING_BOND_CAP)
    return weight*delta**2

  def _residual(self, record):
    """Value of one restraint for the current arrangement."""
    from cctbx import geometry_restraints as gr
    kind, proxy, slot = record
    if(kind == "missing"): return self._missing_bond_term(proxy)
    if(kind == "bond_asu"):
      return gr.bond(sites_cart=self.sites, asu_mappings=self.asu_mappings,
        proxy=proxy).residual()
    restraint = {"bond":gr.bond, "angle":gr.angle, "dihedral":gr.dihedral,
      "chirality":gr.chirality, "planarity":gr.planarity,
      "parallelity":gr.parallelity}[kind]
    return restraint(sites_cart=self.sites, proxy=proxy).residual()

  def _exchange(self, keys, p, q):
    """
    Exchange what slots p and q hold for the keys and return the change of
    the restraint target. Only restraints of slots p and q that involve the
    keys can change, so only these are evaluated.
    """
    affected = set()
    for k in keys: affected.update(self.records_of_key.get(k, []))
    affected = [self.records[i] for i in affected
      if self.records[i][2] in (p, q)]
    before = sum([self._residual(r) for r in affected])
    for k in keys: self.swap(k, p, q)
    return sum([self._residual(r) for r in affected]) - before

  def exchange(self, keys, p, q):
    """
    As _exchange, for the keys whose slots p and q hold different positions.
    """
    return self._exchange([k for k in keys if self.swappable(k, p, q)], p, q)

  def component_target(self, i_component):
    """
    Covalent restraint target of the atoms of one component and the atoms
    restrained with them, for the current arrangement. Used to check
    _exchange(); untangling itself evaluates only what a move changes.
    """
    from cctbx import geometry_restraints
    from scitbx.array_family import flex
    iselection = flex.size_t(sorted(self.component_atoms[i_component]))
    geometry = self.geometry.select(iselection=iselection)
    flags = geometry_restraints.flags.flags(default=True, nonbonded=False)
    result = geometry.energies_sites(sites_cart=self.sites.select(iselection),
      flags=flags, compute_gradients=False).target
    for item in self.component_missing_bonds[i_component]:
      result += self._missing_bond_term(item)
    return result

  def untangle_component(self, i_component, max_sweeps=50):
    """
    Lower the covalent restraint target of one component by exchanging the
    content of two slots for sets of its keys, until no exchange helps.
    Returns the number of accepted moves by kind and the accepted moves as
    a list of (keys, slot, slot).
    """
    moves = self._moves(self.components[i_component])
    accepted = {"cut":0, "residue":0, "piece":0, "atom":0}
    swaps = []
    for sweep in range(max_sweeps):
      changed = False
      for kind, members in moves:
        for p in range(self.n_slots):
          for q in range(p+1, self.n_slots):
            todo = [k for k in members if self.swappable(k, p, q)]
            if(len(todo) == 0): continue
            if(self._exchange(todo, p, q) < -TARGET_TOLERANCE):
              accepted[kind] += 1
              swaps.append((todo, p, q))
              changed = True
            else:
              for k in todo: self.swap(k, p, q) # put back
      if(not changed): break
    return accepted, swaps

  def untangle(self, nproc=1):
    """
    Untangle every component. Components are independent of each other, so
    with nproc > 1 they are done in parallel processes.
    """
    n = len(self.components)
    if(nproc is None or nproc <= 1 or n < 2):
      for i_component in range(n):
        accepted, swaps = self.untangle_component(i_component)
        for kind in accepted: self.n_moves[kind] += accepted[kind]
      return
    from libtbx import easy_mp
    # largest first, for balance
    order = sorted(range(n), key=lambda c: -len(self.components[c]))
    results = easy_mp.pool_map(
      fixed_func = self.untangle_component,
      args       = order,
      processes  = min(nproc, n))
    # the processes worked on copies: repeat their moves here
    for accepted, swaps in results:
      for kind in accepted: self.n_moves[kind] += accepted[kind]
      for keys, p, q in swaps:
        for k in keys: self.swap(k, p, q)

  # ------------------------------------------------- clashes and occupancy

  def _contacts(self):
    """
    Pairs of starting atoms closer than SOFT_CONTACT_DISTANCE, symmetry
    copies included, that are not copies of one key and not in one covalent
    restraint: list of (i, j, distance).
    """
    if(self._contacts_cache is not None): return self._contacts_cache
    xrs = self.model.get_xray_structure()
    unit_cell = xrs.unit_cell()
    sites_frac = unit_cell.fractionalize(self.sites_start)
    pair_sym_table = xrs.pair_asu_table(
      distance_cutoff=SOFT_CONTACT_DISTANCE).extract_pair_sym_table()
    result = []
    for i, pair_sym_dict in enumerate(pair_sym_table):
      ki = self.key_of[i]
      for j, sym_ops in pair_sym_dict.items():
        kj = self.key_of[j]
        if(ki == kj): continue
        if((min(ki,kj), max(ki,kj)) in self.related): continue
        for rt_mx in sym_ops:
          d = unit_cell.distance(sites_frac[i], rt_mx*sites_frac[j])
          result.append((i, j, d))
    self._contacts_cache = result
    return result

  def clash_limit(self, i, j, clash_distance):
    """Distance below which atoms i and j clash."""
    if(self.is_metal[i] or self.is_metal[j]):
      return min(clash_distance, METAL_CLASH_DISTANCE)
    return clash_distance

  def _component_occupancies(self, i_component):
    """Mean occupancy of the keys of a component in each slot."""
    keys = self.components[i_component]
    result = []
    for slot in range(self.n_slots):
      values = [self.occ_start[self.source[self.copies[k][slot]]]
        for k in keys]
      result.append(sum(values)/len(values))
    return result

  def _permute(self, i_component, permutation):
    """Slot s of a component gets what its slot permutation[s] holds."""
    for k in self.components[i_component]:
      held = [(self.source[i], self.sites[i]) for i in self.copies[k]]
      for slot, i in enumerate(self.copies[k]):
        self.source[i], self.sites[i] = held[permutation[slot]]

  def _places(self):
    """place[i]: atom of the split model that holds starting atom i now."""
    result = [None]*self.n_atoms
    for i, s in enumerate(self.source): result[s] = i
    return result

  def order(self, clash_distance):
    """
    Choose the slot permutation of every component: fewest clashes first,
    then occupancy rising from the first slot to the last, then the least
    close contacts. Components that have keys missing in some slots
    (microheterogeneity) are left as they are. Ligands and ions do not
    count here: they follow the polymer, see relabel_non_polymers().
    """
    import itertools
    n = self.n_slots
    free = [c for c, keys in enumerate(self.components)
      if all([None not in self.copies[k] for k in keys])]
    # start from occupancy order
    for c in free:
      occ = self._component_occupancies(c)
      self._permute(c, sorted(range(n), key=lambda s: (round(occ[s],6), s)))
    contacts = self._contacts()
    place = self._places()
    by_component = dict([(c, []) for c in free])
    for i, j, d in contacts:
      if(not (self.is_polymer_atom[i] and self.is_polymer_atom[j])): continue
      ci = self.component_of[self.key_of[i]]
      cj = self.component_of[self.key_of[j]]
      if(ci == cj): continue
      if(ci in by_component): by_component[ci].append((i, j, d))
      if(cj in by_component): by_component[cj].append((j, i, d))
    if(n <= 4): permutations = list(itertools.permutations(range(n)))
    else:
      permutations = [tuple(range(n))]
      for p in range(n):
        for q in range(p+1, n):
          t = list(range(n))
          t[p], t[q] = q, p
          permutations.append(tuple(t))
    def disorder(occ):
      result = 0
      for p in range(n):
        for q in range(p+1, n): result += max(0, occ[p]-occ[q])
      return round(result/OCCUPANCY_TOLERANCE)
    def score(c, permutation, place):
      # a starting atom in slot s moves to the slot that takes from s
      to_slot = [None]*n
      for slot, s in enumerate(permutation): to_slot[s] = slot
      n_clash, soft = 0, 0
      for i, j, d in by_component[c]:
        if(to_slot[self.slot_of[place[i]]] != self.slot_of[place[j]]):
          continue
        if(d < self.clash_limit(i, j, clash_distance)): n_clash += 1
        soft += (SOFT_CONTACT_DISTANCE-d)**2
      occ = self._component_occupancies(c)
      return (n_clash, disorder([occ[s] for s in permutation]),
        round(soft, 6))
    identity = tuple(range(n))
    order = sorted(free, key=lambda c: (len(self.components[c]), c))
    for cycle in range(100):
      changed = False
      for c in order:
        if(len(by_component[c]) == 0): continue
        start = score(c, identity, place)
        best = min([(score(c, p, place), p) for p in permutations])
        if(best[0] < start):
          self._permute(c, best[1])
          for k in self.components[c]:
            for i in self.copies[k]: place[self.source[i]] = i
          if(best[0][0] < start[0]): self.n_clash_permutations += 1
          changed = True
      if(not changed): break

  def relabel_non_polymers(self, clash_distance):
    """
    Choose the labels of ligands and ions for the arrangement of the polymer
    at hand. Their atoms, occupancies and number of copies stay as they are;
    only the labels of their atom groups change, where this is needed:

      - a labelled group gets the label under which it clashes least with
        the atoms of that conformer; of such labels, the one that keeps it
        with the polymer atoms it is bonded to, then the one that keeps it
        with the polymer atoms in contact with it (closer than
        SOFT_CONTACT_DISTANCE) that had its label in the input; otherwise
        it keeps its label. Labelled groups of one residue, and of residues
        bonded to each other, are relabelled together and stay distinct;
      - a blank group with partial occupancy that clashes with atoms of
        some conformers gets the label of the conformer where it clashes
        least: these are alternates encoded by occupancy only.

    Blank groups with full occupancy stay blank. Returns a list of (atom
    index of the first atom of the group, old label, new label).
    """
    import itertools
    n = self.n_slots
    place = self._places()
    groups = {} # (residue group, slot) -> atoms
    for i in range(self.n_atoms):
      if(self.is_polymer_atom[i]): continue
      groups.setdefault((self.rg_of[i], self.slot_of[i]), []).append(i)
    if(len(groups) == 0): return []
    labelled = set([rg for rg, slot in groups if slot is not None])
    parent = dict([(rg, rg) for rg in labelled])
    def find(i):
      while parent[i] != i: i = parent[i]
      return i
    for r1, r2 in self.np_rg_links:
      if(r1 in parent and r2 in parent): parent[find(r1)] = find(r2)
    units = {} # unit id -> atoms as (atom, slot)
    for (rg, slot) in sorted(groups.keys(),
        key=lambda g: (g[0], -1 if g[1] is None else g[1])):
      if(slot is None):
        if(rg in labelled): continue
        if(max([self.occ_start[i] for i in groups[(rg, slot)]]) >=
           1 - OCCUPANCY_TOLERANCE): continue
        uid = (rg, "blank")
      else:
        uid = (find(rg), "labelled")
      units.setdefault(uid, []).extend([(i, slot) for i in groups[(rg, slot)]])
    if(len(units) == 0): return []
    contacts_of, links_of = {}, {}
    for i, j, d in self._contacts():
      for x, y in [(i, j), (j, i)]:
        if(not self.is_polymer_atom[x]):
          contacts_of.setdefault(x, []).append((y, d))
    for i, j in self.np_links: links_of.setdefault(i, []).append(j)
    def score(atoms, in_unit, m):
      n_clash, n_link, n_follow = 0, 0, 0
      for i, s in atoms:
        t = m[s]
        for j, d in contacts_of.get(i, []):
          polymer = self.is_polymer_atom[j]
          if(j in in_unit): tj = m[in_unit[j]]
          elif(polymer):    tj = self.slot_of[place[j]]
          else:             tj = self.new_slot[j]
          if(not (t is None or tj is None or t == tj)): continue
          if(d < self.clash_limit(i, j, clash_distance)): n_clash += 1
          elif(polymer and s is not None and self.origin is not None and
               self.origin[j] == self.labels[s]): n_follow += 1
        if(s is None): continue
        for j in links_of.get(i, []):
          if(self.slot_of[place[j]] == t): n_link += 1
      n_changed = len([s for s in m if m[s] != s])
      return (n_clash, -n_link, -n_follow, n_changed)
    for cycle in range(10):
      changed = False
      for uid in sorted(units.keys()):
        atoms = units[uid]
        in_unit = dict(atoms)
        slots = sorted(set([s for i, s in atoms if s is not None]))
        if(len(slots) == 0):
          candidates = [{None:None}] + [{None:t} for t in range(n)]
        else:
          candidates = [dict(zip(slots, p))
            for p in itertools.permutations(range(n), len(slots))]
        best = None
        for m in candidates:
          value = score(atoms, in_unit, m)
          if(best is None or value < best[0]): best = (value, m)
        for i, s in atoms:
          if(self.new_slot[i] != best[1][s]):
            self.new_slot[i] = best[1][s]
            changed = True
      if(not changed): break
    result = []
    def label(slot): return "" if slot is None else self.labels[slot]
    for key in sorted(groups.keys(),
        key=lambda g: (g[0], -1 if g[1] is None else g[1])):
      i = groups[key][0]
      if(self.new_slot[i] != self.slot_of[i]):
        result.append((i, label(self.slot_of[i]), label(self.new_slot[i])))
    return result

  def clashes(self, clash_distance):
    """
    Clashes of the arrangement and labels at hand: list of (i, j, distance)
    of atoms of the split model that are in the same conformer (a blank
    atom is in every conformer) and closer than the clash distance.
    """
    place = self._places()
    result = []
    seen = set()
    for i, j, d in self._contacts():
      if(d >= self.clash_limit(i, j, clash_distance)): continue
      pi, pj = place[i], place[j]
      si, sj = self.new_slot[pi], self.new_slot[pj]
      if(not (si is None or sj is None or si == sj)): continue
      # copies of unsplit atoms clash in every conformer: report once
      pair = (self.key_of[i], self.key_of[j], round(d, 3))
      if(not self.is_split[pair[0]] and not self.is_split[pair[1]]):
        if(pair in seen): continue
        seen.add(pair)
      result.append((pi, pj, d))
    return result

  def apply(self):
    """
    A copy of the model hierarchy with the arrangement written to its atoms
    (position, B-factor, occupancy) and the new labels of ligands and ions
    written to their atom groups. The order of atoms is that of the model.
    The model itself is not changed.
    """
    from scitbx.array_family import flex
    hierarchy = self.model.get_hierarchy().deep_copy()
    atoms = hierarchy.atoms()
    source = flex.size_t(self.source)
    atoms.set_xyz(atoms.extract_xyz().select(source))
    atoms.set_occ(atoms.extract_occ().select(source))
    atoms.set_b(atoms.extract_b().select(source))
    atoms.set_uij(atoms.extract_uij().select(source))
    i = 0
    for ag in hierarchy.atom_groups():
      if(ag.atoms_size() == 0): continue
      if(self.new_slot[i] != self.slot_of[i]):
        slot = self.new_slot[i]
        ag.altloc = "" if slot is None else self.labels[slot]
      i += ag.atoms_size()
    atoms.reset_serial()
    return hierarchy

  def total_target(self):
    """
    Covalent restraint target of the whole model as it is at the start,
    plus the capped terms of the bonds that some conformers lack.
    """
    from cctbx import geometry_restraints
    flags = geometry_restraints.flags.flags(default=True, nonbonded=False)
    result = self.geometry.energies_sites(sites_cart=self.sites_start,
      flags=flags, compute_gradients=False).target
    for k1, k2, slot, ideal, weight in self.missing_bonds:
      i, j = self.copies[k1][slot], self.copies[k2][slot]
      delta = min(abs(_distance(self.sites_start[i], self.sites_start[j])
        - ideal), MISSING_BOND_CAP)
      result += weight*delta**2
    return result

def conformer_geometry(model, labels):
  """
  (bond rmsd, angle rmsd) of every conformer taken alone, in label order.
  """
  from cctbx import geometry_restraints
  flags = geometry_restraints.flags.flags(default=True, nonbonded=False)
  geometry = model.get_restraints_manager().geometry
  sites_cart = model.get_sites_cart()
  result = []
  for label in labels:
    selection = model.selection("altloc ' ' or altloc '%s'" % label)
    e = geometry.select(selection=selection).energies_sites(
      sites_cart=sites_cart.select(selection), flags=flags,
      compute_gradients=False)
    result.append((e.bond_deviations()[2], e.angle_deviations()[2]))
  return result

def _stretched(model, bond_stretch_threshold):
  """
  Bonds longer than ideal by more than the threshold: list of (atom index,
  atom index, distance, ideal distance).
  """
  from cctbx import geometry_restraints
  geometry = model.get_restraints_manager().geometry
  sites_cart = model.get_sites_cart()
  bond_proxies = geometry.pair_proxies(sites_cart=sites_cart).bond_proxies
  result = []
  def check(i, j, bond):
    if(bond.distance_model - bond.distance_ideal > bond_stretch_threshold):
      result.append((i, j, bond.distance_model, bond.distance_ideal))
  for p in bond_proxies.simple:
    check(p.i_seqs[0], p.i_seqs[1],
      geometry_restraints.bond(sites_cart=sites_cart, proxy=p))
  if(bond_proxies.asu is not None and bond_proxies.asu.size() > 0):
    asu_mappings = bond_proxies.asu_mappings()
    for p in bond_proxies.asu:
      check(p.i_seq, p.j_seq, geometry_restraints.bond(
        sites_cart=sites_cart, asu_mappings=asu_mappings, proxy=p))
  return result

def stretched_bonds(model, bond_stretch_threshold):
  """
  Bonds longer than ideal by more than the threshold: list of (atom id_str,
  atom id_str, distance, ideal distance).
  """
  atoms = model.get_hierarchy().atoms()
  return [(atoms[i].id_str(), atoms[j].id_str(), d, d_ideal)
    for i, j, d, d_ideal in _stretched(model, bond_stretch_threshold)]

def _stop_on_truncated_states(model, is_filled, bond_stretch_threshold):
  """
  Stop if a state of a polymer residue is too far from an atom of the
  residue that it has to use together with the other states: an atom
  without a label, or one that this state lacks and the others have. Such
  a state is incomplete, and what it lacks is not built here.
  """
  atoms = model.get_hierarchy().atoms()
  lines = []
  for i, j, d, d_ideal in _stretched(model, bond_stretch_threshold):
    if(is_filled[i] == is_filled[j]): continue
    if(atoms[i].parent().parent().memory_id() !=
       atoms[j].parent().parent().memory_id()): continue
    if(is_filled[j]): i, j = j, i
    line = "  %s - %s: %.2f A (ideal %.2f)" % (atoms[i].id_str(),
      atoms[j].id_str(), d, d_ideal)
    if(line not in lines): lines.append(line)
  if(len(lines) > 0):
    raise Sorry("Some alternative conformations of residues are incomplete: "
      "they are too far from atoms that they have to share with the other "
      "conformations of the residue (atoms without an alternative "
      "conformation label, or atoms that only the other conformations "
      "have). Bonds to such atoms (listed first):\n%s\n"
      "This cannot be repaired automatically: complete or remove these "
      "alternative conformations." % "\n".join(lines))

def self_symmetry_clashes(model, id_strs, clash_distance):
  """
  Atoms, of those with the given id_str, that are closer than the clash
  distance to their own symmetry copy without being on a special position:
  list of (id_str, distance).
  """
  if(len(id_strs) == 0 or model.crystal_symmetry() is None): return []
  id_strs = set(id_strs)
  atoms = model.get_hierarchy().atoms()
  xrs = model.get_xray_structure()
  unit_cell = xrs.unit_cell()
  sites_frac = xrs.sites_frac()
  pair_sym_table = xrs.pair_asu_table(
    distance_cutoff=clash_distance).extract_pair_sym_table()
  result = []
  for i, pair_sym_dict in enumerate(pair_sym_table):
    if(i not in pair_sym_dict or atoms[i].id_str() not in id_strs): continue
    for rt_mx in pair_sym_dict[i]:
      d = unit_cell.distance(sites_frac[i], rt_mx*sites_frac[i])
      if(SAME_POSITION_TOLERANCE < d < clash_distance):
        result.append((atoms[i].id_str(), d))
  return result

def _sort_atom_groups(hierarchy):
  """
  Put the atom groups of every residue group in the order of their labels
  (after ligands were relabelled).
  """
  for rg in hierarchy.residue_groups():
    groups = list(rg.atom_groups())
    ordered = sorted(groups, key=lambda ag: ag.altloc.strip())
    if([ag.altloc for ag in groups] == [ag.altloc for ag in ordered]):
      continue
    for ag in groups: rg.remove_atom_group(ag)
    for ag in ordered: rg.append_atom_group(ag)
  hierarchy.atoms().reset_serial()

def _stop_on_unlabelled_alternates(arrangement, clashes, id_strs):
  """
  Stop if a polymer atom without alternative conformations and with partial
  occupancy clashes with another atom: these are alternates encoded by
  occupancy only, and copying them into every conformer makes clashes that
  no arrangement can avoid.
  """
  a = arrangement
  lines = []
  for i, j, d in clashes:
    for x in (i, j):
      k = a.key_of[x]
      if(not a.key_polymer[k]): continue
      if(a.is_split[k] or None in a.copies[k]): continue
      total = sum([a.occ_start[c] for c in a.copies[k] if c is not None])
      # each copy may come rounded from a file (3 x 0.33)
      if(total < 1 - OCCUPANCY_TOLERANCE*max(1, a.n_slots)):
        lines.append("  %s - %s: %.2f A" % (id_strs[i], id_strs[j], d))
        break
  if(len(lines) > 0):
    raise Sorry("Polymer atoms with partial occupancy and a single position "
      "(no alternative conformation label, or one label only) are copied into "
      "every conformer and clash there with other atoms:\n%s\n"
      "These look like alternates of each other: give them alternative "
      "conformation labels that say so." % "\n".join(lines))

def run(model, clash_distance=2.1, bond_stretch_threshold=0.5,
        fill_occupancy_to_one=False, untangle=True, nproc=1, log=None):
  """
  Split a model into complete conformers, see the module docstring.

  Parameters
  ----------
  model : mmtbx.model.manager
    Not changed.
  clash_distance : float
    Two non-bonded atoms of one conformer closer than this (A) clash.
  bond_stretch_threshold : float
    Bonds longer than ideal by more than this (A) are reported. Stop if
    a state of a residue makes such a bond with an atom it shares with the
    other states (an incomplete state).
  fill_occupancy_to_one : bool
    Labelled states of a polymer atom that add up to less than 1: if True
    the missing portion is added equally to each state (structure factors
    change), if False they are left as they are. Sums above 1 are scaled
    down to 1 in both cases.
  untangle : bool
    If False the positions stay in the conformers the expansion put them
    in: no exchange for geometry, clashes or occupancy order.
  nproc : int
    Number of processes for untangling.
  log : file-like or None

  Returns
  -------
  libtbx.group_args
    hierarchy: the new iotbx.pdb.hierarchy.root; model: the same as a
    processed mmtbx.model.manager; report: group_args, see show().
    If the model has no alternative conformations (after removing H and
    water) nothing is done: model is the input model and
    report.nothing_done is True.
  """
  if(log is None): log = null_out()
  import time
  t0 = [time.time()]
  def timing(what):
    print("  %-28s %7.1f s" % (what+":", time.time()-t0[0]), file=log)
    t0[0] = time.time()
  if(len(model.get_hierarchy().models()) != 1):
    raise Sorry("One model is required.")
  report = group_args(
    nothing_done      = False,
    labels            = [],
    n_h_removed       = 0,
    n_water_removed   = 0,
    occupancy_changes = [],
    filled            = [],
    bridge_atoms      = [],
    completion        = None,
    n_components      = 0,
    n_moves           = {},
    n_clash_permutations = 0,
    geometry_start    = [],
    geometry_final    = [],
    stretched_bonds   = [],
    broken_bonds      = [],
    clashes           = [],
    relabelled        = [],
    n_atoms           = [])
  stripped = strip_h_and_water(model)
  if(merge_residue_groups(stripped.get_hierarchy()) > 0):
    stripped = _model_from_hierarchy(stripped.get_hierarchy(), template=model)
  report.labels = conformer_labels(stripped.get_hierarchy())
  if(len(report.labels) == 0):
    report.nothing_done = True
    return group_args(hierarchy=model.get_hierarchy(), model=model,
      report=report)
  sel_h = model.selection("element H or element D")
  report.n_h_removed = sel_h.count(True)
  report.n_water_removed = (model.selection("water") & ~sel_h).count(True)
  stripped, report.bridge_atoms = add_bridge_atoms(stripped, log=log)
  timing("bridge atoms")
  bridge_clashes = self_symmetry_clashes(stripped, report.bridge_atoms,
    clash_distance)
  # polymer residues are made complete; ligands and ions stay as they are
  import iotbx.pdb.altloc_completion
  polymer = stripped.get_hierarchy().deep_copy()
  for chain in polymer.chains():
    for rg in chain.residue_groups():
      if(not _is_polymer(rg)): chain.remove_residue_group(rg)
  report.completion = iotbx.pdb.altloc_completion.run(hierarchy=polymer)
  expanded = expand(stripped.get_hierarchy(), report.labels,
    completed=polymer, fill_occupancy_to_one=fill_occupancy_to_one)
  report.occupancy_changes = expanded.occupancy_changes
  report.filled = expanded.filled
  result = _model_from_hierarchy(expanded.hierarchy, template=model)
  _process(result)
  timing("expand, make restraints")
  _stop_on_truncated_states(result,
    _in_model_order(result, expanded.is_filled), bond_stretch_threshold)
  report.geometry_start = conformer_geometry(result, report.labels)
  origin = _in_model_order(result, expanded.origin)
  a = arrangement(model=result, labels=report.labels, origin=origin)
  report.n_components = len(a.components)
  report.n_moves = dict(a.n_moves)
  clashes = []
  if(untangle):
    target = a.total_target()
    for cycle in range(MAX_CYCLES):
      a.untangle(nproc=nproc)
      timing("untangle (cycle %d)" % (cycle+1))
      if(sum(a.n_moves.values()) == 0): break
      candidate = _model_from_hierarchy(a.apply(), template=model)
      _process(candidate)
      timing("make restraints")
      origin_candidate = [origin[i] for i in a.source]
      a_candidate = arrangement(model=candidate, labels=report.labels,
        origin=origin_candidate)
      target_candidate = a_candidate.total_target()
      if(target_candidate > target - TARGET_TOLERANCE):
        # not better with its own restraints: keep what we had
        a = arrangement(model=result, labels=report.labels, origin=origin)
        break
      for kind in a.n_moves: report.n_moves[kind] += a.n_moves[kind]
      result, a, target = candidate, a_candidate, target_candidate
      origin = origin_candidate
    a.order(clash_distance=clash_distance)
    report.n_clash_permutations = a.n_clash_permutations
    relabelled = a.relabel_non_polymers(clash_distance=clash_distance)
    clashes = a.clashes(clash_distance=clash_distance)
    timing("clashes, occupancy order")
    id_start = [x.id_str() for x in result.get_hierarchy().atoms()]
    hierarchy = a.apply()
    id_strs = [x.id_str() for x in hierarchy.atoms()]
    report.relabelled = [(id_start[i], old, new) for i, old, new in relabelled]
    report.clashes = [(id_strs[i], id_strs[j], d) for i, j, d in clashes]
    _stop_on_unlabelled_alternates(a, clashes, id_strs)
    if(a.changed()):
      _sort_atom_groups(hierarchy)
      result = _model_from_hierarchy(hierarchy, template=model)
      _process(result)
      a = arrangement(model=result, labels=report.labels)
  for id_str, d in bridge_clashes:
    report.clashes.append((id_str, id_str + " (symmetry copy)", d))
  atoms = result.get_hierarchy().atoms()
  report.broken_bonds = [(atoms[i].id_str(), atoms[j].id_str(), d)
    for i, j, d in a.broken_bonds()]
  counts = {}
  for ag in result.get_hierarchy().atom_groups():
    key = ag.altloc.strip()
    counts[key] = counts.get(key, 0) + ag.atoms_size()
  report.n_atoms = [counts.get(l, 0) + counts.get("", 0)
    for l in report.labels]
  report.geometry_final = conformer_geometry(result, report.labels)
  report.stretched_bonds = stretched_bonds(result, bond_stretch_threshold)
  return group_args(hierarchy=result.get_hierarchy(), model=result,
    report=report)

def show(report, log, clash_distance=None, bond_stretch_threshold=None):
  """
  Print the report of run() to log.
  """
  def p(*args): print(*args, file=log)
  if(report.nothing_done):
    p("The model has no alternative conformations: nothing is done.")
    return
  p("Removed: %d H/D atoms, %d water atoms" % (report.n_h_removed,
    report.n_water_removed))
  p("Conformers: %d (%s)" % (len(report.labels), ", ".join(report.labels)))
  p("Atoms per conformer: %s" % ", ".join(["%s %d" % (l, n)
    for l, n in zip(report.labels, report.n_atoms)]))
  if(clash_distance is not None):
    p("Clash distance: %.2f A" % clash_distance)
  p("Bridge atoms added: %d" % len(report.bridge_atoms))
  for line in report.bridge_atoms: p("  %s" % line)
  n_completed = 0
  if(report.completion is not None):
    n_completed = len(report.completion.shared) + len(report.completion.copied)
  p("Atoms shared or copied between conformers of a residue: %d" %
    (n_completed + len(report.filled)))
  for id_str, source, target in report.filled:
    p("  %s: %s -> %s (state had no such atom)" % (id_str, source, target))
  p("Occupancy sums changed: %d" % len(report.occupancy_changes))
  for id_str, old, new in report.occupancy_changes:
    p("  %s: %.2f -> %.2f" % (id_str, old, new))
  p("Groups of atoms that can be exchanged between conformers: %d" %
    report.n_components)
  p("Exchanges made for geometry: %s" % ", ".join(["%s %d" % (k,
    report.n_moves[k]) for k in sorted(report.n_moves.keys())]))
  p("Groups rearranged to avoid clashes: %d" % report.n_clash_permutations)
  p("Ligands and ions relabelled: %d" % len(report.relabelled))
  for id_str, old, new in report.relabelled:
    p("  %s: '%s' -> '%s'" % (id_str, old, new))
  p("Geometry of each conformer, start -> final (rmsd bonds, angles):")
  for label, start, final in zip(report.labels, report.geometry_start,
                                 report.geometry_final):
    p("  %s: %.3f %.2f -> %.3f %.2f" % (label, start[0], start[1],
      final[0], final[1]))
  threshold = ""
  if(bond_stretch_threshold is not None):
    threshold = " (longer than ideal by more than %.2f A)" % \
      bond_stretch_threshold
  p("Stretched bonds%s: %d" % (threshold, len(report.stretched_bonds)))
  for id_1, id_2, d, d_ideal in report.stretched_bonds:
    p("  %s - %s: %.2f A (ideal %.2f)" % (id_1, id_2, d, d_ideal))
  p("Broken bonds (bonded in another conformer or consecutive residues): %d"
    % len(report.broken_bonds))
  for id_1, id_2, d in report.broken_bonds:
    p("  %s - %s: %.2f A" % (id_1, id_2, d))
  p("Clashes: %d" % len(report.clashes))
  for id_1, id_2, d in report.clashes:
    p("  %s - %s: %.2f A" % (id_1, id_2, d))
