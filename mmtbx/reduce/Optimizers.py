##################################################################################
#                Copyright 2021-2026 Richardson Lab at Duke University
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function, nested_scopes, generators, division
from __future__ import absolute_import

import argparse, contextlib, itertools, math, re, time

from boost_adaptbx import graph
from boost_adaptbx.graph import connected_component_algorithm as cca

from iotbx.map_model_manager import map_model_manager
from iotbx.data_manager import DataManager
from iotbx.pdb import common_residue_names_get_class
import mmtbx
from scitbx.array_family import flex

# To enable addition of Hydrogens
# @todo See if we can remove the shift and box once reduce_hydrogen is complete
from cctbx.maptbx.box import shift_and_box_model
from mmtbx.hydrogens import reduce_hydrogen

import boost_adaptbx.boost.python as bp
bp.import_ext("mmtbx_probe_ext")
import mmtbx_probe_ext as probeExt
from libtbx.utils import Sorry, Abort
from mmtbx.probe import Helpers
from mmtbx.reduce import Movers, InteractionGraph, HisCharges

bp.import_ext("mmtbx_reduce_ext")
from mmtbx_reduce_ext import OptimizerC, Optimizers_test

##################################################################################
# This file includes a set of functions and classes that implement placement and optimization of
# Reduce's "Movers".

##################################################################################
# Module variables that affect the way computations are handled

_DoCliqueOptimizationInC = True

##################################################################################
# Helper functions

def _VerboseCheck(verbosity, level, message):
  # Returns "" if the level is less than verbosity level, message if it is at or above
  if verbosity >= level:
    return " "*level + message
  else:
    return ""

_lastTime = None
def _ReportTiming(verbosity, message):
  """Use None message to start the timer without printing.
  """
  import time

  global _lastTime
  if message is None:
    _lastTime = time.time()
    return
  curTime = time.time()
  diff = curTime - _lastTime
  _lastTime = curTime
  return _VerboseCheck(verbosity, 2,"Time to {}: {:0.3f}".format(message,diff)+"\n")

class _InfoAccumulator(object):
  """Collects the Optimizer's information string.  += appends a part in O(1),
  str() joins once.  A plain str attribute is copied in full on every
  self._infoString += ..., which is quadratic in the number of Movers."""
  def __init__(self):
    self._parts = []
  def __iadd__(self, s):
    if s:
      self._parts.append(s)
    return self
  def __str__(self):
    return "".join(self._parts)

def AlternatesInModel(model):
  """Returns a set of altloc names of all conformers in all chains.
  :param model: pdb.hierarchy.model to search for alternates.
  :return: Set of strings.  The set is will include only the empty string if no
  chains in the model have alternates.  It has an additional entry for every altloc
  found in every chain.
  """
  ret = set(['']) # We always add this so that it will be present even when there are alternates
  for c in model.chains():
    for alt in c.conformers():
      ret.add(alt.altloc)
  return ret

def ResidueLetters(rg):
  """The altloc letters of residue group rg, sorted; [] for a residue without altlocs."""
  return sorted(set(ag.altloc.strip() for ag in rg.atom_groups()) - set(['']))

def ResidueAtomsForConformer(rg, conf, sameSite=False):
  """The atoms of residue group rg in the pass for the non-blank altloc conf: its
  blank atoms, and
  - its conf atoms if it has conf;
  - otherwise none if sameSite (another residue group of its chain with the same
    resseq and icode has conf: microheterogeneity, e.g. Pro A / Ser B,C; the
    blank atoms stay, as a shared N may be stored in either group);
  - otherwise, with one other letter, that letter's atoms only if every one of
    them has occupancy >= 0.99 (treated as present in every pass, whatever its
    label; a deposited "full" is often 0.99), else none (e.g. a glycan copy
    modelled as altloc B, or a side chain modelled in B only);
  - otherwise (several other letters) the alphabetically first letter's atoms.
  A residue never gives atoms of two letters. The same in every chain, whether or
  not the chain has conf, so a hetero residue gets the same atoms with or without
  TER (PDB) or its own label_asym_id (mmCIF)."""
  letters = ResidueLetters(rg)
  if conf in letters or not letters:
    keep = conf
  elif sameSite:
    keep = None
  elif len(letters) == 1:
    keep = letters[0]
    for ag in rg.atom_groups():
      if ag.altloc.strip() == keep and not all(a.occ >= 0.99 for a in ag.atoms()):
        keep = None
        break
  else:
    keep = letters[0]
  return [a for a in rg.atoms() if a.parent().altloc.strip() in ('', keep)]

def GetAtomsForConformer(model, conf):
  """Returns a list of atoms in the named conformer.  For a model without altlocs,
  it will return all of the atoms.  For a non-blank conf, every residue gives
  ResidueAtomsForConformer(): its blank atoms, its conf atoms, and for a residue
  without conf at most one other conformer as environment (none if a residue
  group at the same site has conf; see there). This also
  holds for the single pass run when the Optimizer is given an altID.
  :param model: pdb.hierarchy.model to search for atoms.
  :param conf: String name of the conformation to find.  Can be "", which finds the
  default conformer; if there is no empty conformation, then it will
  pick the first available conformation for each atom group.
  :return: List of atoms consistent with the specified conformer.
  """
  ret = []
  for ch in model.chains():
    confs = ch.conformers()
    if conf.strip() == '' or (len(confs) == 1 and confs[0].altloc == ''):
      ret += confs[0].atoms()
      continue
    rgs = ch.residue_groups()
    sites = set((rg.resseq, rg.icode) for rg in rgs if conf in ResidueLetters(rg))
    for rg in rgs:
      ret += ResidueAtomsForConformer(rg, conf, (rg.resseq, rg.icode) in sites)
  return ret

def _ReferencePeptideBond(prev, rg, bondTable):
  """The first restrained C-N bond between residue groups prev and rg (any
  letters), as (c, n, distanceIdeal); None if there is none."""
  for c in prev.atoms():
    if c.name.strip() != 'C':
      continue
    for n in rg.atoms():
      if n.name.strip() != 'N':
        continue
      bp = bondTable.lookup(c.i_seq, n.i_seq)
      if bp is not None:
        return (c, n, bp.distance_ideal)
  return None

def CrossLetterPeptideBonds(model, conf, bondTable):
  """The peptide bonds that the pass for the non-blank altloc conf holds between
  two different letters: C of a residue and N of the next one as
  GetAtomsForConformer() gives them (e.g. C of Ala A, its first letter, with N of
  Ser C). Restraints are made per conformer, so no proxy bonds such a pair. It
  counts only if C and N of the same residues are bonded in another letter
  combination (the reference proxy), and is accepted as a bond within 1.25 times
  the reference's ideal distance. Placement, the Optimizer's neighbor lists and
  so its Movers all use this one decision.
  :param model: pdb.hierarchy.model.
  :param conf: String name of the conformation.
  :param bondTable: bond_params_table of the model's geometry restraints.
  :return: List of (c, n, cRef, nRef, distanceIdeal, accepted) in chain order.
  """
  ret = []
  if conf.strip() == '':
    return ret
  for ch in model.chains():
    if len(ch.conformers()) == 1:
      continue
    rgs = ch.residue_groups()
    sites = set((rg.resseq, rg.icode) for rg in rgs if conf in ResidueLetters(rg))
    def held(rg, name):
      for a in ResidueAtomsForConformer(rg, conf, (rg.resseq, rg.icode) in sites):
        if a.name.strip() == name:
          return a
      return None
    # the residue groups at the previous site (several for microheterogeneity)
    lastSite, prevs, here = None, [], []
    for rg in rgs:
      if (rg.resseq, rg.icode) != lastSite:
        prevs, here, lastSite = here, [], (rg.resseq, rg.icode)
      here.append(rg)
      if not ResidueLetters(rg):
        continue
      n = held(rg, 'N')
      if n is None or n.parent().altloc.strip() == '':
        continue
      for prev in prevs:
        if not ResidueLetters(prev):
          continue
        c = held(prev, 'C')
        if c is None:
          continue
        lc = c.parent().altloc.strip()
        if lc in ('', n.parent().altloc.strip()):
          continue
        if bondTable.lookup(c.i_seq, n.i_seq) is not None:
          continue
        ref = _ReferencePeptideBond(prev, rg, bondTable)
        if ref is None:
          continue
        ret.append((c, n, ref[0], ref[1], ref[2],
          c.distance(n) <= 1.25 * ref[2]))
  return ret

def EnvironmentResidues(model, conf, alts):
  """The residue groups of model that are environment in the pass for conf: those
  with altlocs but without conf, with at least one of their letters among the
  passes alts. Their atoms in the pass (ResidueAtomsForConformer(), blank atoms
  included) take part in scoring and get Movers, but the pass writes nothing for
  them: their own passes do. Every other residue is owned by the pass."""
  if conf.strip() == '':
    return []
  passes = set(x.strip() for x in alts)
  ret = []
  for ch in model.chains():
    confs = ch.conformers()
    if len(confs) == 1 and confs[0].altloc == '':
      continue
    for rg in ch.residue_groups():
      letters = ResidueLetters(rg)
      if letters and conf not in letters and passes.intersection(letters):
        ret.append(rg)
  return ret

def _HisAcceptorsInAnyState(mover):
  """i_seq of the ring Nitrogens (NE2, ND1) of a MoverHisFlip that are acceptors in at
  least one of its coarse states."""
  coarse = mover.CoarsePositions()
  ret = set()
  for extras in coarse.extraInfos:
    for i in (0, 4):
      if extras[i].isAcceptor:
        ret.add(coarse.atoms[i].i_seq)
  return ret

def _ProtonationKey(a):
  # (model id, chain id, resseq, icode) of the atom's residue, stripped
  rg = a.parent().parent()
  return (rg.parent().parent().id.strip(), rg.parent().id.strip(), rg.resseq.strip(),
    rg.icode.strip())

def _ResNameAndID(a):
  """Make a string describing the residue and chain for the specfied atom.
  """
  chainID = a.parent().parent().parent().id
  resName = a.parent().resname.strip().upper()
  resID = str(a.parent().parent().resseq_as_int())
  altLoc = a.parent().altloc
  # Don't print the code if it is a space (blank).
  insertionCode = a.parent().parent().icode.strip()
  return "chain "+str(chainID)+" "+altLoc+resName+" "+resID+insertionCode

##################################################################################
# Helper classes

class FlipMoverState(object):
  # Holds information needed to identify a flip Mover within a model file and
  # whether or not it is flipped and has its angles adjusted.
  def __init__(self, moverType, modelId, altId, chain, resName, resIdWithICode, flipped, fixedUp):
    self.moverType = moverType  # String indicating Mover type: AmideFlip or HisFlip
    self.modelId = modelId      # Integer indicating the modelId that the entry corresponds to
    self.altId = altId          # String indicating the altId that the entry corresponds to
    self.chain = chain          # String Chain that the residue is in
    self.resName = resName      # String Name of the residue
    try:                        # String holding the integer ID of the residue and any insertion code
      self.resId = int(resIdWithICode)
      self.iCode = ''
    except Exception:
      self.resId = int(resIdWithICode[:-1])
      self.iCode = resIdWithICode[-1]
    self.flipped = flipped      # Boolean Whether the Mover is flipped in this configuration
    self.fixedUp = fixedUp      # Boolean Whether the fixup has been done on the angles

  def __str__(self):
    return "{} {} '{}' {} {} {}{} {} {}".format(self.moverType, self.modelId, self.altId,
      self.chain, self.resName, self.resId, self.iCode, self.flipped, self.fixedUp)
  def __repr__(self):
      return "Optimizers.FlipMoverState({})".format(str(self))

class FlippedMoverInfo(object):
  def __init__(self, alt, baseAtom):
    self.alt = alt              # String altId of the conformation that is flipped
    self.baseAtom = baseAtom    # Atom that is the base of the flip, which is the N for
                                # an amide flip and the NE for a histidine flip

  def __str__(self):
    return "alt '{}', base atom {}".format(self.alt.strip().upper(), int(self.baseAtom.i_seq))
  def __repr__(self):
      return "Optimizers.FlippedMoverInfo({})".format(str(self))

##################################################################################
# Optimizer, which wraps the OptimizerC class to do the optimization:

def _OperatorKey(rt):
  """Hashable key of a symmetry operator, the same for any denominators."""
  return tuple(int(round(v * 24)) for v in rt.as_double_array())

def _GridKey(xyz, size):
  return (int(math.floor(xyz[0]/size)), int(math.floor(xyz[1]/size)),
          int(math.floor(xyz[2]/size)))

def _AddToGrid(grid, key, value):
  if key not in grid:
    grid[key] = []
  grid[key].append(value)

_HIS_RING_NAMES = ('CG', 'ND1', 'CD2', 'CE1', 'NE2', 'HD1', 'HE2')

def _HisKey(mi, a):
  rg = a.parent().parent()
  return (mi, rg.parent().id.strip(), rg.resseq.strip(), rg.icode.strip())

def _HisMid(m):
  # Mover id: its atoms' i_seqs (a display label can repeat across conformers)
  return tuple(a.i_seq for a in m.CoarsePositions().atoms)

def _HisLetter(atoms):
  """Conformer letter of His ring atoms and H, and whether blank and split mix."""
  alts = [a.parent().altloc.strip() for a in atoms]
  letters = set(x for x in alts if x)
  return (sorted(letters)[0] if letters else ''), bool(letters) and ('' in alts or len(letters) > 1)

def _HisRing(ne2, inPass):
  """(conformer letter, partial, {name: atom}) for the ring and ring H of a His in a
  pass. The ring H name the conformer; a blank ring with split H is a whole one."""
  found = {}
  for a in ne2.parent().parent().atoms():
    if a in inPass and a.name.strip() in ('CG', 'ND1', 'CD2', 'CE1', 'NE2'):
      found[a.name.strip()] = a
  found['NE2'] = ne2
  for n, hn in (('NE2', 'HE2'), ('ND1', 'HD1')):
    if n in found:
      hs = [h for h in inPass[found[n]] if h.element_is_hydrogen()]
      if hs:
        found[hn] = hs[0]
  hs = [found[n] for n in ('HE2', 'HD1') if n in found]
  letter, partial = _HisLetter(hs if hs else found.values())
  return letter, partial, found

def _HisBlankRingLastPassLacksH(recs, key):
  """The letter of the last pass recording His key when that pass sees the blank
  ring without its ring H (it resets the ring and sets the ring N acceptors), else
  None. recs: the run's His records in pass order."""
  mine = [r for r in recs if r['key'] == key]
  if mine and mine[-1].get('hNotInPass'):
    return mine[-1]['passAlt']
  return None

class _HisFirstRunError(Exception):
  """An error in the first run's reconciliation-only code; error: the original."""
  def __init__(self, error):
    Exception.__init__(self, '%s: %s' % (type(error).__name__, error))
    self.error = error

def _HisOutcomeText(run):
  """The model outcome's text, from its code (HisCharges.OUTCOME_TEXT) and detail."""
  text = HisCharges.OUTCOME_TEXT[run['code']]
  return text + (' (%s)' % run['detail'] if run['detail'] else '')

@contextlib.contextmanager
def his_charges_strict(strict=True):
  """Tests: within the block an error while reconciling His charges is re-raised
  after the free result is restored (Optimizer.hisChargesStrict); the previous
  value is always restored."""
  old = Optimizer.hisChargesStrict
  Optimizer.hisChargesStrict = strict
  try:
    yield
  finally:
    Optimizer.hisChargesStrict = old

class _FixedHisOptimizer(object):
  """The OptimizerC of a His that was a one-state Mover: reports its original
  state index, everything else from the OptimizerC."""
  def __init__(self, optC, pin):
    self._optC = optC
    self._pin = pin
  def GetCoarseLocation(self, m):
    return self._pin
  def __getattr__(self, name):
    return getattr(self._optC, name)
  def __eq__(self, other):
    return self._optC == getattr(other, '_optC', other)
  def __ne__(self, other):
    return not self == other
  def __hash__(self):
    return hash(self._optC)

def _SameAtomInfo(a, b):
  return (a.vdwRadius == b.vdwRadius and a.isAcceptor == b.isAcceptor and
          a.isDonor == b.isDonor and a.isDummyHydrogen == b.isDummyHydrogen and
          a.isIon == b.isIon and a.charge == b.charge)

def _SymmetryCopyBonds(originals, copies, neighborLists, crossBonds, unitCell):
  """Bonded neighbors of symmetry copies, and the bonds from originals to copies.
  A copy under R is bonded to its original's neighbors under R and to a partner across
  symmetry under R times the bond's operator; an original to its partner across symmetry
  under the bond's operator. A partner under an operator whose copy was left out as
  coincident (special position) is the atom kept at that site, original or copy.
  :param originals: {i_seq: atom}, the original atoms of the pass.
  :param copies: list of (copy, original, rt_mx), copies at their positions.
  :param neighborLists: bonded neighbors of the originals (within the asymmetric unit).
  :param crossBonds: list of (i_seq, j_seq, rt_mx_ji): j moved by rt_mx_ji is bonded to i.
  :param unitCell: unit cell of the operators.
  :return: ({copy: list of bonded atoms}, list of (atom, bonded atom) to add to the
  originals' and copies' lists: the bonds to copies of originals, and the reverse of the
  bonds of copies to atoms resolved by site).
  """
  from cctbx import sgtbx
  from cctbx.crystal import super_cell
  identity = _OperatorKey(sgtbx.rt_mx())
  tolerance = 2 * super_cell.min_distance_sym_equiv_coincident
  images = {}                     # original i_seq: {operator: copy}
  grid = {}                       # sites of originals and copies, for coincident atoms
  for a in originals.values():
    _AddToGrid(grid, _GridKey(a.xyz, 1.0), a)
  for c, o, rt in copies:
    if o.i_seq not in images:
      images[o.i_seq] = {}
    images[o.i_seq][_OperatorKey(rt)] = c
    _AddToGrid(grid, _GridKey(c.xyz, 1.0), c)

  def image(j, rt):
    # Atom j moved by rt; None outside the copies.
    key = _OperatorKey(rt)
    if key == identity:
      return originals.get(j)
    c = images.get(j, {}).get(key)
    if c is not None or j not in originals:
      return c
    site = Helpers.rvec3(unitCell.orthogonalize(rt * unitCell.fractionalize(originals[j].xyz)))
    element = originals[j].element.strip().upper()
    k = _GridKey(site, 1.0)
    for d in itertools.product((-1, 0, 1), repeat=3):
      for a in grid.get((k[0]+d[0], k[1]+d[1], k[2]+d[2]), []):
        if (a.element.strip().upper() == element and
            (Helpers.rvec3(a.xyz) - site).length() < tolerance):
          return a
    return None

  # Bonds across symmetry by original: (partner i_seq, operator of the partner).
  byOriginal = {}
  for i, j, rt in crossBonds:
    if i in originals and j in originals:
      _AddToGrid(byOriginal, i, (j, rt))
      _AddToGrid(byOriginal, j, (i, rt.inverse()))
  copyNeighbors = {}
  for c, o, rt in copies:
    partners = [image(n.i_seq, rt) for n in neighborLists[o]]
    for j, bondRt in byOriginal.get(o.i_seq, []):
      partners.append(image(j, rt.multiply(bondRt)))
    copyNeighbors[c] = []
    for p in partners:
      if p is not None and p != c and p not in copyNeighbors[c]:
        copyNeighbors[c].append(p)
  originalBonds = []
  for i, partners in byOriginal.items():
    for j, bondRt in partners:
      p = image(j, bondRt)
      if p is not None and p != originals[i]:
        originalBonds.append((originals[i], p))
  for c, o, rt in copies:
    for p in copyNeighbors[c]:
      if p in copyNeighbors:
        reverse = c in copyNeighbors[p]
      else:
        reverse = (p, c) in originalBonds
      if not reverse:
        originalBonds.append((p, c))
  return copyNeighbors, originalBonds

# vertices visited by the cut searches of one count: ~1 s of Python
_CLIQUE_WORK_MAX_STEPS = 5000000
# _CliqueSearchWork() gave up: its cut searches alone passed _CLIQUE_WORK_MAX_STEPS
_CLIQUE_WORK_ABANDONED = -1
# OptimizeCliqueCoarseVertexCut() at verbosity >= 1, when a (sub)Clique falls back.
# These are the only C++ lines below verbosity 3: OptimizerC always gets >= 1.
_CLIQUE_SINGLETONS_LINE = re.compile(r'(Out of time|No time remaining), treating \(sub\)clique of size \d+ as singletons')
# start of the warning about a Clique whose search was limited (reduce2 reads it)
_CLIQUE_WARNING = 'Warning: search limit: '

def _LimitedCliques(warnings):
  """The Cliques (size, model, alternate, first Mover) named by the search-limit warnings in warnings."""
  return set(l[len(_CLIQUE_WARNING):].split('): ')[0] for l in warnings.split('\n')
    if l.startswith(_CLIQUE_WARNING))

def _CliqueSearchWork(n, edges, numStates, limit):
  """Work count of OptimizerC.OptimizeCliqueCoarse() on a Clique, from its graph
  alone: state combinations scored plus subsets tried to find cuts, with the cuts as
  findVertexCut() picks them (fewest Movers, then most even components, first in
  nChooseM order). Each state combination scores all Movers of its (sub)Clique and
  the fixed ones around it. Units, added up for each call of the recursion:
  - no cut (or at most 2 Movers): one per combination of the Movers' states (brute
    force), plus the 2**N - 1 subsets the cut search tried (none at N <= 2);
  - a cut of K Movers: the subsets of 0..K Movers tried, plus for each combination
    of the cut's states one (the rescoring) and the work of each component, which
    the C++ searches again, cut search included, for every such combination.
  :param n: number of Movers, vertices 0..n-1 in the order given to the C++.
  :param edges: list of (i, j) vertex pairs.
  :param numStates: coarse states of each Mover.
  :param limit: stop counting above it.
  :return: the count, None when above limit, or _CLIQUE_WORK_ABANDONED when the cut
  searches would visit more than _CLIQUE_WORK_MAX_STEPS vertices (the C++ would try
  those subsets too, and again for each state of an enclosing cut).
  """
  adj = dict((i, set()) for i in range(n))
  for i, j in edges:
    adj[i].add(j)
    adj[j].add(i)
  steps = [0]   # vertices visited by the cut searches of this count, each subproblem once

  class Abandoned(Exception):
    pass

  def components(keep):
    left = set(keep)
    out = []
    while left:
      todo = [min(left)]
      left.discard(todo[0])
      comp = []
      while todo:
        v = todo.pop()
        comp.append(v)
        for w in adj[v] & left:
          left.discard(w)
          todo.append(w)
      out.append(comp)
    return out

  memo = {}
  def work(nodes):
    # nodes: tuple in graph order, as subsetGraph() keeps them
    if nodes in memo:
      return memo[nodes]
    N = len(nodes)
    w = cut(nodes) if N > 2 else math.prod(numStates[v] for v in nodes)
    if w is None or w > limit:
      return None
    memo[nodes] = w
    return w

  def cut(nodes):
    N = len(nodes)
    best = None
    for k in range(N):
      for removed in itertools.combinations(range(N), k):
        steps[0] += N
        if steps[0] > _CLIQUE_WORK_MAX_STEPS:
          raise Abandoned()
        keep = [nodes[i] for i in range(N) if i not in removed]
        if len(keep) < 2:
          continue
        comps = components(keep)
        if len(comps) > 1:
          sizes = [len(c) for c in comps]
          key = (max(sizes) - min(sizes), max(sizes))
          if best is None or key < best[0]:
            best = (key, [nodes[i] for i in removed], comps)
      if best is not None:
        break
    if best is None:
      w = math.prod(numStates[v] for v in nodes)
      w += 2**N - 1   # every subset of 0..N-1 Movers tried
    else:
      cutStates = math.prod(numStates[v] for v in best[1])
      sub = 0
      for c in best[2]:
        s = work(tuple(sorted(c)))
        if s is None:
          return None
        sub += s
      w = cutStates * (sub + 1) + sum(math.comb(N, k) for k in range(len(best[1]) + 1))
    return w

  try:
    return work(tuple(range(n)))
  except Abandoned:
    return _CLIQUE_WORK_ABANDONED

class Optimizer(object):

  def __init__(self, probePhil, addFlipMovers, model, modelIndex = None, altID = None,
                bondedNeighborDepth = 4,
                useNeutronDistances = False,
                minOccupancy = 0.02,
                preferenceMagnitude = 1.0,
                nonFlipPreference = 0.5,
                skipBondFixup = False,
                flipStates = '',
                verbosity = 1,
                cliqueOutlineFileName = None,
                fillAtomDump = True,
                stopOptimizingAt = time.time() + 1.0e10,
                symmetryCopies = None,
                fixedHisProtonation = None,
                chooseHisProtonation = False,
                reconcileHisCharges = False,
                cliqueSeconds = 120.0,
                maxCliqueWork = 1.0e6,
                reconcileHisStates = False
              ):
    """Constructor.  This is the wrapper class for the C++ OptimizerC and
    it implements the machinery that finds and optimizes Movers.
    :param probePhil: Phil parameters to Probe to be passed on when subroutines are called.
    :param addFlipMovers: Do we add flip Movers along with other types?
    :param model: iotbx model (a group of hierarchy models).  Can be obtained using
    iotbx.map_model_manager.map_model_manager.model().  The model must have Hydrogens,
    which can be added using mmtbx.hydrogens.reduce_hydrogen.place_hydrogens().get_model().
    It must have a valid unit cell, which can be helped by calling
    cctbx.maptbx.box.shift_and_box_model().  It must have had PDB interpretation run on it,
    which can be done using model.process(make_restraints=True) with PDB
    interpretation parameters and hydrogen placement matching the value of the
    useNeutronDistances parameter described below.
    :param modelIndex: Identifies which index from the hierarchy is to be selected.
    If this value is None, optimization will be run sequentially on every model in the
    hierarchy. This is 1-based, the first model is 1.
    :param altID: The conformer alternate location specifier to use.  The value "" will
    cause it to run on the first conformer found in each model.  If this is set to None,
    optimization will be run sequentially for every conformer in the model, starting with
    the last and ending with the first.  This will leave the initial conformer's values as the
    final location for atoms that are not inside a conformer or are in the first conformer;
    the blank atoms of a residue with altlocs keep the values of the last pass that has one
    of its altlocs (the other passes hold it as environment, EnvironmentResidues()).
    :param bondedNeighborDepth: How many hops to ignore bonding when doing Probe calculations.
    A depth of 3 will ignore my bonded neighbors and their bonded
    neighbors and their bonded neighbors.
    :param useNeutronDistances: False will use X-ray/electron cloud distances.  If set to
    True, it will use neutron (nuclear) distances.  This must be set consistently with the
    values used to generate the hydrogens and to run PDB interpretation.
    :param minOccupancy: Minimum occupancy for an atom to be considered in the Probe score.
    :param preferenceMagnitude: Multiplier for the preference energies expressed
    by some Movers for particular orientations.
    :param nonFlipPreference: Preference for not flipping Movers that are flips.  This keeps
    them from being flipped unless their score is significantly better than in the original position.
    :param skipBondFixup: Should we do fixup on Movers or just leave them flipped?  We always do Hydrogen
    removal and fixup for Histidines that are not placed as Movers, but for flips that
    are Movers we don't adjust the bond angles.  This is overridden for Flips that are
    specified in flipStates.
    :param flipStates: String with comma-separated entries. Each entry has the form
    (without the single quotes) '1 . A HIS 11H Flipped AnglesAdjusted'. These are space-separated values.
    The first word is the model number, starting with 1. The second is the lower-case alternate, or
    '.' for all alternates (also use this when there are no alternates in the file).
    The third is the chain ID. The fourth is the residue name. The fifth is the residue id,
    which may include an insertion code as its last character. The sixth is either Flipped or Unflipped.
    If it is Flipped, then another word is added -- AnglesAdjusted or AnglesNotAdjusted,
    specifying whether to do the three-point dock to adjust the bond angles after the flip.
    An example with several entries is (again, no quotes are included:
    '1 a A HIS 11H Unflipped,1 b A ASN 15 Flipped AnglesNotAdjusted,1 . B GLN 27 Flipped AnglesAdjusted'. Any
    Flip Movers that would be placed at the specified location are instead locked in the
    specified configuration.
    :param verbosity: Default value of 1 reports standard information.
    Value of 2 reports timing information.
    Setting it to 0 removes all informational messages.
    Setting it above 1 provides additional debugging information.
    :param cliqueOutlineFileName: Name of file to write Kinemage with outlines of Movers to.
    This file holds spheres showing all possible locations for each atom in each Mover in different
    colors as one master. It shows the outlines expanded by the probe radius, with a single color
    for each clique, as another master. These are useful for determining why the cliques are as
    they are.
    :param fillAtomDump: If true, fill in the atomDump string with the atom information.
    This can take a long time to do, so the caller may want to turn it off if they don't need it.
    :param stopOptimizingAt: Start treating cliques as singletons at this time relative to time.time().
    :param cliqueSeconds: Search time allowed to one multi-Mover Clique (a backstop: brute
    force and the cut search are not checked, so a search can run past it). Out of time,
    the part still to search is treated as singletons. None: only stopOptimizingAt.
    :param maxCliqueWork: A Clique whose search work count (_CliqueSearchWork(): state
    combinations scored plus subsets tried to find cuts) is above this, or too large to
    count, is treated as singletons without searching. The largest normal Cliques count
    about 1.5e4. None: no limit.
    Limited Cliques get a warning that starts with _CLIQUE_WARNING.
    :param symmetryCopies: None, or symmetry copies of atoms near the model as a fixed
    environment: a group_args with crystal_symmetry and models, {model index (0-based):
    list of (copy atom, original atom, rt_mx)}, copy = rt_mx * original in fractional
    coordinates. The copies are scored, never moved; Mover placement sees them (ions,
    candidate acceptors). After the optimization of each alternate they follow their
    originals, and the Cliques near copies that changed are optimized again (at most 3
    passes). Sites still changing after the last pass keep the combination of their
    states that scores best with matching copies.
    :param fixedHisProtonation: None, or the His whose protonation the user set, as
    (model id, chain id, resseq, icode), stripped (reduce_hydrogen.protonation_key()):
    they get no His flip Mover and no ion lock-down, with or without addFlipMovers.
    :param chooseHisProtonation: With addFlipMovers False, each His still chooses HD1,
    HE2 or both, among its unflipped states only, and the Asn/Gln/His flips that would
    score better are listed (getSuggestedFlips()). The ion lock-down then sets a His
    only in its unflipped orientation: a His whose ion binds only the flipped ring is
    listed with the ion. No effect with addFlipMovers. Without this parameter and
    without addFlipMovers the ion lock-down still flips a His to bind its ion.
    :param reconcileHisCharges: His whose conformers end with different charges are
    reconciled to one charge where the scores allow it (mmtbx.reduce.HisCharges):
    the passes run again with those His states fixed, and the result is kept only if
    it passes the rerun checks; otherwise the first result is restored.
    :param reconcileHisStates: with reconcileHisCharges, the target is one ring-H
    state (HD1, HE2 or both) for all conformers instead of one charge, without a cost
    limit (a cost above HisCharges.BUDGET is a warning). Each conformer keeps its ring
    orientation from the first run. The exceptions (metal, ion, user state, partial
    or unscorable ring) and the rerun checks are the same.
    """

    ################################################################################
    # Store the parameters that will be accessed by other methods
    self._probePhil = probePhil
    self._bondedNeighborDepth = bondedNeighborDepth
    self._probeRadius = probePhil.probe_radius
    self._useNeutronDistances = useNeutronDistances
    self._probeDensity = probePhil.density
    self._minOccupancy = minOccupancy
    self._preferenceMagnitude = preferenceMagnitude
    self._nonFlipPreference = nonFlipPreference
    self._skipBondFixup = skipBondFixup
    self._flipStates = flipStates
    self._verbosity = verbosity
    self._cliqueOutlineFileName = cliqueOutlineFileName
    self._symmetryCopies = symmetryCopies
    self._fixedHis = set(fixedHisProtonation or ())
    self._chooseHis = chooseHisProtonation and not addFlipMovers
    self._reconcileHisCharges = reconcileHisCharges
    self._cliqueSeconds = cliqueSeconds
    self._maxCliqueWork = maxCliqueWork
    self._cliqueModel = 0
    self._cliqueAlt = ''
    self._reconcileHisStates = reconcileHisStates
    self._hisBudget = None if reconcileHisStates else HisCharges.BUDGET

    ################################################################################
    # Initialize internal variables.
    self._infoString = _InfoAccumulator()
    self._warningString = ""
    self._atomDump = ""
    self._waterOccCutoff = 0.66 # @todo Make this a parameter, -WaterOCCcutoff param in reduce
    self._waterBCutoff = 40.0   # @todo Make this a parameter, -WaterBcutoff param in reduce
    self._numCalculated = 0
    self._numCached = 0
    self._amideFlipsPerformed = []  # List of amide flips performed, whether by Mover or locking in place.
    self._hisFlipsPerformed = []    # List of histidine flips performed, whether by Mover or locking in place.
    self._suggestedFlips = {}       # flips scored but not made (chooseHisProtonation), by residue
    self._hydrogensToDelete = []    # Each as decided by the last alternate pass that contained it.
    self._environmentAtoms = []     # Symmetry copies and their Phantom Hydrogens in this pass.
    self._passCopies = []           # (copy, original, rt_mx) for the originals in this pass
    self._absentCopies = set()      # copies of deleted originals, not in the spatial query
    self._copySummary = ""
    self._hisNoMover = {}           # (alt, NE2 i_seq): why a His has no Mover
    self._envISeqs = set()          # i_seq of this pass's environment residues (EnvironmentResidues())
    self._envMovers = set()         # this pass's Movers on them
    self._hisCharges = []           # reconcileHisCharges: per model, getHisCharges()

    ################################################################################
    # Get the Cartesian positions of all of the atoms in the entire model and find
    # the bond proxies for all of them.  The proxies will not attempt to bond atoms
    # that are in different model indices.
    _ReportTiming(self._verbosity, None) # Reset timer
    carts = flex.vec3_double()
    for a in model.get_atoms():
      carts.append(a.xyz)
    self._infoString += _ReportTiming(self._verbosity, "get coordinates")
    bondProxies, asuProxies = model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)
    self._infoString += _ReportTiming(self._verbosity, "compute bond proxies")

    ################################################################################
    # Get the bonded neighbor lists for all of the atoms in the model, so we don't get
    # failures when we look up an atom from another in Helpers.getExtraAtomInfo().
    # Only for that: a blank atom here has the neighbors of every conformer (an H copied
    # into A and B counts twice), so each alternate pass below builds its own lists.
    bondedNeighborLists = Helpers.getBondedNeighborLists(model.get_atoms(), bondProxies, asuProxies)
    self._infoString += _ReportTiming(self._verbosity, "compute bonded neighbor lists")

    ################################################################################
    # Get the probeExt.ExtraAtomInfo needed to determine which atoms are potential acceptors.
    # This is done for all atoms in the model.
    # This must operate on the entire model, not just the current model index.
    ret = Helpers.getExtraAtomInfo(
      model = model, bondedNeighborLists = bondedNeighborLists,
      useNeutronDistances=self._useNeutronDistances, probePhil=self._probePhil)
    self._extraAtomInfo = ret.extraAtomInfo
    self._infoString += ret.warnings
    self._infoString += _ReportTiming(self._verbosity, "get extra atom info")
    self._warningString += ret.warnings

    ################################################################################
    # Run optimization for every desired conformer and every desired model, calling
    # placement and then a derived-class single optimization routine for each.  When
    # the modelIndex or altID is None, that means to run over all available cases.
    # For alternates, if there is a non-empty (not "" or " ") case, then we run backwards
    # from the last to the first but do not run for the empty case; if there are only
    # empty cases, then we run just once.  We run the models in order, all of them when
    # None is specified and the specified one if it is specified.

    try:
      # None, not False: callers outside reduce2 keep the scattering-table rule
      model.setup_riding_h_manager(
        use_neutron_distances=True if self._useNeutronDistances else None)
    except Exception:
      # If an optimized has already been run on the model, our riding H manager will
      # have already been set up.  When one has already been set up, this causes
      # an exception and it doesn't make a new one.
      pass
    riding_h_manager = model.get_riding_h_manager()
    h_parameterization = riding_h_manager.h_parameterization

    # Find the single-hydrogen rotators.
    # Dorothee provided this faster approach that uses the riding_h_manager.
    rotatableHydrogens = flex.size_t()
    for p in h_parameterization:
      if p is not None:
        if p.htype  == 'alg1b':
          rotatableHydrogens.append(p.ih)

    self._infoString += _ReportTiming(self._verbosity, "select rotatable hydrogens")

    # Copies take i_seq above the model's and its water phantoms' (the phantoms of a
    # pass all start at getMaxISeq() + 1), one range per model; their phantoms follow.
    if symmetryCopies is not None:
      nextISeq = Helpers.getMaxISeq(model) + 100
      for mi in sorted(symmetryCopies.models):
        for c, o, rt in symmetryCopies.models[mi]:
          probeExt.set_atom_i_seq(c, nextISeq)
          nextISeq += 1
      self._copyPhantomISeq = nextISeq
      # Bonds across symmetry (S-S, links, metal coordination) with their operator:
      # j at rt_mx_ji * j is bonded to i.
      self._crossBonds = []
      if asuProxies.size() > 0:
        # the pair proxies of get_all_bond_proxies() above, cached
        asuMappings = model.get_restraints_manager().geometry.pair_proxies(
          ).bond_proxies.asu_mappings()
        for bp in asuProxies:
          rt = asuMappings.get_rt_mx_ji(pair = bp)
          if not rt.is_unit_mx():
            self._crossBonds.append((bp.i_seq, bp.j_seq, rt))

    startModelIndex = 0
    stopModelIndex = len(model.get_hierarchy().models())
    if modelIndex is not None:
      # The command-line parameter matches the name of the model in the model file, which
      # starts with 1. The internal indexing starts with 0. So we subtract one.
      startModelIndex = (modelIndex - 1)
      stopModelIndex = startModelIndex + 1
    for mi in range(startModelIndex, stopModelIndex):
      # Get the specified model from the hierarchy.
      myModel = model.get_hierarchy().models()[mi]

      ################################################################################
      # Store the states (position and extra atom info) of all of the atoms in this model
      # so that we can restore it for atoms in a given alternate configuration before optimizing
      # each new alternate. Looked up by i_seq.
      initialAtomPositions = {}
      initialExtraAtomInfos = {}
      for a in myModel.atoms():
        initialAtomPositions[a.i_seq] = a.xyz;
        initialExtraAtomInfos[a.i_seq] = probeExt.ExtraAtomInfo(self._extraAtomInfo.getMappingFor(a))

      # Get the list of alternate conformation names present in all chains for this model.
      # If there is more than one result, remove the empty results and then sort them
      # in reverse order so we finalize all non-alternate ones to match the first.
      alts = AlternatesInModel(myModel)
      if len(alts) > 1:
        alts.discard("")
        alts.discard(" ")
      alts = sorted(list(alts), reverse=True)
      self._infoString += _ReportTiming(self._verbosity, "compute alternates")

      # If there is a specified alternate, use it.
      if altID is not None:
        alts = [altID]

      # The passes over the alternates, as a method so that a rerun is a call.
      ctx = dict(model=model, alts=alts, addFlipMovers=addFlipMovers, asuProxies=asuProxies,
        bondProxies=bondProxies, h_parameterization=h_parameterization,
        initialAtomPositions=initialAtomPositions, initialExtraAtomInfos=initialExtraAtomInfos,
        rotatableHydrogens=rotatableHydrogens, stopOptimizingAt=stopOptimizingAt,
        symmetryCopies=symmetryCopies, useNeutronDistances=useNeutronDistances)
      if self._reconcileHisCharges:
        # reruns do not take later models' search time
        stopOptimizingAt += self._ReconcileHisCharges(mi, myModel, ctx)
      else:
        self._OptimizeModel(mi, myModel, ctx)
      if symmetryCopies is not None:
        self._SummarizeCopies(mi)

      #################################################################################
      # Dump information about all of the atoms in the model into a string.
      if fillAtomDump:
        self._atomDump = Helpers.writeAtomInfoToString(myModel.atoms(), self._extraAtomInfo)
        self._infoString += _ReportTiming(self._verbosity, "dump atom info to string")

      #################################################################################
      # Report the fraction of atoms that were calculated and the fraction that were cached.
      if self._numCalculated > 0:
        self._infoString += _VerboseCheck(self._verbosity, 1,
            'Calculated : cached atom scores: {} : {}; fraction calculated {:.2f}\n'.format(
            self._numCalculated, self._numCached,
            self._numCalculated/(self._numCalculated + self._numCached)))

    ################################################################################
    # Release what only the search used: the OptimizerC with its dot cache, the Movers
    # (and the dictionaries and lists keyed by them), the interaction graph, the
    # excluded atoms, the dot scorer and the spatial query. The getters return the
    # strings and lists filled in above. Kept, they stay alive while the caller
    # finalizes and writes the model (4wsc: ~400 MB).
    self._optimizerOf = {}
    self._movers = []
    self._envMovers = set()
    self._envISeqs = set()
    self._moverInfo = {}
    self._initialScore = {}
    self._finalScore = {}
    self._interactionGraph = None
    self._atomMoverLists = None
    self._excludeDict = {}
    self._dotScorer = None
    self._spatialQuery = None

  def _OptimizeModel(self, mi, myModel, ctx, his=None):
    """The passes over the alternates of model index mi, from the state before the
    first pass. ctx: the constructor's values the passes use. his: None, or the
    record of a His charge run (_ReconcileHisCharges): its His states are fixed
    during the search and each pass's His states and scores are recorded."""
    self._cliqueModel = mi    # named in the search-limit warnings
    model = ctx['model']
    alts = ctx['alts']
    addFlipMovers = ctx['addFlipMovers']
    asuProxies = ctx['asuProxies']
    bondProxies = ctx['bondProxies']
    h_parameterization = ctx['h_parameterization']
    initialAtomPositions = ctx['initialAtomPositions']
    initialExtraAtomInfos = ctx['initialExtraAtomInfos']
    rotatableHydrogens = ctx['rotatableHydrogens']
    stopOptimizingAt = ctx['stopOptimizingAt']
    symmetryCopies = ctx['symmetryCopies']
    useNeutronDistances = ctx['useNeutronDistances']

    # Clear the Movers list for each model.  _PlaceMovers() makes a new one for each alternate.
    self._movers = []
    self._copyPassNotes = []
    # i_seq: (atom, [(alt, deleted)]) for the H of each pass; _deleteMes is reset per pass.
    passDeletions = {}
    firstAlt = True
    for alt in alts:
      # Get the atoms from the specified conformer in the model (the empty string is the name
      # of the first conformation in the model; if there is no empty conformation, then it will
      # pick the first available conformation for each atom group.
      passAtoms = GetAtomsForConformer(myModel, alt)

      # Residues without this alternate (EnvironmentResidues()) are its environment:
      # their state as their own passes left it is kept here and put back after the
      # pass, which records nothing for them (deletions, flips, His records, suggested
      # flips, report lines).
      envResidues = EnvironmentResidues(myModel, alt, alts)
      envSaved = [(a, a.xyz, probeExt.ExtraAtomInfo(self._extraAtomInfo.getMappingFor(a)))
        for rg in envResidues for a in rg.atoms()]
      self._envISeqs = set(a.i_seq for a, xyz, e in envSaved)

      # If we are doing the second or later alternate, place all atoms that are in a compatible alternate
      # back into their initial configuration so we start from the same state we would have if this were
      # the only alternate being tested.  This will ensure that we get compatible outputs when run either
      # way.  The environment's atoms in this pass start from their initial state too.
      if firstAlt:
        firstAlt = False
      else:
        envInPass = set(a.i_seq for a in passAtoms if a.i_seq in self._envISeqs)
        for a in myModel.atoms():
          if a.parent().altloc in ['', ' ', alt] or a.i_seq in envInPass:
            a.xyz = initialAtomPositions[a.i_seq]
            self._extraAtomInfo.setMappingFor(a, initialExtraAtomInfos[a.i_seq])
      numAmideFlips = len(self._amideFlipsPerformed)
      numHisFlips = len(self._hisFlipsPerformed)

      self._cliqueAlt = alt     # named in the search-limit warnings
      # Tell about the run we are currently doing.
      self._infoString += _VerboseCheck(self._verbosity, 1,"Running Reduce optimization on model index "+str(mi)+
        ", alternate '"+alt+"'\n")
      self._infoString += _VerboseCheck(self._verbosity, 1,"  bondedNeighborDepth = "+str(self._bondedNeighborDepth)+"\n")
      self._infoString += _VerboseCheck(self._verbosity, 1,"  probeRadius = "+str(self._probeRadius)+"\n")
      self._infoString += _VerboseCheck(self._verbosity, 1,"  useNeutronDistances = "+str(self._useNeutronDistances)+"\n")
      self._infoString += _VerboseCheck(self._verbosity, 1,"  probeDensity = "+str(self._probeDensity)+"\n")
      self._infoString += _VerboseCheck(self._verbosity, 1,"  minOccupancy = "+str(self._minOccupancy)+"\n")
      self._infoString += _VerboseCheck(self._verbosity, 1,"  preferenceMagnitude = "+str(self._preferenceMagnitude)+"\n")

      self._atoms = list(passAtoms)

      # The bonded neighbors within this conformer, for everything in this pass (the
      # shared lists above mix conformers). A new dictionary each pass, so the Phantom
      # Hydrogens added below stay in it.
      passNeighborLists = Helpers.getBondedNeighborLists(self._atoms, bondProxies, asuProxies)
      # A peptide bond the pass holds across two letters has no proxy; bond it here,
      # so that its N gets no single-hydrogen rotator.
      crossBonds = CrossLetterPeptideBonds(myModel, alt,
        model.get_restraints_manager().geometry.bond_params_table)
      if crossBonds:
        passByISeq = dict((a.i_seq, a) for a in self._atoms)
      for c, n, cRef, nRef, ideal, accepted in crossBonds:
        if accepted:
          c, n = passByISeq[c.i_seq], passByISeq[n.i_seq]
          passNeighborLists[c].append(n)
          passNeighborLists[n].append(c)
        self._infoString += _VerboseCheck(self._verbosity, 1,
          "  Peptide bond across letters: %s - %s %.2f A, %s\n" % (
          c.id_str().replace('pdb=', ''), n.id_str().replace('pdb=', ''),
          c.distance(n), "bonded" if accepted else
          "not bonded (over 1.25 x %.2f A)" % ideal))

      ################################################################################
      # Reset the timer
      _ReportTiming(self._verbosity, None)

      ################################################################################
      # Construct the spatial-query information needed to quickly determine which atoms are nearby
      self._spatialQuery = Helpers.createSpatialQuery(self._atoms, self._probePhil)
      self._infoString += _ReportTiming(self._verbosity, "construct spatial query")

      ################################################################################
      # Find the radius of the largest atom we'll have to deal with.
      maxVDWRad = 1
      for a in self._atoms:
        maxVDWRad = max(maxVDWRad, self._extraAtomInfo.getMappingFor(a).vdwRadius)
      self._maximumVDWRadius = maxVDWRad

      ################################################################################
      # Make a set of atoms that are to be deleted based on analysis.  It is initially
      # empty, keeping all of the added Hydrogens in the model.
      self._deleteMes = set()

      ################################################################################
      # Placement of water phantom Hydrogens, including adding them to our 'atoms' list
      # and the spatial query but not adding them to the hierarchy.  This must be done after
      # the bond proxies are constructed and the Movers have been placed so that these
      # fake atoms do not confuse placement.
      # These atoms will not be part of any atom group, and we'll fill in their extra atom
      # info directly.

      # @todo Look up the radius of a water Hydrogen.  This may require constructing a model with
      # a single water in it and asking about the hydrogen radius.  This could also become a
      # Phil parameter.  Also look up the OH bond distance rather than hard-coding it here.
      phantomHydrogenRadius = 1.05
      placedHydrogenDistance = 0.84
      if useNeutronDistances:
        phantomHydrogenRadius = 1.0
        placedHydrogenDistance = 0.98

      # Find every Oxygen that is part of a water and get the Phantom Hydrogens for it
      # unless it has out-of-bounds parameter values.  If the water Oxygen has out-of-bounds
      # parameters, then we remove it from the atoms to be considered (but not from the
      # model) by removing them from the atom list and from the spatial query structure.
      phantoms = [] # List of tuples, the Phantom and its parent Oxygen
      watersToDelete = []
      self._watersWithH = set()   # i_seq of the water Oxygens with H or D in this pass
      maxISeq = Helpers.getMaxISeq(model)
      for a in self._atoms:
        if a.element == 'O' and common_residue_names_get_class(name=a.parent().resname) == "common_water":
          if a.occ >= self._waterOccCutoff and a.b < self._waterBCutoff:
            # We're an acceptor and not a donor.
            ei = self._extraAtomInfo.getMappingFor(a)
            ei.isDonor = False
            ei.isAcceptor = True
            self._extraAtomInfo.setMappingFor(a, ei)

            # Any H or D bonded in this conformer, one as well as two, stops the
            # phantoms: the model's H as given, so the second donor of a one-H
            # water is not represented. Other bonded atoms (a metal linked by
            # the caller's restraints) leave the Oxygen without H.
            numHD = 0
            for n in passNeighborLists[a]:
              if n.element_is_hydrogen():
                numHD += 1
            if numHD > 0:
              self._watersWithH.add(a.i_seq)
              self._infoString += _VerboseCheck(self._verbosity, 3,
                "Not adding phantom Hydrogens on {}: it carries {} H/D\n".format(
                _ResNameAndID(a), numHD))
              continue

            newPhantoms = Helpers.getPhantomHydrogensFor(maxISeq, a, self._spatialQuery, self._extraAtomInfo, self._minOccupancy,
                            False, phantomHydrogenRadius, placedHydrogenDistance)
            if len(newPhantoms) > 0:
              resNameAndID = _ResNameAndID(a)
              self._infoString += _VerboseCheck(self._verbosity, 3,"Added {} phantom Hydrogens on {}\n".format(len(newPhantoms), resNameAndID))
              for p in newPhantoms:
                self._infoString += _VerboseCheck(self._verbosity, 5,"Added phantom Hydrogen at "+str(p.xyz)+"\n")
                phantoms.append( (p,a) )

          else:
            # Occupancy or B factor are out of bounds, so remove this atom from consideration.
            self._infoString += _VerboseCheck(self._verbosity, 3,"Ignoring "+
              a.name.strip()+" "+a.parent().resname.strip()+" "+str(a.parent().parent().resseq_as_int())+
              " "+str(a.parent().parent().parent().id)+
              " with occupancy "+str(a.occ)+" and B factor "+str(a.b)+"\n")
            watersToDelete.append(a)

      if len(self._watersWithH) > 0:
        self._infoString += _VerboseCheck(self._verbosity, 1,"No phantom Hydrogens on "+
          str(len(self._watersWithH))+" waters that carry H or D\n")
      if len(watersToDelete) > 0:
        self._infoString += _VerboseCheck(self._verbosity, 1,"Ignored "+str(len(watersToDelete))+" waters due to occupancy or B factor\n")
        # One pass: list.remove() per water scans the whole list each time.
        # Atom hashing and equality use memory_id(), as list.remove() does.
        ignored = set(watersToDelete)
        self._atoms[:] = [a for a in self._atoms if a not in ignored]
        for a in watersToDelete:
          self._spatialQuery.remove(a)

      if len(phantoms) > 0:

        # Add these atoms to the list of atoms we deal with.
        # Add these atoms to the spatial-query structure.
        # Insert ExtraAtomInfo for each of these atoms, marking each as a dummy and as a donor.
        # Add to the bondedNeighborList with their parent Oxygen as bonded one way so that dots on
        # a Phantom Hydrogen within its Oxygen will be excluded.  Do not mark the Oxygen as being
        # bonded to the Phantom Hydrogen to avoid having it mask collisions between the Oxygen and
        # other atoms.
        origCount = len(self._atoms)
        for p in phantoms:
          self._atoms.append(p[0])
          self._spatialQuery.add(p[0])
          eai = probeExt.ExtraAtomInfo(phantomHydrogenRadius, False, True, True)
          self._extraAtomInfo.setMappingFor(p[0], eai)
          passNeighborLists[p[0]] = [p[1]]

        self._infoString += _VerboseCheck(self._verbosity, 1,"Added "+str(len(phantoms))+" phantom Hydrogens on waters")
        self._infoString += _VerboseCheck(self._verbosity, 1," (Old total "+str(origCount)+", new total "+str(len(self._atoms))+")\n")
      self._infoString += _ReportTiming(self._verbosity, "place water phantom Hydrogens")

      ################################################################################
      # Fix up the donor status for all of the atoms now that we've added the final explicit
      # Phantom Hydrogens.
      Helpers.fixupExplicitDonors(self._atoms, passNeighborLists, self._extraAtomInfo)
      self._infoString += _ReportTiming(self._verbosity, "fixup explicit doners")

      ################################################################################
      # A Histidine ring Nitrogen is an acceptor if and only if it carries no H or D in
      # this conformer (probe typing marks both as acceptors). Before the symmetry copies
      # take their atom info and before placement reads the flags for rotator candidates;
      # placement marks a Nitrogen whose H it deletes as an acceptor, and the His flip
      # Movers set both Nitrogens for each of their states.
      for a in self._atoms:
        if (a.parent().resname.strip().upper() == 'HIS' and
            a.name.strip().upper() in ('ND1', 'NE2')):
          protonated = False
          for n in passNeighborLists[a]:
            if n.element_is_hydrogen():
              protonated = True
          e = self._extraAtomInfo.getMappingFor(a)
          if e.isAcceptor == protonated:
            e.isAcceptor = not protonated
            self._extraAtomInfo.setMappingFor(a, e)

      ################################################################################
      # Symmetry copies as a fixed environment, in the spatial query before Mover
      # placement (ions, candidate acceptors) but not in the atoms it places Movers on.
      if symmetryCopies is not None:
        self._AddSymmetryCopies(mi, maxISeq, passNeighborLists,
          phantomHydrogenRadius, placedHydrogenDistance)
        self._infoString += _ReportTiming(self._verbosity, "add symmetry copies")

      ################################################################################
      # Get the list of Movers using the _PlaceMovers private function.
      # The list of rotatable hydrogens comes from the global model, not just the current
      # model index.  However, we only place on atoms that are also in self._atoms, which
      # only includes those from the current model index.
      # NOTE: We must do this after fixupExplicitDonors() so that the Hydrogens are properly
      # marked as donors.
      deleteAtoms = self._PlaceMovers(self._atoms, rotatableHydrogens, passNeighborLists, h_parameterization,
                         addFlipMovers, alt)
      self._infoString += _VerboseCheck(self._verbosity, 1,"Inserted "+str(len(self._movers))+" Movers\n")
      self._infoString += _VerboseCheck(self._verbosity, 1,'Marked '+str(len(deleteAtoms))+' atoms for deletion\n')
      self._infoString += _ReportTiming(self._verbosity, "place movers")
      # The environment's Movers adapt it within the pass; their results are not reported.
      self._envMovers = set(m for m in self._movers
        if m.CoarsePositions().atoms[0].i_seq in self._envISeqs)
      self._flipHandled = [a for a in self._flipHandled if a.i_seq not in self._envISeqs]
      self._flipCandidates = [c for c in self._flipCandidates if c[1].i_seq not in self._envISeqs]
      if his is not None:
        self._HisStep(his, self._HisPin, mi, his)

      ################################################################################
      # Add the atoms that were unconditionally marked for deletion during placement
      # to the set of atoms to delete.
      self._deleteMes = self._deleteMes.union(deleteAtoms)
      if symmetryCopies is not None:
        self._PlaceSymmetryCopies(passNeighborLists)

      ################################################################################
      # Initialize the Movers to their starting coarse positions.
      for m in self._movers:
        pr = m.CoarsePositions()
        self._setMoverState(pr, 0)
      self._infoString += _ReportTiming(self._verbosity, "initialize Movers")

      ################################################################################
      # Compute the interaction graph, of which each connected component is a Clique.
      # Get a list of singleton Cliques and a list of other Cliques.  Keep separate lists
      # of the singletons and the groups.
      self._interactionGraph, self._atomMoverLists = InteractionGraph.InteractionGraphAllPairs(self._movers,
        self._extraAtomInfo, probeRadius=self._probeRadius)
      components = cca.connected_components( graph = self._interactionGraph )
      maxLen = 0
      singletonCliques = []   # Each entry is a list of integer indices into models with one entry
      groupCliques = []       # Each entry is a list of integer indices into models with >1 entry
      for c in components:
        if len(c) == 1:
          singletonCliques.append(c)
        else:
          groupCliques.append(c)
        if len(c) > maxLen:
          maxLen = len(c)
      self._infoString += _VerboseCheck(self._verbosity, 1,"Found "+str(len(components))+" Cliques ("+
          str(len(singletonCliques))+" are singletons); largest Clique size = "+
          str(maxLen)+"\n")
      self._infoString += _ReportTiming(self._verbosity, "compute interaction graph")
      # Sort the cliques by size, smallest first, so that if we run out of time we will
      # have optimized as many as possible.
      groupCliques = sorted(groupCliques, key=len, reverse=False)

      # If we've been asked to write an interaction graph file, do so.
      if self._cliqueOutlineFileName:
        with open(self._cliqueOutlineFileName, 'w') as f:
          f.write(self._InteractionKinemage(groupCliques))

      ################################################################################
      # Determine excluded atoms to a specified hop count for each atom that will
      # be moved.  Make a dictionary of lists that includes all atoms in all Movers.

      # Get the set of all atoms that can be returned from all conformations of all Movers.
      moverAtoms = set()
      for m in self._movers:
        for a in m.CoarsePositions().atoms:
          moverAtoms.add(a)
        for a in m.FixUp(0).atoms:
          moverAtoms.add(a)

      # Get the excluded list for each atom in the set, making a dictionary.
      # We go at most 3 hops unless one end of the chain has a hydrogen.
      # Look up excluded atoms by i_seq.
      # By bond count only: the lists serve every pose, so a distance taken with the
      # Movers at state 0 does not hold for them.
      self._excludeDict = {}
      for a in moverAtoms:
        self._excludeDict[a.i_seq] = mmtbx.probe.Helpers.getAtomsWithinNBondsAtAnyDistance(a,
          passNeighborLists, self._extraAtomInfo, self._bondedNeighborDepth, 3)
      self._infoString += _ReportTiming(self._verbosity, "determine excluded atoms")

      ################################################################################
      # Construct dot-sphere cache.
      # This must be done after the phantom Hydrogens have been added so that they will be included.
      dotSphereCache = Helpers.createDotSphereCache(self._probePhil)

      ################################################################################
      # Contruct the DotScorer object we'll use to score the dots.
      self._dotScorer = Helpers.createDotScorer(self._extraAtomInfo, self._probePhil)
      self._infoString += _ReportTiming(self._verbosity, "construct dot scorer")

      ################################################################################
      # Construct C++ optimizer.
      optC = OptimizerC(max(self._verbosity, 1), self._preferenceMagnitude,
                        self._maximumVDWRadius, self._minOccupancy, self._probeRadius, self._probeDensity,
                        self._atoms,
                        self._excludeDict, self._dotScorer, dotSphereCache, self._atomMoverLists,
                        self._spatialQuery, self._extraAtomInfo, self._deleteMes)
      self._infoString += _ReportTiming(self._verbosity, "construct OptimizerC")

      ################################################################################
      # Compute and record the initial score for each Mover in its info
      optC.Initialize(self._movers)
      for m in self._movers:
        score = "{:.2f}".format(optC.GetHighScore(m))
        self._moverInfo[m] += " Initial score: " + score
        # Rounded as printed, so the clique totals below match the report.
        self._initialScore[m] = float(score)

      ################################################################################
      # Call internal methods to optimize the single-element Cliques and then to optimize
      # the multi-element Cliques and then to do independent fine adjustment of all
      # Cliques.  Subclasses should overload the called routines, but the global approach
      # taken here will be the same for all of them.  If we want to change the recipe
      # so that we can do global fine optimization, we'll do that here rather than in the
      # subclasses.

      self._OptimizeCliques(optC, singletonCliques, groupCliques, stopOptimizingAt)
      self._OptimizeFine(optC, singletonCliques + groupCliques, stopOptimizingAt)
      self._optimizerOf = dict((m, optC) for m in self._movers)

      ################################################################################
      # Symmetry copies follow their originals; the Cliques near copies that changed
      # are optimized again.
      if symmetryCopies is not None:
        self._RefreshSymmetryCopies(singletonCliques, groupCliques, dotSphereCache,
          stopOptimizingAt, alt)

      ################################################################################
      # Final scores at the final poses: a stored search score is stale when a Clique
      # partner moved after it was set.
      if his is not None:
        self._HisStep(his, self._HisUnpin, his)
      self._finalScore = dict((m, self._ScoreMover(m, self._dotScorer)) for m in self._movers)

      ################################################################################
      # Flips that are not made: each scored alone, with everything else as placed.
      if self._chooseHis:
        self._SuggestFlips(mi, alt, passNeighborLists, dotSphereCache)
      if his is not None:
        self._HisStep(his, lambda: self._HisRecordPass(mi, alt, passNeighborLists, dotSphereCache, his,
          [[self._interactionGraph.vertex_label(v) for v in c] for c in singletonCliques + groupCliques]))

      ################################################################################
      # Print the final state and score for all Movers
      def _scoreMoverReportClash(self, m, index):
        coarse = m.CoarsePositions()
        score = self._preferenceMagnitude * coarse.preferenceEnergies[index]
        clash = False
        # There may not be as many atoms moved as there are atoms, so use the proper length
        for i in range(len(coarse.positions[index])):
          atom = coarse.atoms[i]
          maxRadiusWithoutProbe = self._extraAtomInfo.getMappingFor(atom).vdwRadius + self._maximumVDWRadius
          res = self._dotScorer.score_dots(atom, self._minOccupancy, self._spatialQuery,
            maxRadiusWithoutProbe, self._probeRadius, self._excludeDict[atom.i_seq],
            self._optimizerOf[m].GetDots(atom.i_seq), self._probeDensity, False, False)
          score += res.totalScore()
          if res.hasBadBump:
            clash = True
        return score, clash

      def _printPose(self, m):
        optC = self._optimizerOf[m]
        description = m.PoseDescription(optC.GetCoarseLocation(m), optC.GetFineLocation(m),
                                        not self._skipBondFixup)

        # If the Mover is a flip of some kind, then the substring "lipped " will be present
        # in the description.
        # When that happens, we check the final state and the other flip
        # state (which is half of the coarse states away) to see if both have clashes or
        # if they are close in energy. If so, then we annotate the output.
        # We add the same number of words to the output string in all cases to make things
        # easier for a program to parse.
        # A His Mover with one orientation (HisPlace) has no other flip state.
        oneOrientation = (isinstance(m, Movers.MoverHisFlip) and
          m.EnabledFlipStates() != 3)
        if "lipped " in description and not oneOrientation:
          coarse = m.CoarsePositions()
          numPositions = len(coarse.positions)
          final = optC.GetCoarseLocation(m)
          other = (final + numPositions//2) % numPositions
          self._setMoverState(coarse, other)
          otherScore, otherBump = _scoreMoverReportClash(self, m, other)
          self._setMoverState(coarse, final)
          finalScore, finalBump = _scoreMoverReportClash(self, m, final)
          if otherBump and finalBump:
            description += " BothClash"
          elif "Unflipped" in description and (
              (otherScore > finalScore) and (otherScore - finalScore <= self._nonFlipPreference)):
            description += " Uncertain"
          else:
            description += " ."
        else:
          description += " ."

        # a fixed His is noted at the end: reduce2 reads the flip state by word position
        if his is not None and his['fullAtoms'] and _HisMid(m) in his['fullAtoms']:
          description += " (state {} fixed)".format(his['fullAtoms'][_HisMid(m)])
        self._infoString += _VerboseCheck(self._verbosity, 1,"  {} final score: {:.2f} pose {}\n".format(
          self._moverInfo[m], self._finalScore[m], description))

      self._infoString += _VerboseCheck(self._verbosity, 1,"BEGIN REPORT: Model "+str(mi)+" Alt '"+alt+"':\n")
      # The pass's own Movers only; a Clique's set and totals leave its environment Movers out.
      sortedGroups = sorted(groupCliques, key=len, reverse=True)
      for g in sortedGroups:
        movers = [self._interactionGraph.vertex_label(i) for i in g]
        movers = [m for m in movers if m not in self._envMovers]
        if len(movers) == 0:
          continue
        self._infoString += _VerboseCheck(self._verbosity, 1," Set of "+str(len(movers))+" Movers:")
        # Sum the initial and final scores across all Movers in the group and report this.
        initial = 0.0
        final = 0.0
        for m in movers:
          initial += self._initialScore[m]
          final += self._finalScore[m]
        self._infoString += _VerboseCheck(self._verbosity, 1," Totals: initial score {:.2f}, final score {:.2f}\n".format(initial, final))
        for m in movers:
          _printPose(self, m)
      self._infoString += _VerboseCheck(self._verbosity, 1," Singleton Movers:\n")
      for s in singletonCliques:
        m = self._interactionGraph.vertex_label(s[0])
        if m not in self._envMovers:
          _printPose(self, m)
      if len(self._envMovers) > 0:
        self._infoString += _VerboseCheck(self._verbosity, 1," Environment Movers (residues without alternate '"+
          alt+"'; written by their own passes): "+str(len(self._envMovers))+"\n")
      self._infoString += _VerboseCheck(self._verbosity, 1,"END REPORT\n")

      ################################################################################
      # Do FixUp on the final coarse orientations.  Set the positions, extra atom info
      # and deletion status for all atoms that have entries for each.
      # This must be done after we print the scores because the print methods move the
      # coarse state to see how much it changed.
      if not self._skipBondFixup:
        self._infoString += _VerboseCheck(self._verbosity, 1,"FixUp on all Movers\n")
        for m in self._movers:
          loc = self._optimizerOf[m].GetCoarseLocation(m)
          self._infoString += _VerboseCheck(self._verbosity, 3,"FixUp on {} coarse location {}\n".format(
          self._moverInfo[m],loc))
          self._doFixup(m.FixUp(loc))
        self._infoString += _ReportTiming(self._verbosity, "fix up Movers")

      ################################################################################
      # The environment back to the state its own passes left (after everything in
      # this pass that moves atoms: placement, the search, symmetry refreshes, FixUp).
      for a, xyz, e in envSaved:
        a.xyz = xyz
        self._extraAtomInfo.setMappingFor(a, e)

      ################################################################################
      # Record this pass's deletions. Only H are deleted (Movers and placement).
      for a in passAtoms:
        if a.i_seq in self._envISeqs:
          continue
        deleted = a in self._deleteMes
        if deleted or a.element_is_hydrogen():
          if a.i_seq not in passDeletions:
            passDeletions[a.i_seq] = (a, [])
          passDeletions[a.i_seq][1].append((alt, deleted))

      #################################################################################
      # Record the fraction of atoms that were calculated and the fraction that were cached.
      for optC in set(self._optimizerOf.values()):
        self._numCalculated += optC.GetNumCalculatedAtoms()
        self._numCached += optC.GetNumCachedAtoms()

      #################################################################################
      # Add any flipped amides and histidines to the lists of flipped amides and histidines.
      for m in self._movers:
        coarse_loc = self._optimizerOf[m].GetCoarseLocation(m)

        if isinstance(m, Movers.MoverAmideFlip) and coarse_loc == 1:
          nitrogen = m.CoarsePositions().atoms[0]
          self._amideFlipsPerformed.append(FlippedMoverInfo(alt, nitrogen))

        # The upper half of the coarse states is flipped (both H, HD1 only, HE2 only).
        # Locked (one-orientation) His were recorded when their Mover was added.
        elif (isinstance(m, Movers.MoverHisFlip) and m.EnabledFlipStates() == 3 and
              coarse_loc >= len(m.CoarsePositions().positions)//2):
          nitrogen = m.CoarsePositions().atoms[0]
          self._hisFlipsPerformed.append(FlippedMoverInfo(alt, nitrogen))

      # This pass's flips of environment residues (Movers, ion lock-down, set flip
      # states) are their own passes' to record.
      self._amideFlipsPerformed[numAmideFlips:] = [f for f in self._amideFlipsPerformed[numAmideFlips:]
        if f.baseAtom.i_seq not in self._envISeqs]
      self._hisFlipsPerformed[numHisFlips:] = [f for f in self._hisFlipsPerformed[numHisFlips:]
        if f.baseAtom.i_seq not in self._envISeqs]

    #################################################################################
    # Each atom's deletion is decided by the last pass that contained it as its own
    # (EnvironmentResidues()): an altloc-specific H by its own pass, a blank H by the
    # last pass that has its residue's letter, or the last pass for a residue without
    # altlocs.
    for i in sorted(passDeletions):
      a, decisions = passDeletions[i]
      if decisions[-1][1]:
        self._hydrogensToDelete.append(a)
      if len(set(d for alt, d in decisions)) > 1:
        self._infoString += _VerboseCheck(self._verbosity, 1,
          "Passes disagree about deleting {} {} (deleted in {}; kept in {}); the last pass decides: {}\n".format(
          _ResNameAndID(a), a.name.strip(),
          " ".join("'{}'".format(alt) for alt, d in decisions if d),
          " ".join("'{}'".format(alt) for alt, d in decisions if not d),
          "deleted" if decisions[-1][1] else "kept"))

  ##################################################################################
  # His charges across conformers (policy and checks in HisCharges)

  def _HisRun(self, pins):
    """Record of one run of the passes. pins: {(His key, conformer): fixed state
    (index in the full Mover)}."""
    self._hisNoMover = {}
    return dict(pins=dict(pins), pinned={}, fullAtoms={}, pinMid={}, passes=[], cliques=[], recs=[])

  def _HisStep(self, his, f, *args):
    """f(*args), code that runs only with reconcileHisCharges. In the first run
    (his['first']) an error becomes _HisFirstRunError; MemoryError and Abort
    propagate as they are."""
    try:
      return f(*args)
    except (MemoryError, Abort):
      raise
    except Exception as e:
      if his.get('first'):
        raise _HisFirstRunError(e) from e
      raise

  def _HisPin(self, mi, his):
    # A fixed His becomes a one-state Mover over its full atom list: every array
    # index is 0 during the search; _HisUnpin restores the full Mover.
    his['pinned'] = {}
    for m in self._movers:
      if not isinstance(m, Movers.MoverHisFlip):
        continue
      atoms = m.CoarsePositions().atoms
      kc = (_HisKey(mi, atoms[0]), _HisLetter([atoms[i] for i in (0, 1, 4, 5)])[0])
      pin = his['pins'].get(kc)
      if pin is None:
        continue
      his['pinMid'][kc] = _HisMid(m)
      m._hisFull = (m._coarsePositions, m._extras, m._deleteMes, m._preferenceEnergies, m._fixUpPositions)
      (m._coarsePositions, m._extras, m._deleteMes, m._preferenceEnergies,
        m._fixUpPositions) = [[x[pin]] for x in m._hisFull]
      his['pinned'][m] = pin
      his['fullAtoms'][_HisMid(m)] = pin

  def _HisUnpin(self, his):
    for m, pin in his['pinned'].items():
      optC = self._optimizerOf[m]
      if optC.GetCoarseLocation(m) != 0:
        raise HisCharges.ContractError(("one-state His Mover not at index 0", self._moverInfo[m],
          optC.GetCoarseLocation(m)))
      (m._coarsePositions, m._extras, m._deleteMes, m._preferenceEnergies,
        m._fixUpPositions) = m._hisFull
      del m._hisFull
      pr = m.CoarsePositions()
      for a, xyz in zip(pr.atoms, pr.positions[pin]):
        if (Helpers.rvec3(a.xyz) - Helpers.rvec3(xyz)).length() > 1e-6:
          raise HisCharges.ContractError("His state {} not at its positions for {}".format(pin, self._moverInfo[m]))
      # reports, FixUp and flip lists read the full Mover at the original index
      self._optimizerOf[m] = _FixedHisOptimizer(optC, pin)
    his['pinned'] = {}

  def _HisRecordPass(self, mi, alt, passNeighborLists, dotSphereCache, his, cliques):
    """This pass's fresh final score per Mover id, the Cliques holding a His, and per
    His NE2 the record decide() takes, with each state's scores, deletions and
    positions (ring and ring H, coarse and fix-up). The pass's own Movers and His
    only: residues without its letter are recorded by their own passes. The
    Cliques stay whole (environment Movers included), so an owned Mover's loss next
    to a pinned environment His is attributed to that His's neighbourhood."""
    own = [m for m in self._movers if m not in self._envMovers]
    his['passes'].append((alt, dict((_HisMid(m), round(self._finalScore[m], 3)) for m in own)))
    his['cliques'] += [set(_HisMid(m) for m in c) for c in cliques
      if any(isinstance(m, Movers.MoverHisFlip) for m in c)]
    hisMovers = dict((m.CoarsePositions().atoms[0].i_seq, m) for m in own
      if isinstance(m, Movers.MoverHisFlip))
    environment = set(self._environmentAtoms)
    for a in self._atoms:
      if not (a.name.strip().upper() == 'NE2' and a.parent().resname.strip().upper() == 'HIS'):
        continue
      if a in environment or a.i_seq in self._envISeqs:
        continue
      letter, partial, ring = _HisRing(a, passNeighborLists)
      hs = dict((n, ring.get(n)) for n in ('HE2', 'HD1'))
      rec = dict(passAlt=alt, key=_HisKey(mi, a), conf=letter, partial=partial, ring=ring)
      if hs['HE2'] is None and hs['HD1'] is None and not letter:
        # a pass that sees the blank ring but not the residue's ring H
        rec['hNotInPass'] = any(x.element_is_hydrogen() and x.name.strip()[1:] in ('D1', 'E2')
          and x not in passNeighborLists for x in a.parent().parent().atoms())
      m = hisMovers.get(a.i_seq)
      if m is None:
        rec.update(mover=None, noMover=self._hisNoMover.get((alt, a.i_seq), 'no Mover'))
        his['recs'].append((rec, hs))
        continue
      pr = m.CoarsePositions()
      optC = self._optimizerOf[m]
      chosen = optC.GetCoarseLocation(m)
      half = len(pr.positions)//2 if m.EnabledFlipStates() == 3 else len(pr.positions)
      names = [x.name.strip() for x in pr.atoms]
      states = []
      for i in range(len(pr.positions)):
        self._setMoverState(pr, i)
        score = self._preferenceMagnitude * pr.preferenceEnergies[i] + self._scorePoseAtoms(
          pr.atoms[:len(pr.positions[i])], passNeighborLists, dotSphereCache)
        dels = pr.deleteMes[i]
        fx = m.FixUp(i)
        states.append(dict(index=i, desc=m.PoseDescription(i, None, False), score=round(score, 3),
          charge=int(not dels[1]) + int(not dels[5]) - 1,
          flipped=(m.EnabledFlipStates() == 3 and i >= half),
          dels={'HE2': bool(dels[1]), 'HD1': bool(dels[5])},
          coarse=dict((n, tuple(p)) for n, p in zip(names, pr.positions[i]) if n in _HIS_RING_NAMES),
          fixup=dict((n, tuple(p)) for n, p in zip(names, fx.positions) if n in _HIS_RING_NAMES)))
      self._setMoverState(pr, chosen)
      sug = None
      for k, v in self._suggestedFlips.items():
        if v is not None and v.atom is a:
          sug = dict(gain=None if v.gain is None else round(v.gain, 3), ion=v.ion is not None)
      rec.update(mover=type(m).__name__, mid=_HisMid(m), chosen=chosen, states=states, suggested=sug)
      his['recs'].append((rec, {'HE2': pr.atoms[1], 'HD1': pr.atoms[5]}))

  def _HisByResidue(self, his):
    """{His key: ({conformer: (record, ring H)}, {conformer: final charge, or with
    reconcileHisStates the state label (HisCharges.state_label())})}. A
    conformer's record is the last pass that has its ring H (the pass whose
    deletion decision holds); a record without ring H only fills a gap."""
    deleted = set(self._hydrogensToDelete)
    byRes = {}
    for rec, hs in his['recs']:
      if rec.get('hNotInPass'):
        continue
      confs = byRes.setdefault(rec['key'], {})
      if any(h is not None for h in hs.values()) or rec['conf'] not in confs:
        confs[rec['conf']] = (rec, hs)
    if self._reconcileHisStates:
      # state mode: the final value is the ring-H state label, not the charge
      final = lambda hs: HisCharges.state_label(dict((n, hs.get(n) is None or hs[n] in deleted)
        for n in ('HD1', 'HE2')))
    else:
      final = lambda hs: sum(1 for h in hs.values() if h is not None and h not in deleted) - 1
    # conformers in letter order, whichever pass recorded first
    byRes = dict((key, dict(sorted(confs.items()))) for key, confs in byRes.items())
    return dict((key, (confs, dict((c, final(hs)) for c, (rec, hs) in confs.items())))
      for key, confs in byRes.items())

  def _HisSave(self, myModel):
    """Everything the passes change for this model: the free result when taken
    after the free run, the starting state when taken before it."""
    return dict(ids=[(a.i_seq, a.id_str()) for a in myModel.atoms()],
      xyz=[tuple(a.xyz) for a in myModel.atoms()],
      info=[probeExt.ExtraAtomInfo(self._extraAtomInfo.getMappingFor(a)) for a in myModel.atoms()],
      deleted=list(self._hydrogensToDelete), amides=list(self._amideFlipsPerformed),
      hisFlips=list(self._hisFlipsPerformed), suggested=dict(self._suggestedFlips),
      parts=list(self._infoString._parts), warnings=self._warningString,
      counts=(self._numCalculated, self._numCached), copyNotes=list(getattr(self, '_copyPassNotes', [])),
      passCopies=list(self._passCopies), absent=set(self._absentCopies),
      environment=list(self._environmentAtoms),
      copyXyz=[tuple(c.xyz) for c, o, rt in self._passCopies],
      copyInfo=[probeExt.ExtraAtomInfo(self._extraAtomInfo.getMappingFor(c)) for c, o, rt in self._passCopies])

  def _HisRestore(self, myModel, saved):
    if [(a.i_seq, a.id_str()) for a in myModel.atoms()] != saved['ids']:
      raise HisCharges.ContractError("atoms changed")
    for a, xyz in zip(myModel.atoms(), saved['xyz']):
      a.xyz = xyz
    for a, info in zip(myModel.atoms(), saved['info']):
      self._extraAtomInfo.setMappingFor(a, probeExt.ExtraAtomInfo(info))
    self._hydrogensToDelete[:] = saved['deleted']
    self._amideFlipsPerformed[:] = saved['amides']
    self._hisFlipsPerformed[:] = saved['hisFlips']
    self._suggestedFlips = dict(saved['suggested'])
    self._infoString._parts[:] = saved['parts']
    self._warningString = saved['warnings']
    self._numCalculated, self._numCached = saved['counts']
    self._copyPassNotes = list(saved['copyNotes'])
    self._passCopies = list(saved['passCopies'])
    self._absentCopies = set(saved['absent'])
    self._environmentAtoms = list(saved['environment'])
    for (c, o, rt), xyz in zip(self._passCopies, saved['copyXyz']):
      c.xyz = xyz
    for (c, o, rt), info in zip(self._passCopies, saved['copyInfo']):
      self._extraAtomInfo.setMappingFor(c, probeExt.ExtraAtomInfo(info))

  def _HisFingerprint(self, myModel):
    """The model's result read without _HisSave (it checks the restore): sites and
    atom info of the model atoms and the symmetry copies, deletions, flips,
    suggested flips, report and warnings."""
    def atoms(pairs):
      out = []
      for key, a in pairs:
        e = self._extraAtomInfo.getMappingFor(a)
        out.append((key, tuple(a.xyz), e.vdwRadius, e.isAcceptor, e.isDonor, e.isDummyHydrogen,
          e.isIon, e.charge))
      return out
    return dict(model=atoms(((a.i_seq, a.id_str()), a) for a in myModel.atoms()),
      copies=atoms(((o.i_seq, str(rt)), c) for c, o, rt in self._passCopies),
      absent=sorted(str(x) for x in self._absentCopies), environment=len(self._environmentAtoms),
      deleted=sorted(a.i_seq for a in self._hydrogensToDelete),
      flips=[str(x) for x in self._amideFlipsPerformed] + [str(x) for x in self._hisFlipsPerformed],
      suggested=sorted((str(k), str(v)) for k, v in self._suggestedFlips.items()),
      report=str(self._infoString), warnings=self._warningString,
      copyNotes=[str(x) for x in getattr(self, '_copyPassNotes', [])])

  def _HisRestoreVerified(self, myModel, saved, fingerprint, cause=None):
    """Restore the free result; a restore that cannot be verified stops the run.
    cause: the error that led to the restore, named in the Sorry."""
    after = '' if cause is None else '; restoring after %s: %s' % (type(cause).__name__, cause)
    try:
      self._HisRestore(myModel, saved)
    except HisCharges.ContractError as e:
      raise Sorry("His charges: the first result could not be restored (%s)%s" % (e, after))
    now = self._HisFingerprint(myModel)
    diff = sorted(k for k in fingerprint if now[k] != fingerprint[k])
    if diff:
      raise Sorry("His charges: the restored first result differs from it (%s)%s" % (", ".join(diff), after))

  def _HisWritten(self, his, key, confs):
    """check_written()'s record per conformer from the final state of the model. A
    pass that sees the blank ring without its ring H writes it too: its letter is
    no fixed conformer (check_written reports it)."""
    deleted = set(self._hydrogensToDelete)
    writer = {}
    for rec, hs in his['recs']:
      if rec['key'] == key:
        w = rec['conf'] if not rec.get('hNotInPass') else 'pass %s without the ring H' % (rec['passAlt'] or 'blank')
        for n, x in rec['ring'].items():
          writer[x.i_seq] = w
    out = {}
    for c, (rec, hs) in confs.items():
      ring = rec['ring']
      out[c] = dict(present=dict((n, ring.get(n) is not None and ring[n] not in deleted) for n in ('HE2', 'HD1')),
        xyz=dict((n, tuple(x.xyz)) for n, x in ring.items() if x not in deleted),
        acceptor=dict((n, self._extraAtomInfo.getMappingFor(ring[n]).isAcceptor) for n in ('NE2', 'ND1') if n in ring),
        writer=dict((n, writer[x.i_seq]) for n, x in ring.items()))
    return out

  hisChargesStrict = False    # his_charges_strict(): after the restore, the error is re-raised

  def _ReconcileHisCharges(self, mi, myModel, ctx):
    """The passes, then the His whose conformers end with different charges:
    decide(), reruns with the reconciled states fixed, trial() and plan() until
    DONE (the last run is the result) or RESTORE_FREE (the saved free result is
    written). Any error after the free result is saved (a failed check,
    HisCharges.ContractError, or any other Exception) keeps the free result,
    restored if a rerun started, with a warning (hisChargesStrict: then
    re-raises it). A MemoryError propagates as it is; a user's Stop (Abort)
    propagates after that restore. Errors in the free run propagate.
    Reported in the info string and getHisCharges(). Returns the seconds the
    reruns took (added to later models' deadline). An error in the code that
    runs only with reconciliation before the first result is saved (the saves,
    _HisPin, _HisUnpin, _HisRecordPass): _HisFirstRunFallback."""
    ts = time.time()    # before the first save: a slow failed save is credited too
    try:
      start = self._HisSave(myModel)
    except (MemoryError, Abort):
      raise
    except Exception as e:
      return self._HisFirstRunFallback(mi, myModel, ctx, None, e, ts)
    freeHis = self._HisRun({})
    freeHis['first'] = True
    t0 = time.time()
    try:
      self._OptimizeModel(mi, myModel, ctx, freeHis)
    except _HisFirstRunError as e:
      return self._HisFirstRunFallback(mi, myModel, ctx, start, e.error, t0)
    t1 = time.time()
    try:
      # an error here: the rerun repeats a finished first run (rare; one rule for all)
      free = self._HisSave(myModel)
    except (MemoryError, Abort):
      raise
    except Exception as e:
      return self._HisFirstRunFallback(mi, myModel, ctx, start, e, t0)
    freeRes = None
    first = {}
    cases = []
    run = dict(code=HisCharges.NO_RERUN, detail=None, reruns=0, cap=0, freeSeconds=t1 - t0,
      rerunSeconds=0.0, backstopSeconds=None, unmatched=0)
    st = freePrint = None
    failure = error = None
    try:
      freeRes = self._HisByResidue(freeHis)
      for key, (confs, final) in sorted(freeRes.items()):
        if len(confs) < 2 or len(set(final.values())) < 2:
          continue
        lacking = _HisBlankRingLastPassLacksH([r for r, hs in freeHis['recs']], key)
        if lacking is not None:
          # the final writer of the blank ring is a pass without its ring H
          d = dict(kind='unchanged', decision='unchanged: blank ring last written by pass %s, '
            'which lacks the ring H' % (lacking or 'blank'))
        else:
          d = HisCharges.decide(dict((c, r) for c, (r, hs) in confs.items()), self._hisBudget,
            0.05 * self._preferenceMagnitude, same_state=self._reconcileHisStates)
        cases.append(dict(key=key, free=final, final=dict(final), first=d))
        if d['kind'] == 'reconciled':
          first[key] = d
      if first:
        freePrint = self._HisFingerprint(myModel)
        st = HisCharges.new_plan(sorted(first))
        run.update(cap=HisCharges.run_cap(st['n']),
          backstopSeconds=HisCharges.WALL_BACKSTOP_FACTOR * (t1 - t0))
        out, his, now, pins = self._HisReruns(mi, myModel, ctx, st, start, free, first, freeHis,
          freeRes, cases, t0, t1, run)
        if out == HisCharges.DONE:
          self._HisCheckWritten(st, his, now, pins, freeRes, cases)
          run['code'] = HisCharges.RECONCILED
    except MemoryError:
      raise
    except Abort as e:
      # a user's Stop: a caller that catches it holds the first result
      if run['reruns'] > 0:
        self._HisRestoreVerified(myModel, free, freePrint, e)
      raise
    except Exception as e:
      error = e
      failure = str(e) if isinstance(e, HisCharges.ContractError) else '%s: %s' % (type(e).__name__, e)
      out = HisCharges.RESTORE_FREE
      run.update(code=HisCharges.INTERNAL, detail=failure)
      if run['backstopSeconds'] is None:
        # the description prints it for every outcome but NO_RERUN
        run['backstopSeconds'] = HisCharges.WALL_BACKSTOP_FACTOR * (t1 - t0)
      # His with different charges that decide() did not reach
      done = set(c['key'] for c in cases)
      for key, (confs, final) in sorted((freeRes or {}).items()):
        if key not in done and len(confs) > 1 and len(set(final.values())) > 1:
          cases.append(dict(key=key, free=final, final=dict(final),
            first=dict(kind='reconciled', decision='not decided')))
    if not first and error is None:
      self._HisReportGuarded(mi, cases, None, None, run)
      return 0.0
    if out != HisCharges.DONE and run['reruns'] > 0:
      self._HisRestoreVerified(myModel, free, freePrint, error)
    if out != HisCharges.DONE:
      # nothing of the reruns is written: no His fixed, the free charges
      cases[:] = [c for c in cases if c['first']['kind'] != 'changed']
      for c in cases:
        c['final'] = dict(c['free'])
        for x in ('fixed', 'checked', 'skipped'):
          c.pop(x, None)
    if failure is not None:
      self._warningString += ("Warning: His charges: reconciliation stopped by an error, the "
        "first result is kept: %s: %s\n" % (type(error).__name__, error))
    run['rerunSeconds'] = time.time() - t1
    self._HisReportGuarded(mi, cases, st, out, run)
    if error is not None and self.hisChargesStrict:
      raise error
    return run['rerunSeconds']

  def _HisFirstRunFallback(self, mi, myModel, ctx, start, error, t0):
    """The first run failed in code that runs only with reconciliation, maybe
    partway through a pass (_HisRecordPass sets each His state in turn, so Movers
    and atoms can be left changed): the starting state is restored (start None:
    nothing ran) and the passes run again without reconciliation, with the search
    time the first run was allowed (from the end of the restore; t0 is taken before
    the first save), so the model is
    the one written without it. Warning, code FIRST_RUN_ERROR, no His cases; strict
    mode re-raises error. Returns the seconds of the failed run and the restore
    (added to later models' deadline; the fallback run is the run made without
    reconciliation too)."""
    if start is not None:
      try:
        self._HisRestore(myModel, start)
      except (MemoryError, Abort):
        raise
      except Exception as e:
        raise Sorry("His charges: the starting state could not be restored (%s: %s); restoring after %s: %s"
          % (type(e).__name__, e, type(error).__name__, error))
    tr = time.time()
    # the failed run, save or restore does not use up the search time
    runCtx = dict(ctx, stopOptimizingAt=tr + (ctx['stopOptimizingAt'] - t0))
    self._OptimizeModel(mi, myModel, runCtx)
    detail = '%s: %s' % (type(error).__name__, error)
    self._warningString += ("Warning: His charges: an error in the first run, the model was "
      "optimized again without reconciliation: %s\n" % detail)
    run = dict(code=HisCharges.FIRST_RUN_ERROR, detail=detail, reruns=0, cap=0, freeSeconds=tr - t0,
      rerunSeconds=time.time() - tr, backstopSeconds=None, unmatched=0)
    self._HisReportGuarded(mi, [], None, None, run)
    if self.hisChargesStrict:
      raise error
    return tr - t0

  def _HisReruns(self, mi, myModel, ctx, st, start, free, first, freeHis, freeRes, cases, t0, t1, run):
    """The reruns until plan() ends. Returns (outcome, last run's record, its
    _HisByResidue, its pins); run gets the outcome code, detail and reruns.
    Limits (HisCharges): run_cap(n) reruns; for each rerun, the search time the
    free run was allowed (not its duration); the wall-time backstop. A limit restores the free result,
    as does a rerun whose limited Cliques (_LimitedCliques) differ from the free run's."""
    freeFinal = dict((k, f) for k, (confs, f) in freeRes.items() if len(confs) > 1)
    freeMids = set(m for alt, sc in freeHis['passes'] for m in sc)
    firstCost = dict((k, d['costs'][str(d['target'])]['total']) for k, d in first.items())
    mids = lambda res, k: set(r['mid'] for r, hs in res[k][0].values() if r.get('mid'))
    searchTime = ctx['stopOptimizingAt'] - t0
    backstopEnd = t1 + run['backstopSeconds']
    rerunCost = {}
    his = now = pins = None
    if t1 > ctx['stopOptimizingAt']:
      run.update(code=HisCharges.SEARCH_DEADLINE, detail='the first run reached it')
      return HisCharges.RESTORE_FREE, his, now, pins
    while True:
      if run['reruns'] >= run['cap']:
        run.update(code=HisCharges.RUN_CAP, detail='%d rerun(s)' % run['reruns'])
        out = HisCharges.RESTORE_FREE
        break
      if time.time() > backstopEnd:
        run.update(code=HisCharges.WALL_BACKSTOP, detail='%d x the first run' % HisCharges.WALL_BACKSTOP_FACTOR)
        out = HisCharges.RESTORE_FREE
        break
      keys = sorted(st['accepted'] | set(st['trial']))
      pins = dict(((k, c), p['to']) for k in keys for c, p in first[k]['per'].items())
      run['reruns'] += 1    # from here the model may differ from the free result
      self._HisRestore(myModel, start)
      his = self._HisRun(pins)
      runCtx = dict(ctx, stopOptimizingAt=time.time() + searchTime)
      self._OptimizeModel(mi, myModel, runCtx, his)
      if time.time() > runCtx['stopOptimizingAt']:
        # Cliques may have been optimized as singletons: not comparable
        run.update(code=HisCharges.SEARCH_DEADLINE, detail='rerun %d reached it' % run['reruns'])
        out = HisCharges.RESTORE_FREE
        break
      if _LimitedCliques(self._warningString) != _LimitedCliques(free['warnings']):
        # a fixed His lowers a Clique's work count, and the time allowance varies:
        # a Clique searched differently from the free run is not comparable
        run.update(code=HisCharges.SEARCH_DEADLINE,
          detail='rerun %d limited other Cliques than the first run' % run['reruns'])
        out = HisCharges.RESTORE_FREE
        break
      # contract: same atoms and Mover ids as the free run; each fixed His was a
      # one-state Mover over its own free-run Mover's full atom list
      if [(a.i_seq, a.id_str()) for a in myModel.atoms()] != free['ids']:
        raise HisCharges.ContractError("atoms changed in the rerun")
      if len(his['fullAtoms']) != len(pins) or sorted(his['pinMid']) != sorted(pins):
        raise HisCharges.ContractError(("fixed His without a one-state Mover", sorted(pins), sorted(his['pinMid'])))
      now = self._HisByResidue(his)
      for (k, c), mid in sorted(his['pinMid'].items()):
        if mid not in freeMids or freeRes[k][0].get(c, ({},))[0].get('mid') != mid:
          raise HisCharges.ContractError(("one-state Mover is not the free run's", HisCharges.key_text(k), c, mid))
      for k in keys:
        if k not in now or not mids(now, k):
          raise HisCharges.ContractError("His %s fixed but not recorded in the rerun" % HisCharges.key_text(k))
        confs, final = now[k]
        HisCharges.check_fixed(k, dict((c, p) for (x, c), p in pins.items() if x == k),
          dict((c, r['chosen']) for c, (r, hs) in confs.items() if r.get('states') is not None), final)
      nowFinal = dict((k, f) for k, (confs, f) in now.items() if len(confs) > 1)
      newMismatches = HisCharges.new_mismatches_states if self._reconcileHisStates else HisCharges.new_mismatches
      mm, missing = newMismatches(freeFinal, nowFinal, set(keys))
      ev = HisCharges.trial(freeHis['passes'], his['passes'], dict((k, mids(now, k)) for k in keys),
        freeHis['cliques'] + his['cliques'], st['trial'], dict((k, mids(now, k)) for k in mm),
        missing=missing, budget=float('inf') if self._hisBudget is None else self._hisBudget)
      run['unmatched'] = sum(len(comp['cost']['unmatched']) for comp in ev['components']) + \
        len(ev['unattributedUnmatched'])
      for comp in ev['components']:
        for k in comp['keys']:
          rerunCost[k] = dict(total=comp['cost']['total'], his=comp['cost']['his'],
            neighbours=comp['cost']['neighbours'], keys=comp['keys'])
      out = HisCharges.plan(st, ev, firstCost)
      if out == HisCharges.FAILED:
        # free_saved is True: plan() never returns FAILED
        raise HisCharges.ContractError(("plan() failed with the free result saved", st['reason']))
      if out == HisCharges.RESTORE_FREE:
        run.update(code=HisCharges.REJECTED, detail=st['reason'])
      if out != HisCharges.RERUN:
        break
    for c in cases:
      if c['key'] in rerunCost:
        c['rerun'] = rerunCost[c['key']]
    return out, his, now, pins

  def _HisCheckWritten(self, st, his, now, pins, freeRes, cases):
    """DONE: check_written() for every accepted His; any error raises
    ContractError before anything is marked. Then the marks (fixed, checked,
    skipped), and the final charges of the written run for every reported His
    and every His whose charges changed in it."""
    results = {}
    errors = []
    for k in sorted(st['accepted']):
      confs, final = now[k]
      # expected positions: the free run's full Mover at the original index
      fixed = {}
      for (x, c), p in pins.items():
        if x == k:
          s = freeRes[k][0][c][0]['states'][p]
          fixed[c] = dict(dels=s['dels'], coarse=s['coarse'], fixup=s['fixup'])
      results[k] = r = HisCharges.check_written(fixed, self._HisWritten(his, k, confs), self._skipBondFixup)
      if r['errors']:
        errors.append("His %s written state: %s" % (HisCharges.key_text(k), "; ".join(r['errors'])))
    if errors:
      raise HisCharges.ContractError(" | ".join(errors))
    for c in cases:
      r = results.get(c['key'])
      if r is not None:
        c.update(checked=r['checked'], skipped=r['skipped'], fixed=True)
    reported = set(c['key'] for c in cases)
    for c in cases:
      if c['key'] in now:
        c['final'] = dict(now[c['key']][1])
    for k in sorted(set(now) & set(freeRes) - reported):
      f, n = freeRes[k][1], now[k][1]
      if f != n:
        cases.append(dict(key=k, free=dict(f), final=dict(n), first=dict(kind='changed',
          decision='changed in the rerun (not fixed)')))

  def _HisReport(self, mi, cases, st, out, run):
    """Per His: decision text and category (fixed, metal, cost, rejected,
    restored, unchanged, changed; from the structure, for the qmready hint); the
    model record for getHisCharges()."""
    lines = []
    for c in cases:
      key = c['key']
      d = c['first']
      if d['kind'] == 'kept':
        # decide()'s metal, ion and user-protonation results have no target; a cost has one
        c['category'] = 'metal' if d.get('target') is None else 'cost'
        c['decision'] = d['decision']
      elif d['kind'] != 'reconciled':
        c['category'] = d['kind']
        c['decision'] = d['decision']
      elif c.get('fixed'):
        c['category'] = 'fixed'
        c['decision'] = 'reconciled to %s %s (first stage: %s)' % ('state' if self._reconcileHisStates else 'charge',
          HisCharges.charge_text(list(c['final'].values())[0]), d['decision'])
        cost = d['costs'][str(d['target'])]
        rerun = (c.get('rerun') or {}).get('total')
        if self._reconcileHisStates and max(cost['total'], rerun or 0.0) > HisCharges.BUDGET:
          # no cost limit in state mode: a large screen or rerun cost is reported, not refused
          c['warning'] = ("His %s: one state for all conformers (%s) costs %.2f (%s; rerun %s); the "
            "conformers' own states were %s" % (HisCharges.key_text(key), d['target'], cost['total'],
            ', '.join('%s %.2f' % (x or 'blank', p['cost']) for x, p in sorted(cost['per'].items())),
            'none' if rerun is None else '%.2f' % rerun,
            ', '.join('%s %s' % (x or 'blank', c['free'][x]) for x in sorted(c['free']))))
          self._warningString += 'Warning: ' + c['warning'] + '\n'
      elif st is not None and key in st['rejected']:
        c['category'] = 'rejected'
        c['decision'] = 'kept: rerun failed (%s)' % st['rejected'][key]
      else:
        c['category'] = 'restored'
        c['decision'] = 'kept: %s' % _HisOutcomeText(run)
      lines.append("His charges: His %s conformers %s -> %s\n" % (HisCharges.key_text(key),
        " ".join("%s:%s" % (x or 'blank', HisCharges.charge_text(c['final'][x])) for x in sorted(c['final'])), c['decision']))
    lines.append("His charges: model index %d: %s\n" % (mi, _HisOutcomeText(run)))
    for l in lines:
      self._infoString += _VerboseCheck(self._verbosity, 3, l)
    rec = dict(model=mi, outcome=out, cases=cases, maxRuns=HisCharges.max_runs(st['n']) if st else 0,
      rejected=dict(st['rejected']) if st else {}, reason=_HisOutcomeText(run))
    rec.update(run)
    self._hisCharges.append(rec)

  def _HisReportGuarded(self, mi, cases, st, out, run):
    """_HisReport; the model is final by then, so an error in it (MemoryError
    and Abort apart) only warns, and getHisCharges() gets a record without cases
    (reportError). Strict mode re-raises."""
    n = len(self._hisCharges)
    try:
      self._HisReport(mi, cases, st, out, run)
    except (MemoryError, Abort):
      raise
    except Exception as e:
      if self.hisChargesStrict:
        raise
      detail = '%s: %s' % (type(e).__name__, e)
      self._warningString += "Warning: His charges: the report could not be written: %s\n" % detail
      del self._hisCharges[n:]
      rec = dict(model=mi, outcome=out, cases=[], maxRuns=0, rejected={}, reason='report not written',
        reportError=detail)
      rec.update(run)
      self._hisCharges.append(rec)

  def getHisCharges(self):
    """His whose conformers ended the free run with different charges, per model:
    dict(model, outcome (None, DONE, RESTORE_FREE), code (HisCharges outcome code),
    detail, reason (text from code and detail), reruns, cap, maxRuns, freeSeconds,
    reportError (only when the report failed: the error; cases is then empty),
    rerunSeconds, backstopSeconds, unmatched, rejected {key: why}, cases
    [dict(key, free, final, first (decide()), decision, category, rerun (cost),
    fixed, checked, skipped, warning (reconcileHisStates: a fixed state that costs more
    than HisCharges.BUDGET))]). free and final are {conformer: charge}, or {conformer:
    state label} with reconcileHisStates (HisCharges.state_label())."""
    return self._hisCharges

  ##################################################################################
  # Symmetry copies as a fixed environment

  def _AddSymmetryCopies(self, mi, maxISeq, passNeighborLists, phantomHydrogenRadius,
                         placedHydrogenDistance):
    """Copies of this pass's atoms at their originals' current state, with their atom
    info and bonded neighbors, and Phantom Hydrogens on copied waters, added to the
    spatial query and to self._environmentAtoms (not to self._atoms, so that no Movers
    are placed on them)."""
    self._environmentAtoms = []
    self._passCopies = []
    self._absentCopies = set()
    self._crossNeighbors = []       # (atom, bonded atom), added after Mover placement
    originals = {}
    for a in self._atoms:
      if a.i_seq <= maxISeq:
        originals[a.i_seq] = a
    for c, o, rt in self._symmetryCopies.models.get(mi, []):
      if o.i_seq in originals:
        self._passCopies.append((c, originals[o.i_seq], rt))
    if len(self._passCopies) == 0:
      return
    uc = self._symmetryCopies.crystal_symmetry.unit_cell()
    for c, o, rt in self._passCopies:
      c.xyz = uc.orthogonalize(rt * uc.fractionalize(o.xyz))
      self._extraAtomInfo.setMappingFor(c, probeExt.ExtraAtomInfo(self._extraAtomInfo.getMappingFor(o)))
    copyNeighbors, self._crossNeighbors = _SymmetryCopyBonds(originals, self._passCopies,
      passNeighborLists, self._crossBonds, uc)
    passNeighborLists.update(copyNeighbors)

    for c, o, rt in self._passCopies:
      self._spatialQuery.add(c)
      self._environmentAtoms.append(c)
    # Phantom Hydrogens on copied waters, as on the model's: none on the copy of a
    # water that carries H or D in this pass.
    numSkipped = 0
    for c, o, rt in self._passCopies:
      if c.element == 'O' and common_residue_names_get_class(name=c.parent().resname) == "common_water":
        if o.i_seq in self._watersWithH:
          numSkipped += 1
          self._infoString += _VerboseCheck(self._verbosity, 3,
            "Not adding phantom Hydrogens on the copy of {} ({}): it carries H/D\n".format(
            _ResNameAndID(o), rt))
          continue
        newPhantoms = Helpers.getPhantomHydrogensFor(self._copyPhantomISeq, c, self._spatialQuery,
          self._extraAtomInfo, self._minOccupancy, False, phantomHydrogenRadius, placedHydrogenDistance)
        for p in newPhantoms:
          self._spatialQuery.add(p)
          self._extraAtomInfo.setMappingFor(p, probeExt.ExtraAtomInfo(phantomHydrogenRadius, False, True, True))
          passNeighborLists[p] = [c]
          self._environmentAtoms.append(p)
    if numSkipped > 0:
      self._infoString += _VerboseCheck(self._verbosity, 1,"No phantom Hydrogens on "+
        str(numSkipped)+" symmetry copies of waters that carry H or D\n")

  def _PlaceSymmetryCopies(self, passNeighborLists):
    """After Mover placement: the bonds from originals to copies (and the reverse of
    copies' bonds to atoms at special positions), the copies of atoms
    deleted during placement out of the spatial query, the copies into self._atoms."""
    for o, c in self._crossNeighbors:
      if c not in passNeighborLists[o]:
        passNeighborLists[o].append(c)
    for c, o, rt in self._passCopies:
      if o in self._deleteMes:
        self._spatialQuery.remove(c)
        self._absentCopies.add(c)
    self._atoms.extend(self._environmentAtoms)

  def _FinalOriginalStates(self):
    """Position, atom info and deletion of the atoms that FixUp of the Movers' current
    states changes, by i_seq (the state written out)."""
    xyz = {}
    info = {}
    deleted = {}
    if not self._skipBondFixup:
      for m in self._movers:
        fixUp = m.FixUp(self._optimizerOf[m].GetCoarseLocation(m))
        for i, a in enumerate(fixUp.atoms):
          if i < len(fixUp.positions):
            xyz[a.i_seq] = fixUp.positions[i]
          if i < len(fixUp.extraInfos):
            info[a.i_seq] = fixUp.extraInfos[i]
          if i < len(fixUp.deleteMes):
            deleted[a.i_seq] = fixUp.deleteMes[i]
    return xyz, info, deleted

  def _UpdateCopies(self):
    """Move the copies to their originals' final state.
    :return: list of (previous position, copy) for the copies that changed."""
    xyz, info, deleted = self._FinalOriginalStates()
    uc = self._symmetryCopies.crystal_symmetry.unit_cell()
    changed = []
    for c, o, rt in self._passCopies:
      site = uc.orthogonalize(rt * uc.fractionalize(xyz.get(o.i_seq, o.xyz)))
      e = info.get(o.i_seq)
      if e is None:
        e = self._extraAtomInfo.getMappingFor(o)
      gone = deleted.get(o.i_seq, o in self._deleteMes)
      wasGone = c in self._absentCopies
      if ((Helpers.rvec3(site) - Helpers.rvec3(c.xyz)).length() <= 0.01 and
          _SameAtomInfo(e, self._extraAtomInfo.getMappingFor(c)) and gone == wasGone):
        continue
      changed.append((c.xyz, c))
      if not wasGone:
        self._spatialQuery.remove(c)
      c.xyz = site
      self._extraAtomInfo.setMappingFor(c, probeExt.ExtraAtomInfo(e))
      if gone:
        self._absentCopies.add(c)
      else:
        self._absentCopies.discard(c)
        self._spatialQuery.add(c)
    return changed

  def _CliquesNear(self, cliques, sites):
    """The Cliques with an atom, in any coarse position, within interaction range of
    one of the sites."""
    reach = 2 * self._maximumVDWRadius + 2 * self._probeRadius + 0.5
    grid = {}
    for s in sites:
      _AddToGrid(grid, _GridKey(s, reach), Helpers.rvec3(s))
    result = []
    for c in cliques:
      near = False
      for v in c:
        m = self._interactionGraph.vertex_label(v)
        for pos in m.CoarsePositions().positions:
          for x in pos:
            k = _GridKey(x, reach)
            for d in itertools.product((-1, 0, 1), repeat=3):
              for s in grid.get((k[0]+d[0], k[1]+d[1], k[2]+d[2]), []):
                if (s - Helpers.rvec3(x)).length() < reach:
                  near = True
                  break
              if near: break
            if near: break
          if near: break
        if near: break
      if near:
        result.append(c)
    return result

  def _Poses(self):
    """The current pose of each Mover: (OptimizerC, coarse index, fine index)."""
    result = {}
    for m in self._movers:
      optC = self._optimizerOf[m]
      result[m] = (optC, optC.GetCoarseLocation(m), optC.GetFineLocation(m))
    return result

  def _SetPose(self, m, pose):
    optC, coarse, fine = pose
    self._setMoverState(m.CoarsePositions(), coarse)
    if fine >= 0:
      self._setMoverState(m.FinePositions(coarse), fine)
    self._optimizerOf[m] = optC

  def _ScoreMover(self, m, dotScorer):
    """Score of a Mover's current pose against everything around it, copies included."""
    optC = self._optimizerOf[m]
    coarse = optC.GetCoarseLocation(m)
    fine = optC.GetFineLocation(m)
    pr = m.CoarsePositions()
    index = coarse
    if fine >= 0:
      pr = m.FinePositions(coarse)
      index = fine
    score = self._preferenceMagnitude * pr.preferenceEnergies[index]
    for i in range(len(pr.positions[index])):
      a = pr.atoms[i]
      if a in self._deleteMes:
        continue
      maxRadiusWithoutProbe = self._extraAtomInfo.getMappingFor(a).vdwRadius + self._maximumVDWRadius
      score += dotScorer.score_dots(a, self._minOccupancy, self._spatialQuery,
        maxRadiusWithoutProbe, self._probeRadius, self._excludeDict[a.i_seq],
        optC.GetDots(a.i_seq), self._probeDensity, False, False).totalScore()
    return score

  def _RefreshSymmetryCopies(self, singletonCliques, groupCliques, dotSphereCache,
                             stopOptimizingAt, alt):
    """Copies follow their originals; the Cliques near copies that changed are optimized
    again, in a new OptimizerC, at most 3 passes in all. Movers still changing after the
    last pass keep the combination of their states that scores best with matching copies;
    so do the Movers a pass changed when the time limit overtook it, between their previous
    and new states, and the refresh stops."""
    cliques = singletonCliques + groupCliques
    poses = [self._Poses()]
    passNotes = []
    lastRedo = []
    lastChanged = []
    changed = self._UpdateCopies()
    moved = [(1, self._ChangedSites(changed))]  # (pass after which copies moved, sites)
    while len(changed) > 0 and len(poses) < 3:
      # Past the time limit a pass would solve its Cliques as singletons, in place of
      # the main search's solutions.
      if time.time() > stopOptimizingAt:
        passNotes.append("stopped at pass {} due to time limit".format(len(poses) + 1))
        changed = []
        break
      redo = self._CliquesNear(cliques, moved[-1][1])
      movers = [self._interactionGraph.vertex_label(v) for c in redo for v in c]
      if len(movers) > 0:
        for m in movers:
          self._setMoverState(m.CoarsePositions(), 0)
        optC = OptimizerC(max(self._verbosity, 1), self._preferenceMagnitude,
                          self._maximumVDWRadius, self._minOccupancy, self._probeRadius, self._probeDensity,
                          self._atoms,
                          self._excludeDict, self._dotScorer, dotSphereCache, self._atomMoverLists,
                          self._spatialQuery, self._extraAtomInfo, self._deleteMes)
        optC.Initialize(movers)
        self._OptimizeCliques(optC, [c for c in redo if len(c) == 1],
          [c for c in redo if len(c) > 1], stopOptimizingAt)
        self._OptimizeFine(optC, redo, stopOptimizingAt)
        for m in movers:
          self._optimizerOf[m] = optC
      poses.append(self._Poses())
      lastChanged = [m for m in movers if poses[-1][m][1:] != poses[-2][m][1:]]
      lastRedo = redo
      passNotes.append("pass {}: {} of {} Cliques re-optimized, {} Movers changed".format(
        len(poses), len(redo), len(cliques), len(lastChanged)))
      self._infoString += _VerboseCheck(self._verbosity, 3, "Symmetry copies "+passNotes[-1]+
        "".join("\n  "+self._MoverName(m) for m in lastChanged)+"\n")
      if time.time() > stopOptimizingAt:
        # The time limit overtook the pass, which may have solved Cliques as singletons:
        # the Movers it changed keep the best of their previous and new states.
        if len(lastChanged) > 0:
          text = self._ResolveUnsettled(lastChanged, lastRedo, poses[-2:], cliques)
          passNotes[-1] += ", time limit reached: best of previous and new states of " + text
        changed = []
        break
      changed = self._UpdateCopies()
      moved.append((len(poses), self._ChangedSites(changed)))

    unsettled = ""
    if len(changed) > 0 and len(lastChanged) > 0:
      unsettled = self._ResolveUnsettled(lastChanged, lastRedo, poses, cliques)

    text = "{} passes".format(len(poses)) if len(poses) > 1 else "1 pass"
    if len(passNotes) > 0:
      text += " (" + "; ".join(passNotes) + ")"
    self._copyPassNotes.append((alt, text, unsettled))

  def _ChangedSites(self, changed):
    """Previous and current sites of the copies that _UpdateCopies() changed."""
    sites = []
    for old, c in changed:
      sites.append(old)
      sites.append(c.xyz)
    return sites

  def _MoverName(self, m):
    return self._moverInfo[m].split(" Initial score:")[0]

  def _ResolveUnsettled(self, unsettled, lastRedo, poses, cliques):
    """For Movers whose state changed in the last pass, try the combinations of the
    states they had in the passes of poses (or the passes' combinations if there are more
    than 64), each with copies rebuilt to match, and keep the one with the best summed score
    of the Cliques involved.
    :return: description of the unsettled Movers."""
    options = []
    for m in unsettled:
      seen = []
      for p in poses:
        if p[m][1:] not in [s[1:] for s in seen]:
          seen.append(p[m])
      options.append(seen)
    n = 1
    for o in options:
      n *= len(o)
    if n <= 64:
      combos = list(itertools.product(*options))
    else:
      combos = []
      for p in poses:
        combo = tuple(p[m] for m in unsettled)
        if [x[1:] for x in combo] not in [[x[1:] for x in c] for c in combos]:
          combos.append(combo)
    scored = []
    for c in cliques:
      movers = [self._interactionGraph.vertex_label(v) for v in c]
      if c in lastRedo or any(m in unsettled for m in movers):
        scored.extend(movers)
    best = None
    for combo in combos:
      for m, pose in zip(unsettled, combo):
        self._SetPose(m, pose)
      self._UpdateCopies()
      total = 0.0
      for m in scored:
        total += self._ScoreMover(m, self._dotScorer)
      if best is None or total > best[0]:
        best = (total, combo)
    for m, pose in zip(unsettled, best[1]):
      self._SetPose(m, pose)
    self._UpdateCopies()
    return "{} (best of {} combinations, score {:.2f})".format(
      "; ".join(self._MoverName(m) for m in unsettled), len(combos), best[0])

  def _SummarizeCopies(self, mi):
    n = len(self._symmetryCopies.models.get(mi, []))
    if n == 0:
      return
    parts = []
    for alt, text, unsettled in self._copyPassNotes:
      if len(self._copyPassNotes) > 1:
        text = "alternate '{}': {}".format(alt, text)
      parts.append(text)
    unsettled = "; ".join("'{}' {}".format(alt, u) if len(self._copyPassNotes) > 1 else u
      for alt, text, u in self._copyPassNotes if u)
    line = "Symmetry copies (fixed environment): model {}, {} copy atoms, {}; unsettled: {}\n".format(
      mi + 1, n, "; ".join(parts), unsettled if unsettled else "none")
    self._copySummary += line
    self._infoString += _VerboseCheck(self._verbosity, 1, line)

  def getSymmetryCopies(self):
    """
      :return: the symmetry copies of the last alternate pass as (copy, original, rt_mx,
      present) tuples; present is False for copies of deleted originals.
    """
    return [(c, o, rt, c not in self._absentCopies) for c, o, rt in self._passCopies]

  def getSymmetryCopiesSummary(self):
    """
      :return: one line per model about the symmetry copies: the number of copy atoms,
      the passes and the Movers that did not settle; empty without copies.
    """
    return self._copySummary

  def _OptimizeCliques(self, optC, singletonCliques, groupCliques, stopOptimizingAt):
    """Coarse optimization of the Cliques."""
    # Do coarse optimization on the singleton Movers.  Record the selected coarse
    # index.  We do not check for the time limit here because singleton optimization
    # is very fast.
    for s in singletonCliques:
      mover = self._interactionGraph.vertex_label(s[0])
      (bestScore, infoString) = optC.OptimizeSingleMoverCoarse(mover)
      self._infoString += infoString
      ret = bestScore
      self._infoString += _VerboseCheck(self._verbosity, 1, "Singleton optimized with score {:.2f}\n".format(ret))
    self._infoString += _ReportTiming(self._verbosity, "optimize singletons (coarse)")

    # Do coarse optimization on the multi-Mover Cliques.
    # Sort the interaction-graph edges into their Cliques in one pass, in graph
    # order.  (Building each Clique's subgraph with _subsetGraph() scanned every
    # vertex and edge of the whole graph once per Clique: quadratic in the
    # number of Movers.)
    cliqueOf = {}
    for ci, g in enumerate(groupCliques):
      for v in g:
        cliqueOf[self._interactionGraph.vertex_label(v)] = ci
    cliqueEdges = [[] for g in groupCliques]
    for e in self._interactionGraph.edges():
      source = self._interactionGraph.vertex_label(self._interactionGraph.source(e))
      ci = cliqueOf.get(source)
      if ci is not None:
        cliqueEdges[ci].append((source,
          self._interactionGraph.vertex_label(self._interactionGraph.target(e))))
    for ci, g in enumerate(groupCliques):
      remainingSeconds = stopOptimizingAt - time.time()
      perClique = False
      reason = None
      if remainingSeconds <= 0:
        reason = 'optimized as singletons, the total time limit had passed'
      elif self._cliqueSeconds is not None and self._cliqueSeconds < remainingSeconds:
        remainingSeconds = self._cliqueSeconds
        perClique = True

      # Find all of the Movers in the clique so that we know which ones to capture the state for.
      # connected_components() lists them in graph order, as _subsetGraph() did.
      movers = [self._interactionGraph.vertex_label(i) for i in g]
      index = dict((m, i) for i, m in enumerate(movers))

      # Turn the graph edges into a versa array that holds the edges. The first index is the
      # number of edge and the second is 2D, listing the index of each mover.
      edges = flex.int(flex.grid(len(cliqueEdges[ci]), 2))
      for i, (source, target) in enumerate(cliqueEdges[ci]):
        edges[(i,0)] = index[source]
        edges[(i,1)] = index[target]

      # A search above the work limit is not started: with 0 s the C++ treats the
      # Clique as singletons.
      if reason is None and self._maxCliqueWork is not None:
        work = _CliqueSearchWork(len(movers),
          [(index[s], index[t]) for s, t in cliqueEdges[ci]],
          [len(m.CoarsePositions().positions) for m in movers], self._maxCliqueWork)
        if work is None:
          reason = 'optimized as singletons, search work count above the limit of {:.3g}'.format(
            self._maxCliqueWork)
        elif work == _CLIQUE_WORK_ABANDONED:
          reason = ('optimized as singletons, search work not counted: the cut searches alone '
            'visit more than {:.3g} Movers'.format(_CLIQUE_WORK_MAX_STEPS))
        if reason is not None:
          remainingSeconds = 0.0

      # the count is not charged to the allowance
      startTime = time.time()
      (bestScore, infoString) = optC.OptimizeCliqueCoarse(movers, edges, remainingSeconds)
      # as the C++ wrote it: _VerboseCheck() would indent the next line when it is empty
      if self._verbosity >= 1:
        self._infoString += infoString
      ret = bestScore

      stopTime = time.time()
      cutShort = reason is not None
      if reason is None:
        # Brute force has no time check and can finish late: only the C++ line
        # shows that part of the search fell back to singletons.
        if _CLIQUE_SINGLETONS_LINE.search(infoString):
          reason = 'search ran out of {}, part or all of it optimized as singletons'.format(
            'its {:.3g} s allowance'.format(remainingSeconds) if perClique
            else 'the total time limit')
          cutShort = True

      self._infoString += _VerboseCheck(self._verbosity, 1,
        "Clique of size {:d} optimized with score {:.2f} in {:.2f} seconds\n".format(
          len(movers), ret, stopTime - startTime))
      if reason is not None:
        msg = "{}Clique of {:d} Movers (model {}, {}first {}): {}{}\n".format(_CLIQUE_WARNING,
          len(movers), self._cliqueModel + 1,
          'alt {}, '.format(self._cliqueAlt) if self._cliqueAlt.strip() else '',
          self._MoverName(movers[0]).strip(), reason,
          '; its Movers may not be in their best joint orientations' if cutShort else '')
        if msg not in self._warningString:
          self._warningString += msg
    self._infoString += _ReportTiming(self._verbosity, "optimize cliques (coarse)")

  def _OptimizeFine(self, optC, cliques, stopOptimizingAt):
    """Fine optimization of the Movers of the Cliques, then the polish of each
    multi-Mover Clique, skipped after the time limit."""
    # Do fine optimization on the Movers, one at a time.  A move of a member of
    # a multi-Mover Clique must also raise the scores of its neighbours.
    self._infoString += _VerboseCheck(self._verbosity, 1,"Fine optimization on all Movers\n")
    inCliques = set(self._interactionGraph.vertex_label(v) for c in cliques for v in c)
    neighbours = dict((m, []) for m in inCliques)
    for e in self._interactionGraph.edges():
      source = self._interactionGraph.vertex_label(self._interactionGraph.source(e))
      target = self._interactionGraph.vertex_label(self._interactionGraph.target(e))
      if source in neighbours and target in neighbours:
        neighbours[source].append(target)
        neighbours[target].append(source)
    for m in self._movers:
      if m not in inCliques:
        continue
      (bestScore, infoString) = optC.OptimizeSingleMoverFine(m, neighbours[m])
      self._infoString += infoString
    self._infoString += _ReportTiming(self._verbosity, "optimize all Movers (fine)")

    # Polish each multi-Mover Clique: every coarse state of each member with its
    # own fine search, the other members fixed.
    skipped = 0
    stopped = 0
    for c in cliques:
      if len(c) < 2:
        continue
      movers = [self._interactionGraph.vertex_label(v) for v in c]
      (sweeps, moves, capped, stop, infoString) = optC.PolishClique(movers,
        [neighbours[m] for m in movers], stopOptimizingAt - time.time())
      self._infoString += infoString
      if stop:
        if sweeps == 0:
          skipped += 1
        else:
          stopped += 1
      if capped:
        self._infoString += _VerboseCheck(self._verbosity, 1,
          "Polish of a Clique of size {:d} stopped at the limit of {:d} sweeps\n".format(
            len(movers), sweeps))
    if skipped + stopped > 0:
      self._infoString += _VerboseCheck(self._verbosity, 1,
        "Polish skipped for {:d} and stopped for {:d} Cliques due to time limit\n".format(
          skipped, stopped))
    self._infoString += _ReportTiming(self._verbosity, "polish Cliques")

  def getInfo(self):
    """
      Returns information that the user may care about regarding the processing.  The level of
      detail on this can be set by setting the verbosity parameter during Optimizer construction.
      :return: the information so far collected in the string.  Calling this method also clears
      the information, so that later calls will not repeat it.
      If the object was constructed with fillAtomDump = False, then the atom dump will always be
      empty.
    """
    ret = str(self._infoString)
    self._infoString = _InfoAccumulator()
    return ret

  def getWarnings(self):
    """
      Returns warnings that the user may care about regarding the processing.
      :return: the information so far collected in the string.  Calling this method also clears
      the information, so that later calls will not repeat it.
      If the object was constructed with fillAtomDump = False, then the atom dump will always be
      empty.
    """
    ret = self._warningString
    self._warningString = ""
    return ret

  def getAtomDump(self):
    """
      Returns information about the final status of each atom in each model, including VdW radius
      and various flags.  Useful for debugging and for regression testing.
      :return: the information so far collected in the string.  Calling this method also clears
      the information, so that later calls will not repeat it.
    """
    ret = self._atomDump
    self._atomDump = ""
    return ret

  def getHydrogensToDelete(self):
    """
      Returns a list of the Hydrogens that were marked for deletion during the processing.
      The caller should delete these form the optimized structure to produce a model that
      has the best contacts and fewest clashes.
      Each atom's deletion is decided by the last alternate pass that contained it as its
      own (not as environment, EnvironmentResidues()).
      :return: a list of the Hydrogens that were marked for deletion during the processing,
      in model order (i_seq), so that the deletions are reported in the same order on every run.
    """
    return sorted(self._hydrogensToDelete, key=lambda a: a.i_seq)

  def getFlippedAmides(self):
    """
      Returns a list of the Movers that are amide flips.  The caller may want to use this to report
      on the flips that were done.  This includes both "locked down" Movers and those
      that were chosen to be flipped based on their scores.
      :return: a list with a FlippedMoverInfo from each Mover that is an amide flip.
    """
    return list(self._amideFlipsPerformed)

  def getFlippedHistidines(self):
    """
      Returns a list of the Movers that are histidine flips.  The caller may want to use this to report
      on the flips that were done.  This includes both "locked down" Movers and those
      that were chosen to be flipped based on their scores.
      :return: a list with a FlippedMoverInfo from each Mover that is a histidine flip.
    """
    return list(self._hisFlipsPerformed)

  def getSuggestedFlips(self):
    """
      With chooseHisProtonation and addFlipMovers False: the Asn/Gln amide and His ring
      flips that were scored but not made, each alone with everything else as placed,
      whose best flipped state scores better than the pose kept by more than
      preferenceMagnitude * nonFlipPreference, and the His whose ion binds only the
      flipped ring. A later alternate pass that handles a residue replaces or withdraws
      the entry of an earlier one.
      :return: list of group_args (model, alt, atom (ND2 or NE2), kind 'AmideFlip' or
      'HisFlip', resNameAndID, gain (best flipped state minus the kept pose, without
      the flip penalty; None for a His listed for its ion whose flip could not be
      scored), threshold, chosen and flipped (pose descriptions; flipped None when not
      scored), ion (None,
      or group_args(atom (of the model, also for an ion in a symmetry copy), rt_mx
      (None for the model's own ion), ring (CE1 or CD2, whose site the flipped ring N
      takes), distance (of that ring atom from the ion as modelled)))), in model order.
    """
    result = [v for v in self._suggestedFlips.values() if v is not None]
    return sorted(result, key=lambda v: (v.model, v.atom.i_seq))

  def _scorePoseAtoms(self, atoms, passNeighborLists, dotSphereCache):
    """Score of these atoms where they are now, as _ScoreMover scores a pose, without
    the preference energy; deleted atoms do not count."""
    score = 0.0
    for b in atoms:
      if b in self._deleteMes:
        continue
      if b.i_seq not in self._excludeDict:
        self._excludeDict[b.i_seq] = mmtbx.probe.Helpers.getAtomsWithinNBondsAtAnyDistance(b,
          passNeighborLists, self._extraAtomInfo, self._bondedNeighborDepth, 3)
      rad = self._extraAtomInfo.getMappingFor(b).vdwRadius
      score += self._dotScorer.score_dots(b, self._minOccupancy, self._spatialQuery,
        rad + self._maximumVDWRadius, self._probeRadius, self._excludeDict[b.i_seq],
        dotSphereCache.get_sphere(rad).dots(), self._probeDensity, False, False).totalScore()
    return score

  def _SuggestFlips(self, mi, alt, passNeighborLists, dotSphereCache):
    """Score the flip candidates of this pass (flips off) in the final environment:
    the pose kept (its sites as optimized) against every flipped coarse state of a
    full flip Mover, then back to the kept pose. Copies of the residue do not follow
    its flipped states."""
    from libtbx import group_args
    threshold = self._preferenceMagnitude * self._nonFlipPreference
    copyOf = {}
    for c, o, rt in self._passCopies:
      copyOf[c] = (o, rt)
    def suggestKey(a):
      rg = a.parent().parent()
      return (mi, rg.parent().id, rg.resseq, rg.icode, a.parent().resname, a.parent().altloc)
    # this pass decides for every residue it handled, candidate or not (locked to an
    # ion, linked, not scored); entries of other alternate residues stay
    earlier = {}
    for a in self._flipHandled:
      key = suggestKey(a)
      if self._suggestedFlips.get(key) is not None:
        earlier[key] = self._suggestedFlips[key]
      self._suggestedFlips[key] = None
    for candAlt, a, place, ion in self._flipCandidates:
      # The amide Mover re-places the NH2 Hydrogens when it is built: keep the sites
      # of the residue's atoms in this pass (its H may be in a conformer when its N
      # is not), out of the spatial query while it is built, and put them back.
      saved = [(b, b.xyz) for b in a.parent().parent().atoms() if b in passNeighborLists]
      inQuery = [b for b, xyz in saved if b not in self._deleteMes]
      def restoreSites():
        for b in inQuery:
          self._spatialQuery.remove(b)
        for b, xyz in saved:
          b.xyz = xyz
        for b in inQuery:
          self._spatialQuery.add(b)
      m = None
      pr = None
      chosen = 0
      kind = 'AmideFlip'
      if place is not None:
        chosen = self._optimizerOf[place].GetCoarseLocation(place)
        kind = 'HisFlip'
      stateSet = False
      try:
        for b in inQuery:
          self._spatialQuery.remove(b)
        try:
          if place is None:
            m = Movers.MoverAmideFlip(a, "CA", passNeighborLists, self._nonFlipPreference)
          else:
            m = Movers.MoverHisFlip(a, passNeighborLists, self._extraAtomInfo,
              self._nonFlipPreference)
        finally:
          for b in inQuery:
            self._spatialQuery.add(b)
        pr = m.CoarsePositions()
        restoreSites()
        kept = self._preferenceMagnitude * pr.preferenceEnergies[chosen] + \
          self._scorePoseAtoms(pr.atoms[:len(pr.positions[chosen])], passNeighborLists,
            dotSphereCache)
        half = len(pr.positions)//2
        best = None
        bestScore = None
        for i in range(half, len(pr.positions)):
          stateSet = True
          self._setMoverState(pr, i)
          score = self._preferenceMagnitude * pr.preferenceEnergies[i] + \
            self._scorePoseAtoms(pr.atoms[:len(pr.positions[i])], passNeighborLists,
              dotSphereCache)
          if best is None or score > bestScore:
            best = i
            bestScore = score
      except Exception as e:
        # with flips off the residue is only not scored; the default verbosity (2)
        # would show a level-1 line in every default run
        self._infoString += _VerboseCheck(self._verbosity, 3,"Could not score the flip of "+
          _ResNameAndID(a)+": "+str(e)+"\n")
        m = None
      finally:
        if stateSet:
          self._setMoverState(pr, chosen)
        restoreSites()
      if m is None and ion is None:
        continue
      key = suggestKey(a)
      ionInfo = None
      if ion is not None:
        ionAtom, ring, dist = ion
        rt = None
        if ionAtom in copyOf:
          ionAtom, rt = copyOf[ionAtom]
        ionInfo = group_args(atom = ionAtom, rt_mx = rt, ring = ring, distance = dist)
      if m is None:
        # the ion evidence does not depend on a score
        self._suggestedFlips[key] = group_args(model = mi + 1, alt = candAlt, atom = a,
          kind = kind, resNameAndID = _ResNameAndID(a), gain = None, threshold = threshold,
          chosen = place.PoseDescription(chosen, None, False), flipped = None, ion = ionInfo)
        self._infoString += _VerboseCheck(self._verbosity, 1,"Suggested flip: "+kind+" at "+
          _ResNameAndID(a)+" (for its ion; the flip could not be scored)\n")
        continue
      gain = bestScore - kept + threshold
      if gain > threshold or ionInfo is not None:
        self._suggestedFlips[key] = group_args(model = mi + 1, alt = candAlt, atom = a,
          kind = kind, resNameAndID = _ResNameAndID(a), gain = gain, threshold = threshold,
          chosen = m.PoseDescription(chosen, None, False),
          flipped = m.PoseDescription(best, None, False), ion = ionInfo)
        self._infoString += _VerboseCheck(self._verbosity, 1,"Suggested flip: "+kind+" at "+
          _ResNameAndID(a)+" (flipped minus kept score {:.2f})\n".format(gain))
    for key, entry in earlier.items():
      if self._suggestedFlips[key] is None:
        self._infoString += _VerboseCheck(self._verbosity, 1,"Suggested flip of "+
          entry.resNameAndID+" withdrawn: not suggested in this pass\n")

  def _setMoverState(self, positionReturn, index):
    """
      Move the atoms to their new positions, updating the spatial query structure
      by removing the old and adding the new location.
    """
    # there may be fewer moved atoms than there are total atoms
    numMoved = len(positionReturn.positions[index])
    for i in range(numMoved):
      a = positionReturn.atoms[i]

      self._spatialQuery.remove(a)
      # Make a slice here so that we get a copy of the location rather than a reference to it
      a.xyz = positionReturn.positions[index][i][:]
      self._spatialQuery.add(a)

    # Update the extra atom information associated with each atom.
    # Note that there may be fewer entries than atoms, but they all correspond.
    for i, e in enumerate(positionReturn.extraInfos[index]):
      self._extraAtomInfo.setMappingFor(positionReturn.atoms[i], e)

    # Manage the deletion status of each atom, including ensuring
    # consistency with the spatial-query structure.
    # Note that there may be fewer entries than atoms, but they all correspond.
    for i, doDelete in enumerate(positionReturn.deleteMes[index]):
      if doDelete:
        self._spatialQuery.remove(positionReturn.atoms[i])
        self._deleteMes.add(positionReturn.atoms[i])
        self._infoString += _VerboseCheck(self._verbosity, 10,"Deleting atom\n")
      else:
        self._spatialQuery.add(positionReturn.atoms[i])
        self._deleteMes.discard(positionReturn.atoms[i])
        self._infoString += _VerboseCheck(self._verbosity, 10,"Ensuring deletable atom is present\n")

  def _doFixup(self, fixUp):
    """
      Move the atoms to their fixup positions, updating their extra atom info
      and deletion status
    """
    myAtoms = fixUp.atoms
    for i, p in enumerate(fixUp.positions):
      self._infoString += _VerboseCheck(self._verbosity, 5,"Moving atom to {}\n".format(p))
      myAtoms[i].xyz = p
    for i, e in enumerate(fixUp.extraInfos):
      self._infoString += _VerboseCheck(self._verbosity, 5,"Atom info to {}\n".format(e))
      self._extraAtomInfo.setMappingFor(myAtoms[i], e)
    for i, d in enumerate(fixUp.deleteMes):
      # Either ensure that it is deleted or ensure that it is not depending on the
      # value of the deletion result.
      self._infoString += _VerboseCheck(self._verbosity, 5,"Atom deleted is {}\n".format(d))
      if d:
        self._deleteMes.add(myAtoms[i])
      else:
        self._deleteMes.discard(myAtoms[i])

  ##################################################################################
  # Placement

  def _PlaceMovers(self, atoms, rotatableHydrogenIDs, bondedNeighborLists, hParameters,
                    addFlipMovers, alt):
    """Produce a list of Movers for atoms in a pdb.hierarchy.conformer that has added Hydrogens.
    :param atoms: flex array of atoms to search.  This must have all Hydrogens needed by the
    Movers present in the structure already.
    :param rotatableHydrogenIDs: List of sequence IDs that include those of single hydrogens
    that are rotatable. These are hydrogens that are the only one bound to their neighbor
    and they are rotatable.
    :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
    structure that the atom from the first parameter interacts with that lists all of the
    bonded atoms.  Can be obtained by calling mmtbx.probe.Helpers.getBondedNeighborLists()
    on the atoms of the first parameter, so that only bonds within this conformer count.
    :param hParameters: List indexed by sequence ID that stores the riding
    coefficients for hydrogens that have associated dihedral angles.  This can be
    obtained by calling model.setup_riding_h_manager() and then model.get_riding_h_manager().
    :param addFlipMovers: Do we add flip Movers along with other types?
    :param alt: The alternate conformation identifier for the conformer being processed.  This is
    used for our logging purposes when we are recording flipped movers.
    :return: List of atoms that should be unconditionally marked for deletion as a result
    of the analysis done during placement.
    """

    # Parse the flipStates parameter to get a list of FlipMoverState objects.
    fs = _ParseFlipStates(self._flipStates)

    # List of Movers
    self._movers = []

    # Flip Movers not made with flips off: (alt, ND2/NE2 atom, His Mover or None, ion or None)
    self._flipCandidates = []
    # With flips off, the ND2/NE2 of every Asn, Gln and His of this pass, candidate or not
    self._flipHandled = []

    # List of atoms to delete
    deleteAtoms = []

    # i_seq of the His ring Nitrogens that a coarse state of their His flip Mover makes an
    # acceptor; candidate acceptors for the single-hydrogen rotators below.
    hisStateAcceptors = set()

    # Dictionary mapping from Mover to information about the Mover, so that we can keep track of
    # everything from its placement to its final state in the same string.
    self._moverInfo = {}

    # Initial score of each Mover; the blank chain ID drops a token from its _moverInfo text.
    self._initialScore = {}

    # The radius of a Polar Hydrogen (one that forms Hydrogen bonds)
    # @todo Can we look this up from inside CCTBX?
    polarHydrogenRadius = 1.05

    # For each atom, check to see if it is the main atom from a Mover.  If so, attempt to
    # construct the appropriate Mover.  The construction may fail because not all required
    # atoms are present (especially the Hydrogens).  If the construction succeeds, add the
    # Mover to the list of those to return.  When they fail, add the output to the information
    # string. We check single-hydrogen rotators later.

    for a in atoms:
      # Find the stripped upper-case atom and residue names and the residue ID,
      # and an identifier string
      aName = a.name.strip().upper()
      resName = a.parent().resname.strip().upper()
      resNameAndID = _ResNameAndID(a)
      if self._chooseHis and not addFlipMovers and ((aName == 'ND2' and resName == 'ASN') or
          (aName == 'NE2' and resName in ('GLN', 'HIS'))):
        self._flipHandled.append(a)

      # See if we should construct a MoverNH3Rotator here.
      # Find any Nitrogen that has four total bonded neighbors, three of which are Hydrogens.
      # @todo Is this a valid way to search for them?
      if a.element == 'N' and len(bondedNeighborLists[a]) == 4:
        numH = 0
        for n in bondedNeighborLists[a]:
          if n.element_is_hydrogen():
            numH += 1
        if numH == 3:
          try:
            self._movers.append(Movers.MoverNH3Rotator(a, bondedNeighborLists, hParameters))
            self._infoString += _VerboseCheck(self._verbosity, 1,"Added MoverNH3Rotator "+str(len(self._movers))+" to "+resNameAndID+"\n")
            self._moverInfo[self._movers[-1]] = "NH3Rotator at "+resNameAndID+" "+aName;
          except Exception as e:
            self._infoString += _VerboseCheck(self._verbosity, 0,"Could not add MoverNH3Rotator to "+resNameAndID+": "+str(e)+"\n")

      # See if we should construct a MoverAromaticMethylRotator or MoverTetrahedralMethylRotator here.
      # Find any Carbon that has four total bonded neighbors, three of which are Hydrogens.
      # @todo Is this a valid way to search for them?
      if a.element == 'C' and len(bondedNeighborLists[a]) == 4:
        numH = 0
        neighbor = None
        for n in bondedNeighborLists[a]:
          if n.element_is_hydrogen():
            numH += 1
          else:
            neighbor = n
        if numH == 3:
          # See if the Carbon's other neighbor is attached to two other atoms (3 total).  If so,
          # then insert a MoverAromaticMethylRotator and if not, generate a MoverTetrahedralMethylRotator
          # so that the Hydrogens will be staggered but do not add it to those to be optimized.
          if len(bondedNeighborLists[neighbor]) == 3:
            try:
              self._movers.append(Movers.MoverAromaticMethylRotator(a, bondedNeighborLists, hParameters))
              self._infoString += _VerboseCheck(self._verbosity, 1,"Added MoverAromaticMethylRotator "+str(len(self._movers))+" to "+resNameAndID+" "+aName+"\n")
              self._moverInfo[self._movers[-1]] = "AromaticMethylRotator at "+resNameAndID+" "+aName;
            except Exception as e:
              self._infoString += _VerboseCheck(self._verbosity, 0,"Could not add MoverAromaticMethylRotator to "+resNameAndID+" "+aName+": "+str(e)+"\n")
          else:
            try:
              ignored = Movers.MoverTetrahedralMethylRotator(a, bondedNeighborLists, hParameters)
              self._infoString += _VerboseCheck(self._verbosity, 1,"Used MoverTetrahedralMethylRotator to stagger "+resNameAndID+" "+aName+"\n")
            except Exception as e:
              self._infoString += _VerboseCheck(self._verbosity, 0,"Could not add MoverTetrahedralMethylRotator to "+resNameAndID+" "+aName+": "+str(e)+"\n")

      # See if we should insert a MoverAmideFlip here.
      # @todo Is there a more general way than looking for specific names?
      if addFlipMovers and ((aName == 'XD2' and resName == 'ASX') or (aName == 'XE2' and resName == 'GLX')):
        self._infoString += _VerboseCheck(self._verbosity, 1,"Not attempting to adjust "+resNameAndID+" "+aName+"\n")
      if (addFlipMovers or self._chooseHis) and (
          (aName == 'ND2' and resName == 'ASN') or (aName == 'NE2' and resName == 'GLN')):
        # Find the Oxygen and see if it is within a range of ideal bonding distances to a positive ion.
        # Do this in two steps; find the Carbon bonded to the Nitrogen and then the Oxygen bonded to the
        # Carbon.
        foundIon = False
        oxygen = None
        for b in bondedNeighborLists[a]:
          if b.element.upper() == 'C':
            for b2 in bondedNeighborLists[b]:
              if b2.element.upper() == 'O':
                oxygen = b2
        # If we have a close-enough ion, we skip adding a Mover.
        if oxygen is not None:
          # @todo Check both flips and lock down the appropriate one rather than only checking the first.
          myRad = self._extraAtomInfo.getMappingFor(a).vdwRadius
          minDist = myRad
          maxDist = 0.25 + myRad + self._maximumVDWRadius
          neighbors = self._spatialQuery.neighbors(oxygen.xyz, minDist, maxDist)
          for n in neighbors:
            if n.element_is_positive_ion():
              dist = (Helpers.rvec3(oxygen.xyz) - Helpers.rvec3(n.xyz)).length()
              expected = myRad + self._extraAtomInfo.getMappingFor(n).vdwRadius
              self._infoString += _VerboseCheck(self._verbosity, 5,'Checking AmideFlip against '+n.name.strip()+' at '+str(n.xyz)+' from '+str(oxygen.xyz)+
                ' dist = '+str(dist)+', expected = '+str(expected)+'; N rad = '+str(myRad)+
                ', '+n.name.strip()+' rad = '+str(self._extraAtomInfo.getMappingFor(n).vdwRadius)+'\n')
              # @todo Why are we using -0.65 here and -0.55 for Histidine?
              if dist >= (expected - 0.65) and dist <= (expected + 0.25):
                foundIon = True

        if not foundIon and not addFlipMovers:
          # flips off: scored after the optimization, never made
          self._flipCandidates.append((alt, a, None, None))
        elif not foundIon:
          try:
            # Check to see if the state of this Mover has been specified. If so, place it in
            # the requested state and don't insert the Mover. If not, then insert the Mover.
            s = _FindFlipState(a, fs)
            if s is not None:
              self._infoString += _VerboseCheck(self._verbosity, 1,"Setting MoverAmideFlip for "+resNameAndID+": flipped = "+
                str(s.flipped)+", angles adjusted = "+str(s.fixedUp)+"\n")
              flip = Movers.MoverAmideFlip(a, "CA", bondedNeighborLists, self._nonFlipPreference)
              cp = flip.CoarsePositions()
              if s.flipped:
                index = 1
                # Add to the list of flipped amide movers.
                # We're not appending it to the list of Movers, so it won't be counted later.
                self._amideFlipsPerformed.append(FlippedMoverInfo(alt, a))
              else:
                index = 0
              self._setMoverState(cp, index)
              if s.fixedUp:
                self._doFixup(flip.FixUp(index))
            else:
              self._movers.append(Movers.MoverAmideFlip(a, "CA", bondedNeighborLists, self._nonFlipPreference))
              self._infoString += _VerboseCheck(self._verbosity, 1,"Added MoverAmideFlip "+str(len(self._movers))+" to "+resNameAndID+"\n")
              self._moverInfo[self._movers[-1]] = "AmideFlip at "+resNameAndID+" "+aName
          except Exception as e:
            self._infoString += _VerboseCheck(self._verbosity, 0,"Did not add MoverAmideFlip to "+resNameAndID+": "+str(e)+"\n")

      # A His whose protonation the user set keeps its H as placed: no Mover, no ion
      # lock-down (whose H removal would change the state).
      if (aName == 'NE2' and resName == 'HIS' and self._fixedHis and
          _ProtonationKey(a) in self._fixedHis):
        self._infoString += _VerboseCheck(self._verbosity, 1,"Did not add MoverHisFlip to "+resNameAndID+
          " (protonation set by the user; no ion lock-down)\n")
        self._hisNoMover[(alt, a.i_seq)] = HisCharges.USER_PROTONATION

      # See if we should insert a MoverHisFlip here.
      # @todo Is there a more general way than looking for specific names?
      elif aName == 'NE2' and resName == 'HIS':
        hisNE2 = a    # a is reused below
        try:
          # Get a potential Mover and test both of its Nitrogens in the original and flipped
          # locations.  If one or both of them are near enough to be ionically bound to an
          # ion, then we remove the Hydrogen(s) and lock the Histidine at that orientation
          # rather than inserting the Mover into the list of those to be optimized.
          # @todo Consider checking both configurations to see if either one has two bonds.
          hist = Movers.MoverHisFlip(a, bondedNeighborLists, self._extraAtomInfo, self._nonFlipPreference)

          # Find the four positions to check for Nitrogen ionic bonds
          # The two atoms are NE2 (0th atom with its Hydrogen at atom 1) and
          # ND1 (4th atom with its Hydrogen at atom 5).
          # The first flipped configuration follows the unflipped ones.
          cp = hist.CoarsePositions()
          flippedConfig = len(cp.positions)//2
          ne2Orig = cp.positions[0][0]
          nd1Orig = cp.positions[0][4]
          ne2Flip = cp.positions[flippedConfig][0]
          nd1Flip = cp.positions[flippedConfig][4]

          # See if any are close enough to a positive ion to be ionically bonded.
          # If any are, record whether it is the original or
          # the flipped configuration.  Check the original configuration first.
          # Check out to the furthest distance of any atom's VdW radius.
          myRad = self._extraAtomInfo.getMappingFor(a).vdwRadius
          minDist = myRad
          maxDist = 0.25 + myRad + self._maximumVDWRadius
          bondedConfig = None
          bondedIon = None
          for i,pos in enumerate([ne2Orig, nd1Orig, ne2Flip, nd1Flip]):
            neighbors = self._spatialQuery.neighbors(pos, minDist, maxDist)
            for n in neighbors:
              if n.element_is_positive_ion():
                dist = (Helpers.rvec3(pos) - Helpers.rvec3(n.xyz)).length()
                expected = myRad + self._extraAtomInfo.getMappingFor(n).vdwRadius
                self._infoString += _VerboseCheck(self._verbosity, 5,'Checking Histidine '+str(i)+' against '+n.name.strip()+' at '+str(n.xyz)+' from '+str(pos)+
                  ' dist = '+str(dist)+', expected = '+str(expected)+'; N rad = '+str(myRad)+
                  ', '+n.name.strip()+' rad = '+str(self._extraAtomInfo.getMappingFor(n).vdwRadius)+'\n')
                if dist >= (expected - 0.55) and dist <= (expected + 0.25):
                  # The first two elements in the enumeration come from Histidine configuration 0
                  # and the second two from configuration flippedConfig, so we figure out which one we
                  # should be in.
                  bondedConfig = (i // 2) * flippedConfig
                  bondedIon = n
                  bondedPos = i
                  bondedDist = dist
                  self._infoString += _VerboseCheck(self._verbosity, 5,'Found ionic bond in coarse configuration '+str(bondedConfig)+'\n')
                  break
            if bondedConfig is not None:
              # We want the first configuration that is found, so we don't flip if we don't need to.
              break

          # His states chosen without flips: no heavy atom moves, so a His whose ion
          # binds only its flipped ring is not locked; it is listed as a suggested
          # flip with the ion. Other callers without flips keep the flip.
          ionOnlyFlipped = None
          if bondedConfig == flippedConfig and self._chooseHis:
            # the flipped N sites are the modelled CE1 (NE2) and CD2 (ND1) sites
            ionOnlyFlipped = (bondedIon, ('CE1', 'CD2')[bondedPos - 2], bondedDist)
            bondedConfig = None
            self._infoString += _VerboseCheck(self._verbosity, 1,"Not locking "+resNameAndID+
              ": "+bondedIon.name.strip()+" binds only the flipped ring (flips off)\n")

          # If one of the bonded configurations has at least one Ionic bond, then check each of
          # the Nitrogens in that configuration, removing its Hydrogen if it is bonded to an ion.
          # Set the histidine to that flip state; it will not be inserted as a Mover.
          if bondedConfig is not None:
            # Set the histidine in that flip state
            fixUp = hist.FixUp(bondedConfig)
            coarsePositions = hist.CoarsePositions().positions[bondedConfig]
            for i,a in enumerate(fixUp.atoms):
              if i < len(fixUp.positions):
                a.xyz = fixUp.positions[i]
              if i < len(fixUp.extraInfos):
                self._extraAtomInfo.setMappingFor(a, fixUp.extraInfos[i])

            # If the bonded configuration is the flipped one, then we add it to the list of
            # flipped Histidine Movers so that we can report on it later.
            if bondedConfig == flippedConfig:
              self._hisFlipsPerformed.append(FlippedMoverInfo(alt, a))

            # See if we should remove the Hydrogen from each of the two potentially-bonded
            # Nitrogens and make each an acceptor if we do remove its Hydrogen.  The two atoms
            # are NE2 (0th atom with its Hydrogen at atom 1) and ND1 (4th atom with
            # its Hydrogen at atom 5).
            # NOTE that we must check the position of the atom in its coarse configuration rather
            # than its FixUp configuration because that is the location that we use to determine if
            # we have an ionic bond.
            def _modifyIfNeeded(nitro, coarseNitroPos, hydro):
              # Helper function to check and change things for one of the Nitrogens.
              result = ""
              myRad = self._extraAtomInfo.getMappingFor(nitro).vdwRadius
              minDist = myRad
              maxDist = 0.25 + myRad + self._maximumVDWRadius
              neighbors = self._spatialQuery.neighbors(coarseNitroPos, minDist, maxDist)
              for n in neighbors:
                if n.element_is_positive_ion():
                  dist = (Helpers.rvec3(coarseNitroPos) - Helpers.rvec3(n.xyz)).length()
                  expected = myRad + self._extraAtomInfo.getMappingFor(n).vdwRadius
                  if dist >= (expected - 0.55) and dist <= (expected + 0.25):
                    result += _VerboseCheck(self._verbosity, 1,'Not adding Hydrogen to '+resNameAndID+nitro.name+' and marking as an acceptor '+
                      '(ionic bond to '+n.name.strip()+')\n')
                    extra = self._extraAtomInfo.getMappingFor(nitro)
                    extra.isAcceptor = True
                    self._extraAtomInfo.setMappingFor(nitro, extra)
                    deleteAtoms.append(hydro)
                    break
              return result

            self._infoString += _modifyIfNeeded(fixUp.atoms[0], coarsePositions[0], fixUp.atoms[1])
            self._infoString += _modifyIfNeeded(fixUp.atoms[4], coarsePositions[4], fixUp.atoms[5])

            self._infoString += _VerboseCheck(self._verbosity, 1,"Set MoverHisFlip on "+resNameAndID+" to state "+str(bondedConfig)+"\n")
            # 'a' was reused by the loop above
            self._hisNoMover[(alt, hist.CoarsePositions().atoms[0].i_seq)] = 'ion lock (' + bondedIon.name.strip() + ')'
          elif (addFlipMovers or self._chooseHis) and _HisRingNitrogensWithNonstandardBond(
              a, bondedNeighborLists):
            # Issue #1199: a ring Nitrogen carries a bond outside the standard ring
            # connectivity that the ion lock-down above (bondedConfig) did not handle
            # -- a user-defined/custom bond, or an ion the spatial search missed.
            # Flipping the ring could break that bond, so do not add a flip Mover.
            # Hydrogen placement puts a Hydrogen on both ring Nitrogens, so remove
            # the one on each bonded Nitrogen and mark that Nitrogen as an acceptor
            # (as the ion lock-down does in _modifyIfNeeded); otherwise the Hydrogen
            # is left dangling on an atom that is already bonded to its partner.
            self._infoString += _VerboseCheck(self._verbosity, 1,"Did not add MoverHisFlip to "+resNameAndID+
              " (bond outside standard ring connectivity)\n")
            self._hisNoMover[(alt, a.i_seq)] = 'bond outside the ring'
            if ionOnlyFlipped is not None:
              self._infoString += _VerboseCheck(self._verbosity, 1,"Not listing "+resNameAndID+
                ": a ring Nitrogen has a bond outside the ring\n")
            for nitro, hydro in _HisRingNitrogensWithNonstandardBond(a, bondedNeighborLists):
              extra = self._extraAtomInfo.getMappingFor(nitro)
              extra.isAcceptor = True
              self._extraAtomInfo.setMappingFor(nitro, extra)
              if hydro is not None:
                deleteAtoms.append(hydro)
                self._infoString += _VerboseCheck(self._verbosity, 1,"Not adding Hydrogen to "+resNameAndID+
                  nitro.name+" and marking as an acceptor (bond outside standard ring connectivity)\n")
          elif addFlipMovers: # Add a Histidine flip Mover if we're adding flip Movers
            # Check to see if the state of this Mover has been specified. If so, place it in
            # the requested state and insert the Mover as a non-flipping Histidine.
            # If not, then insert the Histidine flip Mover.
            s = _FindFlipState(a, fs)
            if s is not None:
              if s.flipped:
                enabledFlips = 2
                # Record that we have flipped this Histidine so that we can report on it later.
                # We're not appending it to the list of Movers, so it won't be counted later.
                self._hisFlipsPerformed.append(FlippedMoverInfo(alt, a))
              else:
                enabledFlips = 1
              hist = Movers.MoverHisFlip(a, bondedNeighborLists, self._extraAtomInfo, self._nonFlipPreference,
                enabledFlips, s.fixedUp)
              self._movers.append(hist)
              hisStateAcceptors.update(_HisAcceptorsInAnyState(hist))
              self._infoString += _VerboseCheck(self._verbosity, 1,"Added MoverHisPlace "+str(len(self._movers))+" to "+resNameAndID+"\n")
              # Name is HisPlace rather than HisFlip, so we don't try and report uncertain or clash for it.
              self._moverInfo[self._movers[-1]] = "HisPlace at "+resNameAndID+" "+aName;
              # Set the position to the first coarse position, which will cause a flip if required
              coarse = hist.CoarsePositions()
              self._setMoverState(coarse, 0)
            else:
              self._movers.append(hist)
              hisStateAcceptors.update(_HisAcceptorsInAnyState(hist))
              self._infoString += _VerboseCheck(self._verbosity, 1,"Added MoverHisFlip "+str(len(self._movers))+" to "+resNameAndID+"\n")
              self._moverInfo[self._movers[-1]] = "HisFlip at "+resNameAndID+" "+aName;
          elif self._chooseHis:
            # Flips off: the unflipped states only (HD1, HE2 or both), no fix-up.
            place = Movers.MoverHisFlip(a, bondedNeighborLists, self._extraAtomInfo,
              self._nonFlipPreference, 1)
            self._movers.append(place)
            hisStateAcceptors.update(_HisAcceptorsInAnyState(place))
            self._infoString += _VerboseCheck(self._verbosity, 1,"Added MoverHisPlace "+str(len(self._movers))+" to "+resNameAndID+"\n")
            self._moverInfo[self._movers[-1]] = "HisPlace at "+resNameAndID+" "+aName;
            self._flipCandidates.append((alt, a, place, ionOnlyFlipped))
        except Exception as e:
          self._infoString += _VerboseCheck(self._verbosity, 0,"Did not add MoverHisFlip to "+resNameAndID+": "+str(e)+"\n")
          if self._reconcileHisCharges:
            why = 'failed: ' + str(e)
            try:
              nd1 = [x for x in hisNE2.parent().atoms() if x.name.strip() == 'ND1'] or \
                [x for x in hisNE2.parent().parent().atoms() if x.name.strip() == 'ND1']
              for n in [hisNE2] + nd1[:1]:
                ions = [x for x in self._spatialQuery.neighbors(n.xyz, 0.0, 3.0) if x.element_is_positive_ion()]
                if ions and not [x for x in bondedNeighborLists.get(n, []) if x.element_is_hydrogen()]:
                  # placement put no H on this N (metal-bound)
                  why = 'metal-bound %s without H (%s)' % (n.name.strip(), ions[0].name.strip())
                  break
            except (MemoryError, Abort):
              raise
            except Exception:
              pass    # a diagnostic only: the Mover's own error stands
            self._hisNoMover[(alt, hisNE2.i_seq)] = why

    # The Mover constructors above moved atoms after the spatial query was built
    # (rotators turn their Hydrogens to a starting orientation; tetrahedral methyls are
    # staggered; ion-bound Histidines are locked down).  SpatialQuery.remove() finds an
    # atom by its current position, so these atoms stay listed in the bin of their old
    # position, and whether a neighbor query finds them depends on the bin layout, which
    # depends on the bounding box of all atoms, near or far.  Rebuild it from the current
    # positions before the searches below and the optimization.  Symmetry copies stay in it.
    self._spatialQuery = Helpers.createSpatialQuery(list(atoms) + self._environmentAtoms,
      self._probePhil)

    # Make a dictionary looked up up by i_seq that returns the relevant atom. We use this to
    # look up the single-hydrogen rotators. We only place atoms that are in our current model
    # index and alternate.
    atomDict = {}
    for a in atoms:
      atomDict[a.i_seq] = a

    # A symmetry copy is a candidate acceptor as its original is now: the copies took
    # their atom info before placement (ion lock) and have no His states of their own.
    originalOf = {}
    for c, o, rt in self._passCopies:
      originalOf[c.i_seq] = o

    # Insert all single-hydrogen rotators, which we have a list of. Do this after we've
    # done all of the other Movers so that any lock-down placements will have been done
    # before we check.
    for i in rotatableHydrogenIDs:
      # Find out the associated atom
      try:
        a = atomDict[i]
      except KeyError:
        # This atom is not in our current model index or alternate, so skip it.
        continue

      aName = a.name.strip().upper()
      resName = a.parent().resname.strip().upper()
      resNameAndID = _ResNameAndID(a)
      try:
        # Skip Hydrogens that are not bound to an atom that only has a single other
        # atom bound to it -- we will handle those in other cases.
        neighbor = bondedNeighborLists[a][0]
        if len(bondedNeighborLists[neighbor]) == 2:

          # Get the list of atoms that are bonded in a chain to the hydrogen. They
          # cannot be either hydrogen-bond acceptors nor touches/clashes with this atom.
          # By bond count only, as the Optimizer's exclusions: the H rotates.
          bonded = Helpers.getAtomsWithinNBondsAtAnyDistance(a, bondedNeighborLists,
                                                self._extraAtomInfo, self._bondedNeighborDepth, 3)

          # Construct a list of nearby atoms that are potential acceptors and potential touches.
          # The potential acceptors are also potential touches.
          potentialAcceptors = []
          potentialTouches = []

          # Get the list of nearby atoms.  The center of the search is the atom that
          # the Hydrogen is bound to and its radius is 4 (these values are pulled from
          # the Reduce C++ code). We skip ones whose radius is inside the neighbor so
          # that we don't count the hydrogen itself.
          maxDist = 4.0
          nearby = self._spatialQuery.neighbors(neighbor.xyz,
                                                self._extraAtomInfo.getMappingFor(neighbor).vdwRadius,
                                                maxDist)

          # Check each nearby atom to see if its distance from the neighbor is within
          # the sum of the hydrogen-bond length of the neighbor atom, the radius of
          # a polar Hydrogen, and the radius of the nearby atom, indicating potential
          # overlap.
          # O-H & N-H bond len == 1.0, S-H bond len == 1.3
          XHbondlen = Movers.NominalXHBondLength(neighbor)
          candidates = []
          for n in nearby:
            # We don't count bonded atoms.
            if not n in bonded:
              # We treat them all as potential touches/clashes.
              potentialTouches.append(n)
              d = (Helpers.rvec3(neighbor.xyz) - Helpers.rvec3(n.xyz)).length()
              if d <= XHbondlen + self._extraAtomInfo.getMappingFor(n).vdwRadius + polarHydrogenRadius:
                candidates.append(n)

          # See if each candidate atom is a potential acceptor or a flip partner from
          # Histidine or NH2 flips (which may be moved into that atom's position during a flip).
          # We check the partner (N's for GLN And ASN, C's for HIS) because if it is the O or
          # N atom, we will already be checking it as an acceptor now.
          # In any case, it is a potential touch.
          for c in candidates:
            cName = c.name.strip().upper()
            resName = c.parent().resname.strip().upper()
            flipPartner = (
              (cName == 'ND2' and resName == 'ASN') or
              (cName == 'NE2' and resName == 'GLN') or
              (cName == 'CE1' and resName == 'HIS') or (cName == 'CD2' and resName == 'HIS') )
            orig = originalOf.get(c.i_seq, c)
            acceptor = (self._extraAtomInfo.getMappingFor(orig).isAcceptor or
              orig.i_seq in hisStateAcceptors)
            if acceptor or flipPartner:
              potentialAcceptors.append(c)

          oldXYZ = a.xyz
          self._movers.append(Movers.MoverSingleHydrogenRotator(a, bondedNeighborLists, self._extraAtomInfo,
                                                                hParameters, potentialAcceptors,
                                                                potentialTouches))
          # The constructor turned the Hydrogen to its starting orientation; move its
          # spatial-query entry along with it.
          newXYZ = a.xyz
          if newXYZ != oldXYZ:
            a.xyz = oldXYZ
            self._spatialQuery.remove(a)
            a.xyz = newXYZ
            self._spatialQuery.add(a)
          self._infoString += _VerboseCheck(self._verbosity, 1,"Added MoverSingleHydrogenRotator "+str(len(self._movers))+" to "+resNameAndID+" "+aName+
            " with "+str(len(potentialAcceptors))+" potential nearby acceptors\n")
          self._moverInfo[self._movers[-1]] = "SingleHydrogenRotator at "+resNameAndID+" "+aName;
      except Exception as e:
        self._infoString += _VerboseCheck(self._verbosity, 0,"Could not add MoverSingleHydrogenRotator to "+resNameAndID+" "+aName+": "+str(e)+"\n")

    return deleteAtoms


  def _DescribeMover(self, m):
    type_str = str(type(m))
    type_str = type_str.split("'")[1]  # Extract the type name between the quotes
    type_str = type_str.split(".")[-1]  # Extract the class name from the module path
    info = _ResNameAndID(m.CoarsePositions().atoms[0])
    return type_str + ' ' + info

  def _DescribeMoverAtom(self, m, a, rad, loc):
    return ( ' {' + self._DescribeMover(m) + " " + a.name.strip().upper() + '}' +
            ' r= '+str(rad)+' {:0.3f} {:0.3f} {:0.3f}\n'.format(loc[0], loc[1], loc[2]))

  def _MoverSpheres(self, m, extraRadius=0.0):
    """Returns a Kinemage string listing all of the positions of the movable atoms
    in the specified mover over all possible locations. The radius of each is as
    specified in the atom plus the specified extraRadius.
    """
    ret = ''
    for i, cp in enumerate(m.CoarsePositions().positions):
      for j, loc in enumerate(cp):
        a = m.CoarsePositions().atoms[j]
        rad = str(extraRadius + self._extraAtomInfo.getMappingFor(a).vdwRadius)
        ret += self._DescribeMoverAtom(m, a, rad, loc)
      for fp in m.FinePositions(i).positions:
        for j, loc in enumerate(fp):
          a = m.CoarsePositions().atoms[j]
          rad = str(extraRadius + self._extraAtomInfo.getMappingFor(a).vdwRadius)
          ret += self._DescribeMoverAtom(m, a, rad, loc)
    return ret


  def _InteractionKinemage(self, components):
      # Used to round-robin among a set of differentiable colors
      COLORS = ['gray', 'pink', 'sea', 'sky', 'cyan', 'magenta', 'yellow']
      curColor = 0

      # Write all of the Movers in all of the cliques, with a different color per Mover
      # using the actual atom radius.
      ret = '@kinemage 1\n'
      ret += '@master {movers}\n'
      for c in components:
        movers = [self._interactionGraph.vertex_label(i) for i in c]
        for m in movers:
          color = COLORS[curColor]
          curColor = (curColor + 1) % len(COLORS)
          # Radius will be overridden per atom
          ret += '@spherelist color= '+color+' nobutton master= {movers}\n'
          ret += self._MoverSpheres(m)

      # Write all of the cliques, with a different color per clique
      # using the atom radius plus the probe radius.
      for i, c in enumerate(components):
        cliqueName = 'clique '+str(i)+' size '+str(len(c))
        ret += '@master {'+cliqueName+'}\n'
        color = COLORS[curColor]
        curColor = (curColor + 1) % len(COLORS)
        movers = [self._interactionGraph.vertex_label(i) for i in c]
        for m in movers:
          # Radius will be overridden per atom
          ret += '@spherelist color= '+color+' nobutton master= {'+cliqueName+'}\n'
          ret += self._MoverSpheres(m, self._probeRadius)

      return ret


##################################################################################
# Private helper functions.

def _ParseFlipStates(fs):
  """Produce a list of FlipMoverState objects by parsing a comma-separated flipStates
  string.
  """
  ret = []
  if fs is not None:
    for state in fs.split(','):
      words = state.split()
      if len(words) == 6 or len(words) == 7:
        modelId = int(words[0]) + 1
        altId = words[1].lower()
        if altId == '.':
          altId = ''
        chain = words[2].upper()
        resName = words[3].upper()
        if resName == 'HIS':
          t = 'HisFlip'
        else:
          t = 'AmideFlip'
        resIdWithICode = words[4]
        flipped = words[5] == 'Flipped'
        fixedUp = False
        if len(words) == 7:
          fixedUp = words[6] == 'AnglesAdjusted'
        ret.append(FlipMoverState(t, modelId, altId, chain, resName, resIdWithICode, flipped, fixedUp))
  return ret


def _HisRingNitrogensWithNonstandardBond(ne2Atom, bondedNeighborLists):
  '''Issue #1199: find the Histidine ring Nitrogens (ND1/NE2) that carry a bond
  outside the standard ring connectivity -- i.e. to something other than the two
  ring Carbons and the Nitrogen's single Hydrogen (a user-defined/custom bond, or
  a metal). Such a bond would be broken by flipping the ring; and unlike the ring
  Carbons, the ring Nitrogens are not protected by the MoverHisFlip
  neighbour-count checks, so we detect the extra bond here from the restraint bond
  list.
  :param ne2Atom: NE2 atom of the Histidine.
  :param bondedNeighborLists: dictionary from probe.Helpers.getBondedNeighborLists().
  :return: list of (nitrogen, hydrogen) tuples, one per ring Nitrogen that has a
  nonstandard bond. hydrogen is that Nitrogen's bonded Hydrogen, or None if the
  input model does not have it (we do not assume the Hydrogens are present).
  '''
  # ND1 through CE1, by the bonds of this conformer: it may be in another atom_group
  # than NE2 when the ring is split into conformers
  ring = [ne2Atom]
  for c in bondedNeighborLists[ne2Atom]:
    if c.name.strip() != 'CE1':
      continue
    for n in bondedNeighborLists[c]:
      if n.name.strip() == 'ND1' and n not in ring:
        ring.append(n)
  ret = []
  for atom in ring:
    hydrogen = None
    nonstandard = False
    for neighbor in bondedNeighborLists[atom]:
      if neighbor.element_is_hydrogen():
        hydrogen = neighbor
      elif neighbor.element.strip().upper() != 'C':
        nonstandard = True
    if nonstandard:
      ret.append((atom, hydrogen))
  return ret

def _FindFlipState(a, flipStates):
  '''Is an atom in the list of residues that have flip states? If so, return the associated state
  :param a: atom
  :param flipStates: List of FlipMoverState objects
  :return: FlipMoverState if the residue is in the list, None if not
  '''
  ag = a.parent()
  rg = ag.parent()
  chain = rg.parent()
  for fs in flipStates:
    modelId = chain.parent().id
    # We must offset the index of the model by 1 to get to the 1-based model ID
    if ( (modelId == fs.modelId + 1 or modelId == '') and
          (chain.id == fs.chain) and
          (rg.resseq_as_int() == fs.resId) and
          (ag.altloc.strip() == '' or fs.altId == '.' or ag.altloc.strip().lower() == fs.altId.lower()) and
          (rg.icode.strip() == fs.iCode.strip())
        ):
      return fs
  return None


def _subsetGraph(g, keepLabels):
  """
  Return a new graph that is a subset of g that includes only
  the vertices whose labels are listed, and edges between these.
  :param g: Graph to be subsetted.
  :param keepLabels: Labels on the vertices that should end up in the resulting
  graph.  Must be a subset of the vertices in the graph, but it can include all
  of the vertices, in which case it will make a full copy of the graph.
  :return: Boost graph that is a subset of the original graph.
  """

  # We keep a dictionary from labels to the vertex in the new graph that points to
  # that label so that we can construct edges in the new graph.
  vertexForLabel = {}

  # Construct a subgraph that consists only of vertices to be kept and
  # edges both of whose ends are on these vertices.
  ret = graph.adjacency_list(
        vertex_type = "list",   # List so that deletions do not invalidate iterators and descriptors
        )
  for v in g.vertices():
    label = g.vertex_label(v)
    if label in keepLabels:
      vertexForLabel[label] = ret.add_vertex(label)
  for e in g.edges():
    sourceLabel = g.vertex_label( g.source(e) )
    targetLabel = g.vertex_label( g.target(e) )
    if sourceLabel in keepLabels and targetLabel in keepLabels:
      ret.add_edge(vertex1 = vertexForLabel[sourceLabel], vertex2 = vertexForLabel[targetLabel])

  return ret


##################################################################################
# Test function and associated data and helpers to verify that all functions behave properly.

# Class to pass default Probe parameters as if they were in a probePhil structure
class _philLike:
  def __init__(self):
    self.probe_radius = 0.25
    self.density = 16.0
    self.worse_clash_cutoff = 0.5
    self.clash_cutoff = 0.4
    self.contact_cutoff = 0.25
    self.uncharged_hydrogen_cutoff = 0.6
    self.charged_hydrogen_cutoff = 0.8
    self.bump_weight = 10.0
    self.hydrogen_bond_weight = 4.0
    self.gap_weight = 0.25
    self.allow_weak_hydrogen_bonds = False
    self.ignore_ion_interactions = False
    self.set_polar_hydrogen_radius = True
    self.excluded_bond_chain_length = 4

def _optimizeFragment(pdb_raw, bondedNeighborDepth = 4):
  """Returns an optimizer constructed based on the raw PDB data snippet passed in.
  :param pdb_raw: A string that includes a snippet of a PDB file. Should have no alternates.
  :param bondedNeighborDepth: How many levels of bonded neighbors to include in the
  neighbor lists for each atom. Some tests were done with a value of 3, which differs
  from the newer default of 4.
  :return: Optimizer initialized based on the snippet.
  """
  dm = DataManager(['model'])
  dm.process_model_str("my_data.pdb", pdb_raw)
  model = dm.get_model()

  # Add Hydrogens to the model
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model)
  reduce_add_h_obj.run()
  model = reduce_add_h_obj.get_model()

  # Interpret the model after adding Hydrogens to it so that
  # all of the needed fields are filled in when we use them below.
  # @todo Remove this once place_hydrogens() does all the interpretation we need.
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.allow_polymer_cross_special_position=True
  p.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
  p.pdb_interpretation.proceed_with_excessive_length_bonds=True
  p.pdb_interpretation.disable_uc_volume_vs_n_atoms_check=True
  model.process(make_restraints=True,pdb_interpretation_params=p) # make restraints

  # Optimization will place the movers.
  probePhil = _philLike()
  return Optimizer(probePhil, True, model, bondedNeighborDepth=bondedNeighborDepth)

def _numMovers(opt):
  """The number of Movers of the last alternate pass, from the info string (the
  Optimizer releases its Movers at the end of its construction); verbosity >= 1."""
  return int(re.findall(r'Inserted (\d+) Movers', str(opt._infoString))[-1])

def Test(inFileName = None, dumpAtoms = False):

  """Test function for all functions provided above.
  :param inFileName: Name of a PDB or CIF file to load (default makes a small molecule)
  :return: Empty string on success, string describing the problem on failure.
  """

  #========================================================================
  # Test the imported C++ functions
  ret = Optimizers_test()
  if len(ret) > 0:
    return "Optimizers.Test(): " + ret

  #========================================================================
  # Test model to use for validating the alternate/conformer-selection functions.

  alternates_test = (
"""
ATOM      1  N   HIS A  1       26.965  32.911   7.593  1.00  7.19           N
ATOM      2  N  AHIS A  2       26.965  32.911   7.593  1.00  7.19           N
ATOM      3  N  BHIS A  2       26.965  32.911   7.593  1.00  7.19           N
ATOM      4  N  CHIS A  2       26.965  32.911   7.593  1.00  7.19           N
END
"""
    )

  #========================================================================
  # Test the AlternatesInModel() function.
  # Test the GetAtomsForConformer() function.
  dm = DataManager(['model'])
  dm.process_model_str("alternates_test.pdb",alternates_test)
  model = dm.get_model()
  model.process(make_restraints=False)
  model = dm.get_model().get_hierarchy().models()[0]
  alts = AlternatesInModel(model)
  if alts != set(['','A','B','C']):
    return "Optimizers.Test(): Incorrect set of alternates for AlternatesInModel() test: " + str(alts)

  for a in alts:
    count = len(GetAtomsForConformer(model, a))
    if count != 2:
      return "Optimizers.Test(): Incorrect atom count for GetAtomsForConformer() test: " + str(count)

  #========================================================================
  # Unit tests for each type of Optimizer.

  ################################################################################
  # Test using a snippet from 7c31 to make sure the single-hydrogen rotator sets
  # the orientation of the Hydrogen to make good contact with a nearby Oxygen. The B
  # alternate has been removed and the A moved to the main alternate. This is not
  # close enough for a good Hydrogen bond, so it does not point straight at the Oxygen.

  pdb_7c31_two_residues = """\
CRYST1   27.854   27.854   99.605  90.00  90.00  90.00 P 43          8
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.035901  0.000000  0.000000        0.00000
SCALE2      0.000000  0.035901  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010040        0.00000
ATOM     68  N   SER A   5     -31.155  49.742   0.887  1.00 10.02           N
ATOM     69  CA  SER A   5     -32.274  48.937   0.401  1.00  9.76           C
ATOM     70  C   SER A   5     -33.140  49.851  -0.454  1.00  9.47           C
ATOM     71  O   SER A   5     -33.502  50.939  -0.012  1.00 10.76           O
ATOM     72  CB  SER A   5     -33.086  48.441   1.599  1.00 12.34           C
ATOM     73  OG  SER A   5     -34.118  47.569   1.179  1.00 19.50           O
ATOM    758  N   VAL B   2     -34.430  42.959   3.043  1.00 19.95           N
ATOM    759  CA  VAL B   2     -32.994  42.754   2.904  0.50 18.09           C
ATOM    760  C   VAL B   2     -32.311  44.083   2.629  1.00 17.65           C
ATOM    761  O   VAL B   2     -32.869  44.976   1.975  1.00 18.97           O
ATOM    762  CB  VAL B   2     -32.703  41.738   1.778  0.50 19.70           C
ATOM    763  CG1 VAL B   2     -33.098  42.307   0.419  0.50 19.61           C
ATOM    764  CG2 VAL B   2     -31.224  41.334   1.755  0.50 16.49           C
TER    1447      CYS B  47
END
"""
  opt = _optimizeFragment(pdb_7c31_two_residues)
  numMovers = _numMovers(opt)
  if numMovers != 3:
    return "Optimizers.Test(): Incorrect number of Movers for single-hydrogen rotator test: " + str(numMovers)

  # See what the pose angle is on the Mover. It should be 62 degrees, and is reported
  # after 'pose Angle '.
  # NOTE: This test was generated when the SingleHydrogenRotator was testing all coarse
  # angles in addition to those towards acceptors. When we switched to the original Reduce
  # behavior of only doing a single non-clashing angle along with the acceptor angles, this
  # test no longer passes. Leaving it in for now in case we want to switch back to the
  # original behavior.
  #angle = int(re.search(r'(?<=pose Angle )[-+]?\d+', opt.getInfo()).group(0))
  #if angle != 62:
  #  return "Optimizers.Test(): Unexpected angle ("+str(angle)+") for single-hydrogen rotator, expected 62"

  ################################################################################
  # Test using a modified snippet from 7c31 to make sure the single-hydrogen rotator sets
  # the orientation of the Hydrogen to make a good hydrogen bond with a nearby Oxygen.
  # The above example has been stripped down and has the Oxygen moved closer.

  pdb_7c31_two_residues_close_O = """\
CRYST1   27.854   27.854   99.605  90.00  90.00  90.00 P 43          8
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.035901  0.000000  0.000000        0.00000
SCALE2      0.000000  0.035901  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010040        0.00000
ATOM     68  N   SER A   5     -31.155  49.742   0.887  1.00 10.02           N
ATOM     69  CA  SER A   5     -32.274  48.937   0.401  1.00  9.76           C
ATOM     70  C   SER A   5     -33.140  49.851  -0.454  1.00  9.47           C
ATOM     71  O   SER A   5     -33.502  50.939  -0.012  1.00 10.76           O
ATOM     72  CB  SER A   5     -33.086  48.441   1.599  1.00 12.34           C
ATOM     73  OG  SER A   5     -34.118  47.569   1.179  1.00 19.50           O
ATOM    761  O   VAL B   2     -33.462  45.185   1.977  1.00 18.97           O
TER    1447      CYS B  47
END
"""
  opt = _optimizeFragment(pdb_7c31_two_residues_close_O)
  numMovers = _numMovers(opt)
  if numMovers != 2:
    return "Optimizers.Test(): Incorrect number of Movers for single-hydrogen rotator H-bond test: " + str(numMovers)

  # See what the pose angle is on the Mover. It is reported after 'pose Angle '.
  expected = 113
  angle = int(re.findall(r'(?<=pose Angle )[-+]?\d+', opt.getInfo())[1])
  if angle != expected:
    return "Optimizers.Test(): Unexpected angle ("+str(angle)+") for single-hydrogen rotator H-bond, expected "+str(expected)+", found "+str(angle)

  ################################################################################
  # Test using a modified snippet from 7c31 to make sure the single-hydrogen rotator sets
  # the orientation of the Hydrogen touching a nearby Carbon.
  # The above example has the Oxygen replaced by a Carbon.

  pdb_7c31_two_residues_close_C = """\
CRYST1   27.854   27.854   99.605  90.00  90.00  90.00 P 43          8
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.035901  0.000000  0.000000        0.00000
SCALE2      0.000000  0.035901  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010040        0.00000
ATOM     68  N   SER A   5     -31.155  49.742   0.887  1.00 10.02           N
ATOM     69  CA  SER A   5     -32.274  48.937   0.401  1.00  9.76           C
ATOM     70  C   SER A   5     -33.140  49.851  -0.454  1.00  9.47           C
ATOM     71  O   SER A   5     -33.502  50.939  -0.012  1.00 10.76           O
ATOM     72  CB  SER A   5     -33.086  48.441   1.599  1.00 12.34           C
ATOM     73  OG  SER A   5     -34.118  47.569   1.179  1.00 19.50           O
ATOM    761  C   VAL B   2     -33.462  45.185   1.977  1.00 18.97           C
TER    1447      CYS B  47
END
"""
  opt = _optimizeFragment(pdb_7c31_two_residues_close_C)
  numMovers = _numMovers(opt)
  if numMovers != 2:
    return "Optimizers.Test(): Incorrect number of Movers for single-hydrogen rotator clash test: " + str(numMovers)

  # See what the pose angle is on the Mover. It is reported after 'pose Angle '.
  # HG touches VAL B 2 C (2.62 A). It was -64, touching its own C (2.67 A, 4 bonds
  # away), which the excluded atoms missed when they depended on the first pose.
  expected = -163
  angle = int(re.findall(r'(?<=pose Angle )[-+]?\d+', opt.getInfo())[1])
  if angle != expected:
    return "Optimizers.Test(): Unexpected angle ("+str(angle)+") for single-hydrogen rotator clash, expected "+str(expected)+", found "+str(angle)

  ################################################################################
  # Test using a modified snippet from 1ubq to make sure the NH3 rotator sets
  # the orientation as expected to align with a nearby Oxygen and waters.

  pdb_1ubq_two_residues_close_O = """\
CRYST1   50.840   42.770   28.950  90.00  90.00  90.00 P 21 21 21    4
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.019670  0.000000  0.000000        0.00000
SCALE2      0.000000  0.023381  0.000000        0.00000
SCALE3      0.000000  0.000000  0.034542        0.00000
ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00  9.67           N
ATOM      2  CA  MET A   1      26.266  25.413   2.842  1.00 10.38           C
ATOM      3  C   MET A   1      26.913  26.639   3.531  1.00  9.62           C
ATOM      4  O   MET A   1      27.886  26.463   4.263  1.00  9.62           O
ATOM      5  CB  MET A   1      25.112  24.880   3.649  1.00 13.77           C
ATOM      6  CG  MET A   1      25.353  24.860   5.134  1.00 16.29           C
ATOM      7  SD  MET A   1      23.930  23.959   5.904  1.00 17.17           S
ATOM      8  CE  MET A   1      24.447  23.984   7.620  1.00 16.11           C
ATOM    127  N   VAL A  17      30.310  25.458   5.384  1.00  8.99           N
ATOM    128  CA  VAL A  17      30.288  24.245   6.193  1.00  8.85           C
ATOM    129  C   VAL A  17      29.279  23.227   5.641  1.00  8.04           C
ATOM    130  O   VAL A  17      28.478  23.522   4.725  1.00  8.99           O
ATOM    131  CB  VAL A  17      29.903  24.590   7.665  1.00  9.78           C
ATOM    132  CG1 VAL A  17      30.862  25.496   8.389  1.00 12.05           C
ATOM    133  CG2 VAL A  17      28.476  25.135   7.705  1.00 10.54           C
TER     603      GLY A  76
HETATM  625  O   HOH A  98      25.928  21.774   2.325  1.00 13.70           O
HETATM  637  O   HOH A 110      28.824  25.094   0.886  0.77 36.99           O
END
"""
  # This test is done with a bonded-neighbor depth of 3 because that is how we were
  # doing the calculations when we set it up.
  opt = _optimizeFragment(pdb_1ubq_two_residues_close_O, 3)
  numMovers = _numMovers(opt)
  if numMovers != 2:
    return "Optimizers.Test(): Incorrect number of Movers for NH3 rotator test: " + str(numMovers)

  # See what the pose angle is on the Mover. It should be 163 degrees, and is reported
  # after 'pose Angle '.
  angle = int(re.findall(r'(?<=pose Angle )[-+]?\d+', opt.getInfo())[0])
  if angle != 163:
    return "Optimizers.Test(): Unexpected angle ("+str(angle)+") for NH3 rotator, expected 163"

  ################################################################################
  # Test using snippet from 1dfu to ensure that the Histidine optimzization code will
  # properly skip hydrogens that are going to be deleted when scoring.
  pdb_1dfu_his80 = (
"""
CRYST1   75.600   76.600   95.100  90.00  90.00  90.00 C 2 2 21      8
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.013257  0.000000  0.000000        0.00000
SCALE2      0.000000  0.013063  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010515        0.00000
ATOM   1443  N   HIS P  80      27.889  23.654 -11.281  1.00 41.74           N
ATOM   1444  CA  HIS P  80      28.738  24.315 -10.299  1.00 45.36           C
ATOM   1445  C   HIS P  80      29.467  25.487 -10.941  1.00 46.80           C
ATOM   1446  O   HIS P  80      29.927  25.397 -12.080  1.00 47.32           O
ATOM   1447  CB  HIS P  80      29.765  23.342  -9.714  1.00 47.26           C
ATOM   1448  CG  HIS P  80      30.457  23.865  -8.492  1.00 48.91           C
ATOM   1449  ND1 HIS P  80      29.947  23.705  -7.221  1.00 49.73           N
ATOM   1450  CD2 HIS P  80      31.593  24.589  -8.351  1.00 49.34           C
ATOM   1451  CE1 HIS P  80      30.739  24.306  -6.351  1.00 49.69           C
ATOM   1452  NE2 HIS P  80      31.745  24.850  -7.011  1.00 50.50           N
ATOM   1453  N   PRO P  81      29.581  26.606 -10.211  1.00 48.50           N
ATOM   1454  CA  PRO P  81      30.255  27.810 -10.698  1.00 49.93           C
ATOM   1455  C   PRO P  81      31.662  27.557 -11.236  1.00 51.19           C
ATOM   1456  O   PRO P  81      32.041  28.115 -12.265  1.00 51.83           O
ATOM   1457  CB  PRO P  81      30.264  28.714  -9.468  1.00 50.10           C
ATOM   1458  CG  PRO P  81      29.003  28.336  -8.774  1.00 49.68           C
ATOM   1459  CD  PRO P  81      29.026  26.831  -8.864  1.00 49.20           C
ATOM   1460  N   TYR P  82      32.432  26.716 -10.549  1.00 51.94           N
ATOM   1461  CA  TYR P  82      33.797  26.440 -10.985  1.00 52.99           C
ATOM   1462  C   TYR P  82      34.257  24.985 -10.876  1.00 52.34           C
ATOM   1463  O   TYR P  82      35.456  24.712 -10.890  1.00 52.96           O
ATOM   1464  CB  TYR P  82      34.769  27.349 -10.225  1.00 54.88           C
ATOM   1465  CG  TYR P  82      34.757  27.162  -8.722  1.00 56.85           C
ATOM   1466  CD1 TYR P  82      35.404  26.078  -8.127  1.00 57.60           C
ATOM   1467  CD2 TYR P  82      34.096  28.070  -7.895  1.00 57.47           C
ATOM   1468  CE1 TYR P  82      35.394  25.903  -6.743  1.00 58.41           C
ATOM   1469  CE2 TYR P  82      34.080  27.904  -6.510  1.00 58.47           C
ATOM   1470  CZ  TYR P  82      34.730  26.819  -5.942  1.00 58.59           C
ATOM   1471  OH  TYR P  82      34.718  26.650  -4.577  1.00 59.33           O
HETATM 1782  O   HOH P 134      28.550  22.420  -4.280  1.00 33.59           O
END
"""
    )
  opt = _optimizeFragment(pdb_1dfu_his80)
  numMovers = _numMovers(opt)
  if numMovers != 3:
    return "Optimizers.Test(): Incorrect number of Movers for His placement test: " + str(numMovers)
  info = opt.getInfo()
  if not 'HD1NotPlaced' in info:
    return "Optimizers.Test(): Unexpected HD1 placement for His."
  if not 'HE2Placed' in info:
    return "Optimizers.Test(): Missing HE2 placement for His."


  ################################################################################
  # Test using snippet from 1xso to ensure that the Histidine placement code will lock down the
  # Histidine, set its Nitrogen states, and mark its Hydrogens for deletion.  To
  # ensure that hydrogen placement puts them there in the first place, we first
  # move the CU and ZN far from the Histidine before adding Hydrogens, then move
  # them back before building the spatial hierarchy and testing.
  pdb_1xso_his_61_and_ions = (
"""
ATOM    442  N   HIS A  61      26.965  32.911   7.593  1.00  7.19           N
ATOM    443  CA  HIS A  61      27.557  32.385   6.403  1.00  7.24           C
ATOM    444  C   HIS A  61      28.929  31.763   6.641  1.00  7.38           C
ATOM    445  O   HIS A  61      29.744  32.217   7.397  1.00  9.97           O
ATOM    446  CB  HIS A  61      27.707  33.547   5.385  1.00  9.38           C
ATOM    447  CG  HIS A  61      26.382  33.956   4.808  1.00  8.78           C
ATOM    448  ND1 HIS A  61      26.168  34.981   3.980  1.00  9.06           N
ATOM    449  CD2 HIS A  61      25.174  33.397   5.004  1.00 11.08           C
ATOM    450  CE1 HIS A  61      24.867  35.060   3.688  1.00 12.84           C
ATOM    451  NE2 HIS A  61      24.251  34.003   4.297  1.00 11.66           N
HETATM 2190 CU    CU A   1      22.291  33.388   3.996  1.00 13.22          CU
HETATM 2191 ZN    ZN A 152      27.539  36.010   2.881  1.00  9.34          ZN
END
"""
    )
  dm = DataManager(['model'])
  dm.process_model_str("1xso_snip.pdb",pdb_1xso_his_61_and_ions)
  model = dm.get_model()

  # Make sure we have a valid unit cell.  Do this before we add hydrogens to the model
  # to make sure we have a valid unit cell.
  model = shift_and_box_model(model = model)

  # Find and move the Copper and Zinc atoms far from the Histidine so that
  # Hydrogen placement and bond proxies won't consider them to be bonded.
  # This tests the algorithm behavior for the case of all Hydrogens present
  # in all cases.
  for a in model.get_hierarchy().models()[0].atoms():
    if a.element.upper() == "CU":
      origPositionCU = a.xyz
      a.xyz = (1000, 1000, 1000)
    if a.element.upper() == "ZN":
      origPositionZN = a.xyz
      a.xyz = (1000, 1000, 1000)

  # Add Hydrogens to the model
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model)
  reduce_add_h_obj.run()
  model = reduce_add_h_obj.get_model()

  # Interpret the model after adding Hydrogens to it so that
  # all of the needed fields are filled in when we use them below.
  # @todo Remove this once place_hydrogens() does all the interpretation we need.
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.allow_polymer_cross_special_position=True
  p.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
  p.pdb_interpretation.proceed_with_excessive_length_bonds=True
  model.process(make_restraints=True,pdb_interpretation_params=p) # make restraints

  # Get the first model in the hierarchy.
  firstModel = model.get_hierarchy().models()[0]

  # Get the list of alternate conformation names present in all chains for this model.
  alts = AlternatesInModel(firstModel)

  # Get the atoms from the first conformer in the first model (the empty string is the name
  # of the first conformation in the model; if there is no empty conformation, then it will
  # pick the first available conformation for each atom group.
  atoms = GetAtomsForConformer(firstModel, "")

  # Get the Cartesian positions of all of the atoms we're considering for this alternate
  # conformation.
  carts = flex.vec3_double()
  for a in atoms:
    carts.append(a.xyz)

  # Get the bond proxies for the atoms in the model and conformation we're using and
  # use them to determine the bonded neighbor lists.
  bondProxies, asuProxies = model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)
  bondedNeighborLists = Helpers.getBondedNeighborLists(atoms, bondProxies, asuProxies)

  # Get the probeExt.ExtraAtomInfo needed to determine which atoms are potential acceptors.
  probePhil = _philLike()
  ret = Helpers.getExtraAtomInfo(model = model, bondedNeighborLists = bondedNeighborLists,
      useNeutronDistances=False, probePhil=probePhil)
  extra = ret.extraAtomInfo

  # Also compute the maximum VDW radius among all atoms.
  maxVDWRad = 1
  for a in atoms:
    maxVDWRad = max(maxVDWRad, extra.getMappingFor(a).vdwRadius)

  # Put the Copper and Zinc back in their original positions before we build the
  # spatial-query structure.  This will make them close enough to be bonded to
  # the Nitrogens and should cause Hydrogen removal and marking of the Nitrogens
  # as acceptors.
  for a in model.get_hierarchy().models()[0].atoms():
    if a.element.upper() == "CU":
      a.xyz = origPositionCU
    if a.element.upper() == "ZN":
      a.xyz = origPositionZN

  # Optimization will place the movers, which should be none because the Histidine flip
  # will be constrained by the ionic bonds. There is a rotatable hydrogen placed at the
  # terminus, but no Flip Mover.
  probePhil = _philLike()
  opt = Optimizer(probePhil, True, model)
  numMovers = _numMovers(opt)
  if numMovers != 1:
    return "Optimizers.Test(): Incorrect number of Movers for 1xso Histidine test: " + str(numMovers)

  # Make sure that the two ring Nitrogens have been marked as acceptors.
  # Make sure that the two hydrogens have been marked for deletion.
  for a in model.get_hierarchy().models()[0].atoms():
    name = a.name.strip()
    if name in ["ND1", "NE2"]:
      if not opt._extraAtomInfo.getMappingFor(a).isAcceptor:
        return 'Optimizers.Test(): '+name+' in 1xso Histidine test was not an acceptor'
    if name in ["HD1", "HE2"]:
      if not a in opt._deleteMes:
        return 'Optimizers.Test(): '+name+' in 1xso Histidine test was not set for deletion'

  #========================================================================
  # Check a clique with multiple elements to be sure that it was properly
  # globally optimized. This is a set of ACT residues that are offset
  # such that they want to line up the same way.  A carbon is placed to
  # force the oriention away from the initial solution.
  pdb_multi_act = (
"""
CRYST1   93.586  127.886  251.681  90.00  90.00  90.00 I 2 2 2
SCALE1      0.010685  0.000000  0.000000        0.00000
SCALE2      0.000000  0.007819  0.000000        0.00000
SCALE3      0.000000  0.000000  0.003973        0.00000
HETATM    1  C   ACT A   1       6.494 -47.273 -37.006  1.00 16.65           C
HETATM    2  O   ACT A   1       5.654 -47.981 -37.645  1.00 15.33           O
HETATM    3  CH3 ACT A   1       6.790 -47.410 -35.548  1.00 17.13           C
HETATM    4  OXT ACT A   1       7.110 -46.425 -37.699  1.00 17.05           O
HETATM  101  C   ACT A 101       8.994 -49.773 -37.606  1.00 16.65           C
HETATM  102  O   ACT A 101       8.154 -50.481 -38.245  1.00 15.33           O
HETATM  103  CH3 ACT A 101       9.290 -49.910 -36.148  1.00 17.13           C
HETATM  104  OXT ACT A 101       9.610 -48.925 -38.299  1.00 17.05           O
HETATM  111  C   ACT A 111      11.494 -52.273 -38.206  1.00 16.65           C
HETATM  112  O   ACT A 111      10.654 -52.981 -38.845  1.00 15.33           O
HETATM  113  CH3 ACT A 111      11.790 -52.410 -36.748  1.00 17.13           C
HETATM  114  OXT ACT A 111      12.110 -51.425 -38.899  1.00 17.05           O
HETATM  121  C   ACT A 121      13.994 -54.773 -38.806  1.00 16.65           C
HETATM  122  O   ACT A 121      13.154 -55.481 -39.445  1.00 15.33           O
HETATM  123  CH3 ACT A 121      14.290 -54.910 -37.348  1.00 17.13           C
HETATM  124  OXT ACT A 121      14.610 -53.925 -39.499  1.00 17.05           O
HETATM  131  C   ACT A 131      16.494 -57.273 -39.406  1.00 16.65           C
HETATM  132  O   ACT A 131      15.654 -57.981 -40.045  1.00 15.33           O
HETATM  133  CH3 ACT A 131      16.790 -57.410 -37.948  1.00 17.13           C
HETATM  134  OXT ACT A 131      17.110 -56.425 -40.099  1.00 17.05           O
HETATM  141  C   ACT A 141      18.994 -59.773 -40.006  1.00 16.65           C
HETATM  142  O   ACT A 141      18.154 -60.481 -40.645  1.00 15.33           O
HETATM  143  CH3 ACT A 141      19.290 -59.910 -38.548  1.00 17.13           C
HETATM  144  OXT ACT A 141      19.610 -58.925 -40.699  1.00 17.05           O
HETATM  200  C   ACT A 200       4.500 -45.500 -35.000  1.00 16.65           C
END
"""
    )
  dm = DataManager(['model'])
  dm.process_model_str("pdb_multi_act.pdb",pdb_multi_act)
  model = dm.get_model()

  # Make sure we have a valid unit cell.  Do this before we add hydrogens to the model
  # to make sure we have a valid unit cell.
  model = shift_and_box_model(model = model)

  # Add Hydrogens to the model
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model)
  reduce_add_h_obj.run()
  model = reduce_add_h_obj.get_model()

  # Interpret the model after adding Hydrogens to it so that
  # all of the needed fields are filled in when we use them below.
  # @todo Remove this once place_hydrogens() does all the interpretation we need.
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.allow_polymer_cross_special_position=True
  p.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
  p.pdb_interpretation.proceed_with_excessive_length_bonds=True
  model.process(make_restraints=True,pdb_interpretation_params=p) # make restraints

  # Get the first model in the hierarchy.
  firstModel = model.get_hierarchy().models()[0]

  # Get the list of alternate conformation names present in all chains for this model.
  alts = AlternatesInModel(firstModel)

  # Get the atoms from the first conformer in the first model (the empty string is the name
  # of the first conformation in the model; if there is no empty conformation, then it will
  # pick the first available conformation for each atom group.
  atoms = GetAtomsForConformer(firstModel, "")

  # Get the Cartesian positions of all of the atoms we're considering for this alternate
  # conformation.
  carts = flex.vec3_double()
  for a in atoms:
    carts.append(a.xyz)

  # Get the bond proxies for the atoms in the model and conformation we're using and
  # use them to determine the bonded neighbor lists.
  bondProxies, asuProxies = model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)
  bondedNeighborLists = Helpers.getBondedNeighborLists(atoms, bondProxies, asuProxies)

  # Get the probeExt.ExtraAtomInfo needed to determine which atoms are potential acceptors.
  probePhil = _philLike()
  ret = Helpers.getExtraAtomInfo(model = model, bondedNeighborLists = bondedNeighborLists,
      useNeutronDistances=False,probePhil=probePhil)
  extra = ret.extraAtomInfo

  # Also compute the maximum VDW radius among all atoms.
  maxVDWRad = 1
  for a in atoms:
    maxVDWRad = max(maxVDWRad, extra.getMappingFor(a).vdwRadius)

  # Optimization will place the movers. Make sure we got as many as we expected.
  # Make sure that the orientation for all of the movers is correct.
  # Test with each type of optimizer, from the base to the more derived, so
  # that we find out about failures on the base classes first.
  probePhil = _philLike()

  global _DoCliqueOptimizationInC

  print('Testing Optimizer')
  opt = Optimizer(probePhil, True, model, modelIndex = None, altID = None,
                bondedNeighborDepth = 4,
                useNeutronDistances = False,
                minOccupancy = 0.02,
                preferenceMagnitude = 1.0,
                nonFlipPreference = 0.5,
                skipBondFixup = False,
                flipStates = '',
                verbosity = 1)
  numMovers = _numMovers(opt)
  if numMovers != 6:
    return "Optimizers.Test(): Incorrect number of Movers for Optimizer multi-ACT test: " + str(numMovers)
  res = re.findall(r'pose Angle [-+]?\d+', opt.getInfo())
  for r in res:
    if not 'pose Angle 90' in r:
      return "Optimizers.Test(): Unexpected angle ("+str(r)+") for Optimizer multi-ACT test"

  #========================================================================
  # Check a case where an AmideFlip would be locked down and have its Hydrogen removed.
  # @todo

  #========================================================================
  # Check that the occupancy and B-factor cut-offs for water Oxygens are causing them
  # to be ignored in the calculations by putting an atom in the way and making sure it
  # is ignored.
  pdb_b_factor = (
"""
CRYST1   93.586  127.886  251.681  90.00  90.00  90.00 I 2 2 2
SCALE1      0.010685  0.000000  0.000000        0.00000
SCALE2      0.000000  0.007819  0.000000        0.00000
SCALE3      0.000000  0.000000  0.003973        0.00000
HETATM    1  C   ACT A   1       6.494 -47.273 -37.006  1.00 16.65           C
HETATM    2  O   ACT A   1       5.654 -47.981 -37.645  1.00 15.33           O
HETATM    3  CH3 ACT A   1       6.790 -47.410 -35.548  1.00 17.13           C
HETATM    4  OXT ACT A   1       7.110 -46.425 -37.699  1.00 17.05           O
HETATM  200  O   HOH A 200       8.200 -48.610 -35.148  1.00 50.00           O
HETATM  201  O   HOH A 201       5.800 -46.810 -34.548  1.00 20.00           O
END
"""
    )
  dm = DataManager(['model'])
  dm.process_model_str("pdb_b_factor.pdb",pdb_b_factor)
  model = dm.get_model()

  # Make sure we have a valid unit cell.  Do this before we add hydrogens to the model
  # to make sure we have a valid unit cell.
  model = shift_and_box_model(model = model)

  # Add Hydrogens to the model
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model)
  reduce_add_h_obj.run()
  model = reduce_add_h_obj.get_model()

  # Interpret the model after adding Hydrogens to it so that
  # all of the needed fields are filled in when we use them below.
  # @todo Remove this once place_hydrogens() does all the interpretation we need.
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.allow_polymer_cross_special_position=True
  p.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
  p.pdb_interpretation.proceed_with_excessive_length_bonds=True
  model.process(make_restraints=True,pdb_interpretation_params=p) # make restraints

  # Get the first model in the hierarchy.
  firstModel = model.get_hierarchy().models()[0]

  # Get the list of alternate conformation names present in all chains for this model.
  alts = AlternatesInModel(firstModel)

  # Get the atoms from the first conformer in the first model (the empty string is the name
  # of the first conformation in the model; if there is no empty conformation, then it will
  # pick the first available conformation for each atom group.
  atoms = GetAtomsForConformer(firstModel, "")

  # Get the Cartesian positions of all of the atoms we're considering for this alternate
  # conformation.
  carts = flex.vec3_double()
  for a in atoms:
    carts.append(a.xyz)

  # Get the bond proxies for the atoms in the model and conformation we're using and
  # use them to determine the bonded neighbor lists.
  bondProxies, asuProxies = model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)
  bondedNeighborLists = Helpers.getBondedNeighborLists(atoms, bondProxies, asuProxies)

  # Get the probeExt.ExtraAtomInfo needed to determine which atoms are potential acceptors.
  probePhil = _philLike()
  ret = Helpers.getExtraAtomInfo(model = model, bondedNeighborLists = bondedNeighborLists,
      useNeutronDistances=False,probePhil=probePhil)
  extra = ret.extraAtomInfo

  # Also compute the maximum VDW radius among all atoms.
  maxVDWRad = 1
  for a in atoms:
    maxVDWRad = max(maxVDWRad, extra.getMappingFor(a).vdwRadius)

  # Optimization will place the movers. Make sure we got as many as we expected.
  # Make sure that the orientation for all of the movers is correct.
  probePhil = _philLike()
  opt = Optimizer(probePhil, True, model, modelIndex = None, altID = None,
                bondedNeighborDepth = 4,
                useNeutronDistances = False,
                minOccupancy = 0.02,
                preferenceMagnitude = 1.0,
                nonFlipPreference = 0.5,
                skipBondFixup = False,
                flipStates = '',
                verbosity = 1)
  numMovers = _numMovers(opt)
  if numMovers != 1:
    return "Optimizers.Test(): Incorrect number of Movers for B-factor test: " + str(numMovers)
  res = re.findall(r'pose Angle [-+]?\d+', opt.getInfo())
  for r in res:
    if not 'pose Angle 90' in r:
      return "Optimizers.Test(): Unexpected angle ("+str(r)+") for B-factor test"

  #========================================================================
  # Several ignored waters (B factor or occupancy) between kept ones: they leave the
  # Optimizer's atom list, the remaining atoms keep their order, and the model keeps them.
  pdb_waters = (
"""
CRYST1   93.586  127.886  251.681  90.00  90.00  90.00 I 2 2 2
HETATM    1  C   ACT A   1       6.494 -47.273 -37.006  1.00 16.65           C
HETATM    2  O   ACT A   1       5.654 -47.981 -37.645  1.00 15.33           O
HETATM    3  CH3 ACT A   1       6.790 -47.410 -35.548  1.00 17.13           C
HETATM    4  OXT ACT A   1       7.110 -46.425 -37.699  1.00 17.05           O
HETATM  200  O   HOH A 200      12.200 -48.610 -35.148  1.00 50.00           O
HETATM  201  O   HOH A 201       1.800 -46.810 -34.548  1.00 20.00           O
HETATM  202  O   HOH A 202      12.200 -42.610 -35.148  0.50 20.00           O
HETATM  203  O   HOH A 203       1.800 -52.810 -34.548  1.00 20.00           O
HETATM  204  O   HOH A 204       6.200 -52.610 -41.148  1.00 45.00           O
END
"""
    )
  dm = DataManager(['model'])
  dm.process_model_str("pdb_waters.pdb", pdb_waters)
  model = shift_and_box_model(model = dm.get_model())
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model)
  reduce_add_h_obj.run()
  model = reduce_add_h_obj.get_model()
  model.process(make_restraints=True, pdb_interpretation_params=p)
  expected = [a.memory_id() for a in GetAtomsForConformer(model.get_hierarchy().models()[0], "")
    if a.parent().resname.strip() != 'HOH' or a.parent().parent().resseq_as_int() in (201, 203)]
  def waterOxygens():
    return [a.memory_id() for a in model.get_hierarchy().atoms()
      if a.parent().resname.strip() == 'HOH' and a.element.strip().upper() == 'O']
  watersBefore = waterOxygens()
  opt = Optimizer(_philLike(), True, model, modelIndex = None, altID = None, verbosity = 1)
  # Phantom Hydrogens are appended after the model's atoms
  got = [a.memory_id() for a in opt._atoms if not a.name.strip().startswith('H?')]
  if got != expected:
    return "Optimizers.Test(): Wrong atoms or order after ignoring waters"
  if not "Ignored 3 waters due to occupancy or B factor" in opt.getInfo():
    return "Optimizers.Test(): Missing ignored-waters line"
  if waterOxygens() != watersBefore:
    return "Optimizers.Test(): Ignored waters were removed from the model"

  #========================================================================
  # Generate an example data model with a small molecule in it or else read
  # from the specified file.
  infoString = ""
  if inFileName is not None and len(inFileName) > 0:
    # Read a model from a file using the DataManager
    print('Reading model from',inFileName)
    dm = DataManager()
    dm.process_model_file(inFileName)
    model = dm.get_model(inFileName)
  else:
    # Generate a small-molecule model using the map model manager
    print('Generating model')
    mmm=map_model_manager()         #   get an initialized instance of the map_model_manager
    mmm.generate_map()              #   get a model from a generated small library model and calculate a map for it
    model = mmm.model()             #   get the model

  # Add Hydrogens to the model
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model)
  reduce_add_h_obj.run()
  model = reduce_add_h_obj.get_model()

  # Interpret the model after shifting and adding Hydrogens to it so that
  # all of the needed fields are filled in when we use them below.
  # @todo Remove this once place_hydrogens() does all the interpretation we need.
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.allow_polymer_cross_special_position=True
  p.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
  p.pdb_interpretation.proceed_with_excessive_length_bonds=True
  model.process(make_restraints=True, pdb_interpretation_params=p) # make restraints

  # Run the optimizer on the model to make sure it doesn't crash.
  opt = Optimizer(probePhil, True, model)

  # Write debugging output if we've been asked to
  if dumpAtoms:
    fn = "deleteme.pdb"
    fn = model.get_hierarchy().write_pdb_or_mmcif_file(
        target_format = 'pdb',
        target_filename = fn)
    print("Wrote model to '%s'" %fn)
    f = open("atomDump.pdb","w")
    f.write(opt.getAtomDump())

  #========================================================================
  # @todo Unit test a multi-model case, a multi-alternate case, and singles of each.

  return ""

##################################################################################
# If we're run on the command line, test our classes and functions.
if __name__ == '__main__':

  #==============================================================
  # Parse command-line arguments.  The 0th argument is the name
  # of the script. There can be the name of a PDB file to read.
  parser = argparse.ArgumentParser(description='Test mmtbx.reduce.Optimizers.')
  parser.add_argument("--dumpAtoms", help="dump the atoms into PDB files to help debug", action="store_true")
  parser.add_argument('inputFile', nargs='?', default="")
  args = parser.parse_args()

  ret = Test(args.inputFile, args.dumpAtoms)
  if len(ret) == 0:
    print('Success!')
  else:
    print(ret)

  assert (len(ret) == 0)
