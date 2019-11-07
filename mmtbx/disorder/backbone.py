
"""
Tools for analyzing backbone motion in ensembles and statically disordered
structures.
"""

# derived from BackrubFinder2.java by Ian Davis

from __future__ import absolute_import, division, print_function
from libtbx import slots_getstate_setstate
import math
import sys
from six.moves import range

def get_calphas(pdb_hierarchy):
  n_models = pdb_hierarchy.models_size()
  calphas = []
  for n in range(n_models) : calphas.append([])
  for i_mod, model in enumerate(pdb_hierarchy.models()):
    for chain in model.chains():
      if (not chain.is_protein()) : continue
      for residue_group in chain.residue_groups():
        atom_group = residue_group.only_atom_group()
        c_alpha = atom_group.get_atom("CA")
        if (c_alpha is not None):
          #print "%d: %s" % (i_mod, c_alpha.id_str())
          calphas[i_mod].append(c_alpha)
        else :
          calphas[i_mod].append(None)
  for calphas_model in calphas :
    assert len(calphas_model) == len(calphas[0])
  return calphas

def get_cbetas(pdb_hierarchy):
  n_models = pdb_hierarchy.models_size()
  cbetas = [ [] ] * n_models
  for i_mod, model in enumerate(pdb_hierarchy.models()):
    for chain in model.chains():
      if (not chain.is_protein()) : continue
      for residue_group in chain.residue_groups():
        atom_group = residue_group.only_atom_group()
        c_beta = atom_group.get_atom("CB")
        if (c_beta is None):
          c_beta = atom_group.get_atom("2HA")
        if (c_beta is not None):
          cbetas[i_mod].append(c_beta)
        else :
          cbetas[i_mod].append(None)
  for cbetas_model in cbetas :
    assert len(cbetas_model) == len(cbetas[0])
  return cbetas

def circ_stddev(t, deg=True):
  assert (len(t) > 0)
  from scitbx.array_family import flex
  if (not deg):
    t = t * 180 / math.pi
  mean = flex.mean(t)
  a = sa = 0
  # FIXME use C++ array operations
  for i in range(len(t)):
    a = abs(mean - t[i]) % 360
    if (a > 180.0):
      a = 360.0 - a
    sa += a*a
  return math.sqrt(sa / len(t))

def circ_len(t, deg=True):
  assert (len(t) > 0)
  from scitbx.array_family import flex
  if (deg):
    t = math.pi * (t/180)
  sx = flex.sum(flex.cos(t)) / len(t)
  sy = flex.sum(flex.sin(t)) / len(t)
  return math.sqrt(sx**2 + sy**2)

def circ_mean(t, deg=True):
  assert (len(t) > 0)
  from scitbx.array_family import flex
  if (deg):
    t = math.pi * (t/180)
  sx = flex.sum(flex.cos(t)) / len(t)
  sy = flex.sum(flex.sin(t)) / len(t)
  return math.degrees(math.atan2(sy, sx))

class backrub_residue(slots_getstate_setstate):
  __slots__ = ["id_str", "i_mod", "j_mod", "rmsd", "backrub_angle", "xyz"]
  def __init__(self, calpha, i_mod, j_mod, rmsd, backrub_angle):
    residue_group = calpha.parent().parent()
    self.id_str = residue_group.id_str()
    self.i_mod = i_mod
    self.j_mod = j_mod
    self.rmsd = rmsd
    self.backrub_angle = backrub_angle
    self.xyz = calpha.xyz

  def show(self, out=sys.stdout, prefix=""):
    print(prefix+"backrub %s (%s,%s): angle=%.1f" % \
      (self.id_str, self.i_mod, self.j_mod, self.backrub_angle), file=out)

def evaluate_backrub_pair_impl(
    calphas_A,
    calphas_B,
    labels=(),
    max_calpha_sep=5.0,
    rmsd_limit=0.1,
    backrub_angle_limit=10.0) : # FIXME is this an appropriate cutoff?
  assert (len(calphas_A) == len(calphas_B) == 5)
  if (None in calphas_A) or (None in calphas_B):
    return None
  for k_res in range(0, 4):
    dist = calphas_A[k_res].distance(calphas_A[k_res+1])
    if (dist > max_calpha_sep):
      return None
  from scitbx.array_family import flex
  from scitbx.math import superpose
  from scitbx.matrix import col
  import scitbx.math
  sites_A = flex.vec3_double([ calphas_A[k].xyz for k in [0,1,3,4] ])
  sites_B = flex.vec3_double([ calphas_B[k].xyz for k in [0,1,3,4] ])
  lsq_fit = superpose.least_squares_fit(
    reference_sites=sites_A,
    other_sites=sites_B)
  sites_B_new = lsq_fit.other_sites_best_fit()
  rmsd = sites_B_new.rms_difference(sites_A)
  ca2 = (col(sites_A[1]) + col(sites_B_new[1])) / 2
  ca3r = col(calphas_A[2].xyz)
  ca3m = lsq_fit.rt() * calphas_B[2].xyz
  ca4 = (col(sites_A[2]) + col(sites_B_new[2])) / 2
  backrub_angle = scitbx.math.dihedral_angle(
    sites=[ca3r.elems, ca2.elems, ca4.elems, ca3m.elems],
    deg=True)
  if ((rmsd <= rmsd_limit) and
      (abs(backrub_angle) >= backrub_angle_limit)):
    if (len(labels) == 0):
      labels = (calphas_A[2].fetch_labels().altloc,
                calphas_B[2].fetch_labels().altloc)
    return backrub_residue(
      calpha=calphas_A[2],
      i_mod=labels[0],
      j_mod=labels[1],
      rmsd=rmsd,
      backrub_angle=backrub_angle)
  return None

def find_backrubs(
    pdb_hierarchy=None,
    residue_group=None,
    max_calpha_sep=5.0,
    rmsd_limit=0.1,
    backrub_angle_limit=10.0):
  assert ([pdb_hierarchy, residue_group].count(None) == 1)
  if (residue_group is not None):
    return find_ensemble_backrubs(
      residue_group=residue_group,
      max_calpha_sep=max_calpha_sep,
      rmsd_limit=rmsd_limit,
      backrub_angle_limit=backrub_angle_limit)
  backrubs = []
  for chain in pdb_hierarchy.only_model().chains():
    if (not chain.is_protein()):
      continue
    residue_groups = chain.residue_groups()
    for residue_group in residue_groups[2:-2] :
      if (residue_group.atom_groups_size() == 1):
        continue
      br = find_ensemble_backrubs(
        residue_group=residue_group,
        max_calpha_sep=max_calpha_sep,
        rmsd_limit=rmsd_limit,
        backrub_angle_limit=backrub_angle_limit)
      if (br is not None):
        backrubs.extend(br)
  return backrubs

def find_ensemble_backrubs(
    pdb_hierarchy=None,
    residue_group=None,
    max_calpha_sep=5.0,
    rmsd_limit=0.1,
    backrub_angle_limit=10.0):
  assert ([pdb_hierarchy, residue_group].count(None) == 1)
  from scitbx.array_family import flex
  from scitbx.math import superpose
  from scitbx.matrix import col
  import scitbx.math
  if (residue_group is not None):
    pdb_hierarchy = extract_backrub_residue_groups(residue_group)
    if (pdb_hierarchy is None):
      return []
  models = pdb_hierarchy.models()
  n_models = len(models)
  assert (n_models > 1)
  calphas = get_calphas(pdb_hierarchy)
  assert (len(calphas) == n_models)
  backrubs = []
  def get_labels(i_mod, j_mod):
    return (models[i_mod].id.strip(), models[j_mod].id.strip())
  for k_res in range(2, len(calphas[0]) - 2):
    for i_mod in range(n_models-1):
      for j_mod in range(i_mod+1, n_models):
        calphas_i = [ calphas[i_mod][k] for k in range(k_res-2, k_res+3) ]
        calphas_j = [ calphas[j_mod][k] for k in range(k_res-2, k_res+3) ]
        if (not None in calphas_i) and (not None in calphas_j):
          br = evaluate_backrub_pair_impl(
            calphas_A=calphas_i,
            calphas_B=calphas_j,
            labels=get_labels(i_mod, j_mod),
            max_calpha_sep=max_calpha_sep,
            rmsd_limit=rmsd_limit,
            backrub_angle_limit=backrub_angle_limit)
          if (br is not None):
            backrubs.append(br)
  return backrubs

#-----------------------------------------------------------------------
# UTILITY FUNCTIONS
def extract_backrub_residue_groups(residue_group,
    as_ensemble=True):
  import iotbx.pdb.hierarchy
  chain = residue_group.parent()
  root = chain.parent().parent()
  assert (len(root.models()) == 1)
  residue_groups = chain.residue_groups()
  backrub_residues = None
  for i_res, other_rg in enumerate(residue_groups):
    if (other_rg == residue_group):
      if (i_res < 2) or (i_res > len(residue_groups) - 3):
        return None
      backrub_residues = [ residue_groups[k] for k in range(i_res-2,i_res+3) ]
      break
  assert (backrub_residues is not None)
  root = iotbx.pdb.hierarchy.root()
  model = iotbx.pdb.hierarchy.model()
  root.append_model(model)
  new_chain = iotbx.pdb.hierarchy.chain(id=chain.id)
  model.append_chain(new_chain)
  for rg in backrub_residues :
    new_chain.append_residue_group(rg.detached_copy())
  if (as_ensemble):
    return alternate_conformations_as_multiple_models(root)
  return root

def alternate_conformations_as_multiple_models(pdb_hierarchy):
  from mmtbx.command_line import altloc_remediate
  import iotbx.pdb.hierarchy
  altlocs = set()
  for chain in pdb_hierarchy.only_model().chains():
    if (not chain.is_protein()) and (not chain.is_na()):
      continue
    for residue_group in chain.residue_groups():
      atom_groups = residue_group.atom_groups()
      if (len(atom_groups) > 1):
        if (atom_groups[0].altloc.strip() == ''):
          altloc_remediate.spread_to_residue(residue_group)
      altlocs.update(set([ a.altloc for a in residue_group.atom_groups() ]))
  if ('' in altlocs) : altlocs.remove('')
  n_confs = len(altlocs)
  if (n_confs <= 1):
    return pdb_hierarchy
  root = iotbx.pdb.hierarchy.root()
  altlocs = sorted(altlocs)
  models = {}
  for i_mod in range(n_confs):
    model = iotbx.pdb.hierarchy.model(id=altlocs[i_mod])
    models[altlocs[i_mod]] = model
    root.append_model(model)
  for chain in pdb_hierarchy.only_model().chains():
    if (not chain.is_protein()) and (not chain.is_na()):
      continue
    chains = {}
    for altloc in altlocs :
      new_chain = iotbx.pdb.hierarchy.chain(id=chain.id)
      models[altloc].append_chain(new_chain)
      chains[altloc] = new_chain
    for residue_group in chain.residue_groups():
      atom_groups = residue_group.atom_groups()
      if (len(atom_groups) == 1):
        for altloc in altlocs :
          chains[altloc].append_residue_group(residue_group.detached_copy())
      else :
        def add_atom_group(ag, new_chain):
          new_ag = ag.detached_copy()
          new_ag.altloc = ''
          new_rg = iotbx.pdb.hierarchy.residue_group(
            resseq=residue_group.resseq,
            icode=residue_group.icode)
          new_rg.append_atom_group(new_ag)
          new_chain.append_residue_group(new_rg)
        have_altlocs = set()
        for atom_group in atom_groups :
          altloc = atom_group.altloc
          add_atom_group(atom_group, chains[altloc])
          have_altlocs.add(altloc)
        for altloc in altlocs :
          if (not altloc in have_altlocs):
            add_atom_group(atom_groups[0], chains[altloc])
  return root
