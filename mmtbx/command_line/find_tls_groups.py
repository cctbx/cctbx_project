# LIBTBX_SET_DISPATCHER_NAME phenix.find_tls_groups

import sys
import iotbx.pdb
from mmtbx.tls import tools
from scitbx.array_family import flex
import mmtbx.secondary_structure
import libtbx.math_utils
import random

##### PERMTOOLS
def consequtive_permutations(iterable, r=None):
  pool = tuple(iterable)
  n = len(pool)
  r = n if r is None else r
  if r > n: return
  indices = range(n)
  cycles = range(n, n-r, -1)
  yield tuple(pool[i] for i in indices[:r])
  while n:
    for i in reversed(range(r)):
      cycles[i] -= 1
      if cycles[i] == 0:
        indices[i:] = indices[i+1:] + indices[i:i+1]
        cycles[i] = n - i
      else:
        j = cycles[i]
        indices[i], indices[-j] = indices[-j], indices[i]
        tmp = []
        good=True
        for k in indices[:r]:
          x = pool[k]
          ltmp = len(tmp)
          if(ltmp>0 and tmp[ltmp-1]-x!=-1):
            good=False
            break
          tmp.append(x)
        if(good): yield tuple(pool[i] for i in indices[:r])
        break
    else:
      return

def all_permutations(x):
  result = [x[:]]
  perm = []
  for ngr in range(1,len(x)):
    r = list(consequtive_permutations(x, ngr))
    #print ngr,r
    perm.append(r)
  #print
  for r in perm[1:]:
    for r_ in r:
      res = [r_]
      for x_ in x:
        if(not x_ in r_): res.append(x_)
      result.append(res)
  #for r in result:
  #  print r
  #print
  unique_set = []
  for r in result:
    for rr in r:
      if(not rr in unique_set):
        unique_set.append(rr)
  tmp = []
  for r in unique_set:
    if(isinstance(r,tuple)): tmp.append(list(r))
    else: tmp.append([r])
  unique_set = tmp[:]
  #print unique_set, len(unique_set)
  def reinitialize(x, unique_set):
    tmp1 = []
    for i in unique_set:
      for j in i:
        if(not j in tmp1): tmp1.append(j)
    tmp1.sort()
    tmp2 = []
    for i in x:
      for j in i:
        if(not j in tmp2): tmp2.append(j)
    tmp2.sort()
    if(tmp1==tmp2): return True
    else: return False
  result = []
  for i in unique_set:
    #print
    #print "i:", i
    res = [i[:]]
    for j in unique_set:
      if(i!=j):
        #print "  j:", j, res
        good=True
        for jj in j:
          if(jj in i):
            good=False
            break
        if(good and not j in res):
          for jjj in j:
            for res_ in res:
              for res__ in res_:
                if(jjj == res__):
                  good=False
                  break
          if(good):
            res.append(j)
            res.sort()
            if(reinitialize(res, unique_set)):
              if(not res in result): result.append(res[:])
              #print "---res:", res
              res = [i[:]]
  print "  All permutations: %s"%str(len(result))
  #for r in result:
  #  print "    ",r
  return result
###############

def group_residues(residues):
  sels, sel = [], []
  cntr = 0
  chs1=0
  cntr3 = 0
  for i, r in enumerate(residues):
    chs1 += r[0].size()
    next = i+1
    if(next<len(residues)):
      if(r[1] == residues[next][1]):
        sel.append(r)
        #print i, r[1], cntr
      else:
        sel.append(r)
        sels.append(sel)
        for ss in sel: cntr3 += ss[0].size()
        sel = []
        #print i, r[1], cntr
        cntr += 1
    else:
      sel.append(r)
      sels.append(sel)
      for ss in sel: cntr3 += ss[0].size()
      sel = []
      #print i, r[1], cntr
      cntr += 1
  print "    Number of group candidates: %d"%len(sels)
  print "    Number of residues in group candidates:"
  chs2=0
  for i, s in enumerate(sels):
    print "      group %d: %d"%(i, len(s))
    for s_ in s: chs2 += s_[0].size()
  assert min([chs1,chs2,cntr3]) == max([chs1,chs2,cntr3])
  return sels

def regroup_groups(sels, residues, fragment_size):
  new_sels = []
  sel = []
  for i in xrange(len(sels)):
    sel.extend(sels[i])
    if(len(sel)>=fragment_size):
      new_sels.append(sel[:])
      sel = []
  if(len(sel)>0):
    if(len(sel)>=fragment_size):
      new_sels.append(sel[:])
    else:
      new_sels[len(new_sels)-1].extend(sel[:])
  check_sum_size1 = 0
  print "    Number of groups: %d"%len(new_sels)
  print "    Number of residues in groups:"
  for i, s in enumerate(new_sels):
    print "      group %d: %d"%(i, len(s))
    for s_ in s:
      check_sum_size1 += s_[0].size()
  check_sum_size2 = 0
  for r in residues:
    check_sum_size2 += r[0].size()
  assert check_sum_size1 == check_sum_size2,[check_sum_size1, check_sum_size2]
  return new_sels


def get_model_partitioning(residues, secondary_structure_selection,
                           fragment_size=5):
  print "  Grouping residues by secondary structure..."
  sels = group_residues(residues)
  print "  Re-grouping to achieve minimum requested fragment size (%s residues)..."%\
    str(fragment_size)
  fragment_size = 10
  new_sels = sels[:]
  while len(new_sels)>10:
    print "   Trying fragment size:", fragment_size
    new_sels = regroup_groups(sels, residues, fragment_size)
    fragment_size += 1
  perms = all_permutations(x = list(xrange(len(new_sels))))
  return new_sels, perms

def chains_and_atoms(pdb_hierarchy, secondary_structure_selection):
  new_secondary_structure_selection = flex.bool()
  get_class = iotbx.pdb.common_residue_names_get_class
  chains_and_residue_selections = []
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      result = []
      for rg in chain.residue_groups():
        result_ = flex.size_t()
        is_secondary_structure = False
        for ag in rg.atom_groups():
          if(get_class(name=ag.resname) == "common_amino_acid" or
             get_class(name=ag.resname) == "common_rna_dna"):
            for atom in ag.atoms():
              result_.append(atom.i_seq)
              if(not is_secondary_structure):
                is_secondary_structure = \
                  secondary_structure_selection[atom.i_seq]
              new_secondary_structure_selection.append(
                secondary_structure_selection[atom.i_seq])
        if(result_.size()>0):
          result.append([result_, is_secondary_structure, rg.resseq,
                         rg.unique_resnames()])
      if(len(result)>0):
        chains_and_residue_selections.append([chain.id, result])
  print "Considering these chains:"
  for ch in chains_and_residue_selections:
    print "  chain %s (number of residues selected: %d)" % (ch[0], len(ch[1]))
  return chains_and_residue_selections, new_secondary_structure_selection

def tls_group_selections(groups, perm):
  result = []
  for p in perm:
    one_group = flex.size_t()
    #print "p:", p
    for p_ in p:
      #print "p_",groups[p_]
      for g in groups[p_]:
        one_group.extend(g[0])
    result.append(one_group)
  #print result
  return result

def tls_refinery(u_cart, sites_cart, selection, max_iterations=50):
  sites_cart_ = sites_cart.select(selection)
  cm = sites_cart_.mean_weighted(weights=flex.double(selection.size(),1))
  return tools.tls_from_uaniso_minimizer(
    uaniso         = u_cart.select(selection),
    T_initial      = [0,0,0,0,0,0],
    L_initial      = [0,0,0,0,0,0],
    S_initial      = [0,0,0,0,0,0,0,0,0],
    refine_T       = True,
    refine_L       = True,
    refine_S       = True,
    origin         = cm,
    sites          = sites_cart_,
    max_iterations = max_iterations)

def chunks(size, n_groups):
  rp = list(xrange(size))
  chunk_size = size/n_groups
  nc = chunk_size
  counter = 0
  sum_size = 0
  res = []
  while nc <= size:
    next = nc+chunk_size
    if(next>size or next+chunk_size>size): next = size
    if(counter==n_groups and nc+chunk_size>size): break
    if(len(res)>n_groups-1): break
    r = random.randrange(nc,next)
    try: ev = size-1-r>1 and r-max(res)>1
    except: ev = size-1-r>1 and not r in res
    if(ev):
      res.append(r)
      nc+=chunk_size
      if(len(res)>n_groups-2): break
  result = []
  for i, r in enumerate(res):
   if i==0: result.append([0,r])
   elif(i==len(res)): result.append([r,size])
   else: result.append([res[i-1]+1,res[i]])
  result.append([res[len(res)-1]+1,size-1])
  tmp = []
  for r in result:
    r_ = flex.size_t(range(r[0],r[1]))
    assert r_.size() > 0, [result, res]
    tmp.append(r_)
  return tmp

def tls_refinery_random_groups(u_cart, sites_cart, n_groups, n_runs=20):
  t = 0
  for tr in xrange(n_runs):
    selections = chunks(size=u_cart.size(), n_groups=n_groups)
    for selection in selections:
      mo = tls_refinery(u_cart=u_cart, sites_cart=sites_cart, selection=selection)
      t += mo.f
  return t/n_runs

def chain_selection_from_residues(residues):
  chain_selection = flex.size_t()
  for r in residues:
    chain_selection.extend(r[0])
  return chain_selection

def permutations_as_atom_selection_string(groups, perm):
  result = []
  for p in perm:
    one_group = []
    #print "p:", p
    for p_ in p:
      #print "p_",groups[p_]
      for g in groups[p_]:
        one_group.append(g[2])
    resseq = "resseq %s:%s"%(one_group[0].strip(),
      one_group[len(one_group)-1].strip())
    result.append(resseq)
    #print one_group
  #print result
  return result

def run(args):
  default_message="""\

phenix.find_tls_groups: Tool for automated partitioning a model into TLS groups.

Usage:
  phenix.find_tls_groups model.pdb
  """
  if(len(args) != 1):
    print default_message
    return
  pdb_file_name = args[0]
  if(not iotbx.pdb.is_pdb_file(pdb_file_name)):
    print "A PDB file is required."
    return
  pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  pdb_atoms = pdb_hierarchy.atoms()
  pdb_atoms.reset_i_seq()
  #
  xray_structure = pdb_inp.xray_structure_simple()
  xray_structure.convert_to_anisotropic()
  sites_cart = xray_structure.sites_cart()
  cm = sites_cart.mean_weighted(weights=flex.double(sites_cart.size(),1))
  unit_cell = xray_structure.unit_cell()
  u_cart = xray_structure.scatterers().extract_u_cart(unit_cell)
  #
  ssm = mmtbx.secondary_structure.manager(
    pdb_hierarchy                = pdb_hierarchy,
    xray_structure               = xray_structure,
    sec_str_from_pdb_file        = None,
    params                       = None,
    assume_hydrogens_all_missing = None,
    tmp_dir                      = None)
  ssm.find_automatically()
  alpha_h_selection = ssm.alpha_selection()
  secondary_structure_selection = ssm.alpha_selection() | \
      ssm.beta_selection() | ssm.base_pair_selection()
  assert secondary_structure_selection.size() == u_cart.size()
  ssm.show_summary()
  chains_and_residue_selections, secondary_structure_selection = chains_and_atoms(
    pdb_hierarchy                 = pdb_hierarchy,
    secondary_structure_selection = secondary_structure_selection)
  chains_and_permutations = []
  chains_and_atom_selection_strings = []
  for crs in chains_and_residue_selections:
    print "Processing chain %s:"%crs[0]
    chain_selection = chain_selection_from_residues(crs[1])
    groups, perms = get_model_partitioning(residues = crs[1],
      secondary_structure_selection = secondary_structure_selection)
    print "  Fitting TLS matrices..."
    dic = {}
    target_best = 1.e+9
    for perm in perms:
      selections = tls_group_selections(groups, perm)
      target = 0
      for selection in selections:
        mo = tls_refinery(
          u_cart     = u_cart,
          sites_cart = sites_cart,
          selection  = selection)
        target += mo.f
      dic.setdefault(len(perm), []).append([target,perm])
      print "    target=%10.3f (TLS groups: %s), permutation:"%(
        target,len(perm)),perm
    print "    Best fits:"
    print "      No. of         Targets"
    print "      groups   best   rand.pick diff.  score permutation"
    score_best = -1.e+9
    perm_choice = None
    for k, v in zip(dic.keys(),dic.values()):
      t_best = v[0][0]
      perm_best = v[0][1]
      for v_ in v:
        if(v_[0]<t_best):
          t_best = v_[0]
          perm_best = v_[1]
      r = tls_refinery_random_groups(
        u_cart     = u_cart.select(chain_selection),
        sites_cart = sites_cart.select(chain_selection),
        n_groups   = k)
      score = (r-t_best)/(r+t_best)*100.
      print "         %3d   %6.1f   %6.1f %6.1f %6.1f"%(
        k,t_best, r, r-t_best, score), perm_best
      if(score > score_best):
        score_best = score
        perm_choice = perm_best[:]
    #
    chains_and_permutations.append([crs[0],perm_choice])
    chains_and_atom_selection_strings.append([crs[0],
      permutations_as_atom_selection_string(groups, perm_choice)])
    #
  print
  print "%sSUMMARY%s"%("-"*36,"-"*37)
  print
  print "Optimal TLS groups:"
  for chain_and_permutation in chains_and_permutations:
    print chain_and_permutation
  print
  print "TLS groups (atom selection strings):"
  for r in chains_and_atom_selection_strings:
    prefix = "chain %s and "%r[0]
    for r_ in r[1:]:
      for r__ in r_:
        print prefix+"(%s)"%r__
  print

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
