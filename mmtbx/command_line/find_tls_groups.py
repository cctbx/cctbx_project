# LIBTBX_SET_DISPATCHER_NAME phenix.find_tls_groups

import sys, time
import iotbx.pdb
import mmtbx.tls
from mmtbx.tls import tools
from scitbx.array_family import flex
import mmtbx.secondary_structure
import random
from copy import deepcopy
from cctbx import adptbx

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

def all_permutations(N):
  unique_set = list(xrange(N))
  array = [[i] for i in unique_set]
  gss = list(range(2,N))
  gss.reverse()
  def is_all(x, us):
    tmp = []
    for i in x:
      for j in i:
        if(not j in tmp): tmp.append(j)
    tmp.sort()
    return tmp == us
  def not_in(v, t):
    for i in v:
      if i in t:
        return False
    return True
  result = [array[:]]
  for gs in gss:
    #print "group by:",gs
    for start_index in xrange(N):
      tmp = []
      tmp_ = []
      for i,a in enumerate(array):
        if(i>=start_index and i<=start_index+gs-1):
          tmp_.append(a[0])
        else: tmp.append(a)
        if(len(tmp_)==gs):
          tmp.append(tmp_)
          tmp_ = []
      tmp__=[]
      if(is_all(tmp, unique_set)):
        #print "  tmp=",tmp
        result.append(tmp)
        tmp__.append(tmp)
        for i,ti in enumerate(tmp):
          ztmp = []
          for j,tj in enumerate(tmp):
            if(j>i and abs(j-i)==1): ztmp.append([ti[0],tj[0]])
            elif(i!=j): ztmp.append(tj)
          if(not ztmp in result and is_all(ztmp, unique_set)):
            result.append(ztmp)
            #print "     ",ztmp
            tmp__.append(ztmp)
  def h1(x):
    tmp = []
    for i in x:
      for j in i:
        if not j in tmp: tmp.append(j)
    return tmp
  for i, xtmp_ in enumerate(result):
    for j, ytmp_ in enumerate(result):
      if(j>i):
        for xtmp__ in xtmp_:
          tmp = []
          tmp.append(xtmp__)
          tmpl = h1(tmp)
          for ytmp__ in ytmp_:
            if(ytmp__ != xtmp__ and not_in(ytmp__, tmpl)):
              tmp.append(ytmp__)
              if(is_all(tmp, unique_set)):
                tmp.sort()
                if(not tmp in result):
                  result.append(tmp)
                  tmp=[]
                  break
  ##########
  print "All permutations:",len(result)
  n = []
  for r_ in result:
    n_=[]
    #print r_
    for r__ in r_:
      #print "  ", r__
      r__.sort()
      n_.append(r__)
    n_.sort()
    n.append(n_)
  result = n
  result.sort()
  print "finally:", len(result)


  ### DEBUG
  for i,ri in enumerate(result):
    for j,rj in enumerate(result):
      if(i!=j): assert ri != rj

  ###
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
  chs2=0
  for i, s in enumerate(sels):
    for s_ in s: chs2 += s_[0].size()
  assert min([chs1,chs2,cntr3]) == max([chs1,chs2,cntr3])
  return sels

def regroup_groups(sels, residues, fragment_size, max_sels=15):
  sel = []
  min_size = 1.e+6
  i_min = None
  chsum1 = 0
  for i_seq, s in enumerate(sels):
    chsum1 += len(s)
    if(len(s)<min_size):
      i_min = i_seq
      min_size = len(s)
  while (len(sels)>max_sels or min_size < fragment_size):
    new_sels = []
    min_size = 1.e+6
    i_min = None
    for i_seq, s in enumerate(sels):
      if(len(s)<min_size):
        i_min = i_seq
        min_size = len(s)
    if(min_size < fragment_size or len(sels)>max_sels):
      l,r = None,None
      if(i_min-1>=0): l=sels[i_min-1]
      if(i_min+1<len(sels)): r=sels[i_min+1]
      if([l,r].count(None)==0):
        if(len(l)<len(r)):
          x = deepcopy(sels[i_min-1])
          y = deepcopy(sels[i_min])
          x.extend(y)
          sels[i_min-1]=x
          sels = sels[:i_min]+sels[i_min+1:]
        else:
          x = deepcopy(sels[i_min])
          y = deepcopy(sels[i_min+1])
          x.extend(y)
          sels[i_min]=x
          sels = sels[:i_min+1]+sels[i_min+2:]
      elif(l is not None):
          x = deepcopy(sels[i_min-1])
          y = deepcopy(sels[i_min])
          x.extend(y)
          sels[i_min-1]=x
          sels = sels[:i_min]+sels[i_min+1:]
      elif(r is not None):
          x = deepcopy(sels[i_min])
          y = deepcopy(sels[i_min+1])
          x.extend(y)
          sels[i_min]=x
          sels = sels[:i_min+1]+sels[i_min+2:]
  chsum2 = 0
  for i_seq, s in enumerate(sels):
    chsum2 += len(s)
  assert chsum1 == chsum2, [chsum1, chsum2]
  return sels


def get_model_partitioning(residues, secondary_structure_selection):
  fragment_size = 5
  print "  Grouping residues by secondary structure..."
  sels = group_residues(residues)
  print "  Fragment size:", fragment_size
  for i, s in enumerate(sels):
    print "      group %d: %d"%(i, len(s))
  print "  Re-grouping to achieve maximum possible nuber of groups..."
  new_sels = regroup_groups(sels, residues, fragment_size)
  print "  Fragment size:", fragment_size
  for i, s in enumerate(new_sels):
    print "      group %d: %d"%(i, len(s))
  len_new_sels = len(new_sels)
  if(len_new_sels==10):
    from mmtbx.tls import perm10
    perms = perm10.res
  elif(len_new_sels==11):
    from mmtbx.tls import perm11
    perms = perm11.res
  elif(len_new_sels==12):
    from mmtbx.tls import perm12
    perms = perm12.res
  elif(len_new_sels==13):
    from mmtbx.tls import perm13
    perms = perm13.res
  elif(len_new_sels==14):
    from mmtbx.tls import perm14
    perms = perm14.res
  elif(len_new_sels==15):
    from mmtbx.tls import perm15
    perms = perm15.res
  elif(len_new_sels<10):
    perms = all_permutations(len(list(xrange(len_new_sels))))
  else: raise RuntimeError("Too many permutations.")

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
    for p_ in p:
      for g in groups[p_]:
        one_group.extend(g[0])
    result.append(one_group)
  return result

def tls_refinery(sites_cart, selection, u_cart=None, u_iso=None,  max_iterations=30):
  sites_cart_ = sites_cart.select(selection)
  cm = sites_cart_.mean_weighted(weights=flex.double(selection.size(),1))
  assert [u_cart, u_iso].count(None)==1
  if(u_cart is not None):
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
  else:
    return tools.tls_from_uiso_minimizer(
      uiso         = u_iso.select(selection),
      T_initial      = [0],
      L_initial      = [0,0,0,0,0,0],
      S_initial      = [0,0,0],
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

def tls_refinery_random_groups(sites_cart, n_groups, u_cart=None, u_iso=None, n_runs=20):
  assert [u_cart, u_iso].count(None)==1
  t = 0
  for tr in xrange(n_runs):
    selections = chunks(size=sites_cart.size(), n_groups=n_groups)
    for selection in selections:
      mo = tls_refinery(u_cart=u_cart, u_iso=u_iso, sites_cart=sites_cart, selection=selection)
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
    for p_ in p:
      for g in groups[p_]:
        one_group.append(g[2])
    resseq = "resseq %s:%s"%(one_group[0].strip(),
      one_group[len(one_group)-1].strip())
    result.append(resseq)
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
  u_cart = None#xray_structure.scatterers().extract_u_cart(unit_cell)
  u_iso  = xray_structure.extract_u_iso_or_u_equiv()#*adptbx.u_as_b(1.)
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
  if(u_cart is not None):
    assert secondary_structure_selection.size() == u_cart.size()
  else:
    assert secondary_structure_selection.size() == u_iso.size()
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
    for i_perm, perm in enumerate(perms):
      if i_perm%100==0:
        print "    ...perm %d of %d"%(i_perm, len(perms))
      selections = tls_group_selections(groups, perm)
      target = 0
      for selection in selections:
        mo = tls_refinery(
          u_cart     = u_cart,
          u_iso      = u_iso,
          sites_cart = sites_cart,
          selection  = selection)
        target += mo.f
      dic.setdefault(len(perm), []).append([target,perm])
      #print "    perm %d of %d: target=%8.3f (TLS groups: %s), permutation:"%(
      #  i_perm, len(perms),target,len(perm)),perm
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
      if(u_cart is not None):
        u_cart_ = u_cart.select(chain_selection)
      else: u_cart_ = None
      if(u_iso is not None):
        u_iso_ = u_iso.select(chain_selection)
      else: u_iso_ = None
      r = tls_refinery_random_groups(
        u_cart     = u_cart_,
        u_iso      = u_iso_,
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
  t0 = time.time()
  run(args=sys.argv[1:])
  print "Time: %10.3f"%(time.time()-t0)
