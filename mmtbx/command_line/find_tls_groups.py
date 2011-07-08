# LIBTBX_SET_DISPATCHER_NAME phenix.find_tls_groups

from mmtbx.tls import tools
from mmtbx.refinement import print_statistics
import mmtbx.secondary_structure
import iotbx.pdb
from cctbx import adptbx
from scitbx.array_family import flex
import scitbx.linalg
import libtbx.phil
from libtbx.utils import Sorry
from libtbx import Auto
from copy import deepcopy
import cStringIO
import random
import os
import time
import sys


master_phil = libtbx.phil.parse("""
  pdb_file = None
    .type = path
  nproc = 1
    .type = int
  random_seed = 4865136
    .type = int
""")

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
  n = []
  for r_ in result:
    n_=[]
    for r__ in r_:
      r__.sort()
      n_.append(r__)
    n_.sort()
    n.append(n_)
  result = n
  result.sort()
  #print "  Number of all permutations:",len(result)
  ### DEBUG
  for i,ri in enumerate(result):
    for j,rj in enumerate(result):
      if(i!=j): assert ri != rj
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

def regroup_groups(sels, residues, fragment_size, max_sels):
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

def split_groups(sels, fragment_size):
  new_sels = []
  for sel in sels:
    if(len(sel)>fragment_size*3):
      #print "len(sel):", len(sel)
      is_ss_cntr=0
      for s_ in sel:
        #print s_
        if(s_[1]): is_ss_cntr += 1
      if(is_ss_cntr*100./len(sel)<50.):
        nc = 0
        new_sel = []
        while nc < len(sel):
          new_sels.append(sel[nc:nc+fragment_size])
          nc += fragment_size
      else:
        new_sels.append(sel)
    else:
      new_sels.append(sel)
  return new_sels

def show_groups(sels, out=None):
  if (out is None) :
    out = sys.stdout
  min_group_size = 1.e+6
  if 0: print >> out, "          Residues  Resseq  Sec.Structure"
  for i, s in enumerate(sels):
    is_ss_cntr=0
    for s_ in s:
      if(s_[1]): is_ss_cntr += 1
    if 0: print >> out, "      #%d: %8d  %s-%s %9d"%(i, len(s), s[0][2].strip(),
      s[len(s)-1][2].strip(), is_ss_cntr)
    if(len(s)<min_group_size): min_group_size = len(s)
  return min_group_size

def sels_as_selection_arrays(sels):
  result = []
  for s in sels:
    r_ = flex.size_t()
    for s_ in s:
      r_.extend(s_[0])
    if(r_.size()>0): result.append(r_)
  return result

def get_model_partitioning(residues,
                           secondary_structure_selection,
                           max_sels=13,
                           out=None) :
  if (out is None) :
    out = sys.stdout
  fragment_size = 5
  print >> out, "  Grouping residues by secondary structure..."
  sels = group_residues(residues)
  print >> out, "  Fragment size:", fragment_size
  print >> out, "    Initial groups..."
  min_group_size = show_groups(sels=sels)
  print >> out, "      n_groups=", len(sels)
  ###
  print >> out, "    Splitting groups is necesary..."
  print >> out, "      n_groups=", len(sels)
  sels = split_groups(sels = sels, fragment_size = fragment_size)
  show_groups(sels=sels)
  ###
  if(len(sels) > max_sels or min_group_size < fragment_size and len(sels)>1):
    print >> out, "  Re-grouping to achieve maximum possible nuber of groups..."
    sels = regroup_groups(sels, residues, fragment_size, max_sels)
    print >> out, "    n_groups=", len(sels)
    show_groups(sels=sels)
  len_new_sels = len(sels)
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
  elif(len_new_sels==1):
    perms = [[[0]]]
  elif(len_new_sels<10):
    perms = all_permutations(len(list(xrange(len_new_sels))))
  else: raise RuntimeError("Too many permutations.")
  return sels, perms

def chains_and_atoms(pdb_hierarchy, secondary_structure_selection,
    out=None) :
  if (out is None) :
    out = sys.stdout
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
          print >> out, ag.resname, get_class(name=ag.resname)
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
          result.append(
            [result_, is_secondary_structure, rg.resseq, rg.unique_resnames()])
      if(len(result)>0):
        chains_and_residue_selections.append([chain.id, result])
  print >> out, "Considering these chains:"
  for ch in chains_and_residue_selections:
    print >> out, "  chain '%s' (number of residues selected: %d)" % (ch[0], len(ch[1]))
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

def tls_refinery(sites_cart, selection, u_cart=None, u_iso=None,
                 use_minimizer=False, max_iterations=100):
  sites_cart_ = sites_cart.select(selection)
  cm = sites_cart_.mean_weighted(weights=flex.double(sites_cart_.size(),1))
  assert [u_cart, u_iso].count(None)==1
  if(not use_minimizer):
    obj = tools.tls_ls_derivative_coefficients(
      origin     = cm,
      sites_cart = sites_cart_,
      u_iso      = u_iso.select(selection))
    def s1(use_generalized_inverse=True):
      if(not use_generalized_inverse):
        obj.a.matrix_inversion_in_place()
        res = obj.a.matrix_multiply(obj.b)
      else:
        es = scitbx.linalg.eigensystem.real_symmetric(
          m=obj.a,
          relative_epsilon=1.e-12,
          absolute_epsilon=0)
        a = es.generalized_inverse_as_packed_u().matrix_packed_u_as_symmetric()
        res = a.matrix_multiply(obj.b)
      return res
    result = s1()
    target = tools.ls_target_from_iso_tls(
      t = result[0],
      l = tuple(result[1:7]),
      s = tuple(result[7:]),
      origin = cm,
      sites_cart = sites_cart_,
      u_isos = u_iso.select(selection))
    #print "target:",target
    class foo: pass
    foo.f = target
    return foo
  else:
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
      minimized = tools.tls_from_uiso_minimizer(
        uiso           = u_iso.select(selection),
        T_initial      = [0],
        L_initial      = [0,0,0,0,0,0],
        S_initial      = [0,0,0],
        refine_T       = True,
        refine_L       = True,
        refine_S       = True,
        origin         = cm,
        sites          = sites_cart_,
        max_iterations = max_iterations)
      if(0): # DEBUG
        print "Minimization:", xxx.f
        print "T_min:", minimized.T_min
        print "L_min:", minimized.L_min
        print "S_min:", minimized.S_min
        print
      return minimized

def chunks(size, n_groups):
  rp = list(xrange(size))
  chunk_size = size//n_groups
  nc = chunk_size
  counter = 0
  sum_size = 0
  res = []
  check = 0
  while nc <= size:
    check += 1
    chunk_size_ = chunk_size + random.randrange(-1,2)*int(0.5*chunk_size)
    next = nc+chunk_size_
    if check >100 or nc>=next or next>=size:
      ###
      chunk_size = size//n_groups
      nc = chunk_size
      counter = 0
      sum_size = 0
      res = []
      check = 0
      ###
      check=0
    if(next>size or next+chunk_size_>size): next = size
    if(counter==n_groups and nc+chunk_size_>size): break
    if(len(res)>n_groups-1): break
    r = random.randrange(nc,next)
    try: ev = size-1-r>1 and r-max(res)>1
    except Exception: ev = size-1-r>1 and not r in res
    if(ev):
      res.append(r)
      nc+=chunk_size_
      if(len(res)>n_groups-2): break
  result = []
  for i, r in enumerate(res):
   if i==0: result.append([0,r])
   elif(i==len(res)): result.append([r,size])
   else: result.append([res[i-1]+1,res[i]])
  result.append([res[len(res)-1]+1,size-1])
  tmp = []
  for r in result:
    a = r[0]
    b = r[1]
    if(a!=0): a = a
    b = b+1
    r_ = flex.size_t(range(a,b))
    assert r_.size() > 0, [result, res]
    tmp.append(r_)
  # DEBUG
  cntr = 0
  for s in tmp:
    cntr += s.size()
  assert cntr == size
  #
  return tmp

def tls_refinery_random_groups(sites_cart, n_groups, u_cart=None, u_iso=None, n_runs=50):
  assert [u_cart, u_iso].count(None)==1
  t = 0
  for tr in xrange(n_runs):
    while True:
      selections = chunks(size=sites_cart.size(), n_groups=n_groups)
      #print [(min(s),max(s)) for s in selections]
      if(len(selections) == n_groups): break
    assert len(selections) == n_groups
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

# XXX for multiprocessing
class analyze_permutations (object) :
  def __init__ (self, groups, sites_cart, u_cart, u_iso) :
    self.groups = groups
    self.sites_cart = sites_cart
    self.u_cart = u_cart
    self.u_iso = u_iso

  def __call__ (self, perm) :
    selections = tls_group_selections(self.groups, perm)
    target = 0
    for selection in selections:
      mo = tls_refinery(
        u_cart     = self.u_cart,
        u_iso      = self.u_iso,
        sites_cart = self.sites_cart,
        selection  = selection)
      target += mo.f
    return target

def run (args=(), params=None, pdb_hierarchy=None, xray_structure=None,
    out=None) :
  if (out is None) :
    out = sys.stdout
  print_statistics.make_header("phenix.find_tls_groups", out=out)
  default_message="""\

phenix.find_tls_groups: Tool for automated partitioning a model into TLS groups.

Usage:
  phenix.find_tls_groups model.pdb [nproc=...]
"""
  if(len(args) == 0):
    print default_message
    return
  cmdline_phil = []
  for arg in args :
    if os.path.isfile(arg):
      if iotbx.pdb.is_pdb_file(arg):
        pdb_phil = libtbx.phil.parse("pdb_file=%s" % os.path.abspath(arg))
        cmdline_phil.append(pdb_phil)
      else:
        try: file_phil = libtbx.phil.parse(file_name=arg)
        except Exception: raise Sorry("Bad parameter file: %s"%arg)
        cmdline_phil.append(file_phil)
    else:
      try: arg_phil = libtbx.phil.parse(arg)
      except Exception: raise Sorry("Bad parameter: %s"%arg)
      cmdline_phil.append(arg_phil)
  working_phil = master_phil.fetch(sources=cmdline_phil)
  params = working_phil.extract()
  pdb_file_name = params.pdb_file
  if (pdb_file_name is None) or (not iotbx.pdb.is_pdb_file(pdb_file_name)) :
    print "A PDB file is required."
    return
  if (params.nproc is None) :
    params.nproc = 1
  pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  pdb_atoms = pdb_hierarchy.atoms()
  pdb_atoms.reset_i_seq()
  #
  xray_structure = pdb_inp.xray_structure_simple()
  return find_tls(
    params         = params,
    pdb_inp        = pdb_inp,
    pdb_hierarchy  = pdb_hierarchy,
    xray_structure = xray_structure,
    out            = out)

def total_score(pdb_hierarchy, sites_cart, u_iso, selection_strings):
  assert sites_cart.size() == u_iso.size()
  target = 0
  for sel_str in selection_strings:
    sel = pdb_hierarchy.atom_selection_cache().selection(
      string = sel_str.replace('"',""))
    assert sel.size() == u_iso.size()
    target += tls_refinery(sites_cart=sites_cart, selection=sel, u_iso=u_iso).f
  return target

def external_tls(pdb_inp, pdb_hierarchy, sites_cart, u_iso, out=None) :
  if (out is None) :
    out = sys.stdout
  pdb_inp_tls = mmtbx.tls.tools.tls_from_pdb_inp(
  remark_3_records = pdb_inp.extract_remark_iii_records(3),
  pdb_hierarchy    = pdb_hierarchy)
  print_statistics.make_header("TLS groups from PDB file header",
    out = out)
  selection_strings = []
  if(len(pdb_inp_tls.tls_params)>0):
    for tp in pdb_inp_tls.tls_params:
      print >> out, "  ", tp.selection_string
      selection_strings.append(tp.selection_string)
  else:
    print >> out, "  ... none found."
  if(len(selection_strings)>0):
    total_target = total_score(
      pdb_hierarchy     = pdb_hierarchy,
      sites_cart        = sites_cart,
      u_iso             = u_iso,
      selection_strings = selection_strings)
    print >> out
    print >> out, "Total target for groups from PDB file header: %10.1f"%total_target

def check_adp(u_iso, step=10, out=None) :
  if (out is None) :
    out = sys.stdout
  min_adp = flex.min(u_iso)
  if(min_adp<=0):
    raise Sorry("Negative or zero isotropic B-factors found in input file. "+
      "Run 'phenix.pdbtools --show-adp-statistics model.pdb' to identify "+
      "the problem atoms.")
  i = 0
  while i < u_iso.size():
    if(i+step < u_iso.size()):
      u_iso_i = u_iso[i:i+step]
    else:
      u_iso_i = u_iso[i:]
    if(u_iso_i.size() >= step//2):
      min_adp = flex.min(u_iso)
      max_adp = flex.max(u_iso)
      if(abs(min_adp-max_adp)<0.1):
        raise Sorry("At least 10 bonded atoms have identical ADPs.")
    i+=step
  return True

def merge_groups_by_connectivity(pdb_hierarchy, xray_structure,
                                 selection_strings=None, selection_arrays=None):
  assert [selection_strings, selection_arrays].count(None)==1
  if(selection_strings is None): selections = selection_arrays
  else:
    selections = []
    for ss in selection_strings:
      sa = pdb_hierarchy.atom_selection_cache().selection(string = ss.replace('"',""))
      selections.append(sa)
  for i_seq, si in enumerate(selections):
    for j_seq, sj in enumerate(selections):
      if(i_seq < j_seq):
        xi = xray_structure.select(si)
        xj = xray_structure.select(sj)
        if(xi.scatterers().size() > xj.scatterers().size()):
          distances = xi.closest_distances(xj.sites_frac(), distance_cutoff=6).smallest_distances
          cnt = ((distances > 0) & (distances < 3)).count(True)
          assert distances.size() == xj.scatterers().size()
          distances = distances.select(distances > 0)
          p = cnt*100./xj.scatterers().size()
          if(p>=1):
            print
            if(selection_strings is not None):
              print sj
              print si
            print i_seq,j_seq, p, flex.min_default(distances,0), flex.mean_default(distances,0)
        else:
          distances = xj.closest_distances(xi.sites_frac(), distance_cutoff=6).smallest_distances
          cnt = ((distances > 0) & (distances < 3)).count(True)
          assert distances.size() == xi.scatterers().size()
          distances = distances.select(distances > 0)
          p = cnt*100./xi.scatterers().size()
          if(p>=1):
            print
            if(selection_strings is not None):
              print sj
              print si
            print i_seq,j_seq, p, flex.min_default(distances,0), flex.mean_default(distances,0)

  #
  print

def find_tls (params,
              pdb_inp,
              pdb_hierarchy,
              xray_structure,
              return_as_list=False,
              out=None) :
  if (out is None) :
    out = sys.stdout
  print_statistics.make_header("Analyzing inputs", out=out)
  if (params.random_seed is None) :
    params.random_seed = flex.get_random_seed()
  random.seed(params.random_seed)
  flex.set_random_seed(params.random_seed)
  xray_structure.convert_to_isotropic()
  sites_cart = xray_structure.sites_cart()
  u_cart = None
  u_iso  = xray_structure.extract_u_iso_or_u_equiv()#*adptbx.u_as_b(1.) # ?
  if(not check_adp(u_iso=u_iso, out=out)): return None
  #
  ssm = mmtbx.secondary_structure.manager(
    pdb_hierarchy                = pdb_hierarchy,
    xray_structure               = xray_structure,
    sec_str_from_pdb_file        = None,
    params                       = None,
    assume_hydrogens_all_missing = None,
    tmp_dir                      = None)
  ssm.find_automatically(log=out)
  alpha_h_selection = ssm.alpha_selection()
  secondary_structure_selection = ssm.alpha_selection() | \
      ssm.beta_selection() | ssm.base_pair_selection()
  if(u_cart is not None):
    assert secondary_structure_selection.size() == u_cart.size()
  else:
    assert secondary_structure_selection.size() == u_iso.size()
  ssm.show_summary(out=out)
  chains_and_residue_selections, secondary_structure_selection = chains_and_atoms(
    pdb_hierarchy                 = pdb_hierarchy,
    secondary_structure_selection = secondary_structure_selection,
    out                           = out)
  chains_and_permutations = []
  chains_and_atom_selection_strings = []
  print_statistics.make_header("Processing chains", out=out)
  if (params.nproc is None) :
    params.nproc = 1
  for crs in chains_and_residue_selections:
    print_statistics.make_sub_header("Processing chain '%s'"%crs[0],
      out=out)
    chain_selection = chain_selection_from_residues(crs[1])
    groups, perms = get_model_partitioning(residues = crs[1],
      secondary_structure_selection = secondary_structure_selection,
      out = out)
    #
    #print
    #selection_arrays = sels_as_selection_arrays(sels = groups)
    #merge_groups_by_connectivity(
    #  pdb_hierarchy     = pdb_hierarchy,
    #  xray_structure    = xray_structure,
    #  selection_arrays  = selection_arrays)
    #assert 0
    #
    if(len(perms)==1):
      print >> out, "  Whole chain is considered as one TLS group."
      chains_and_atom_selection_strings.append([crs[0],[]])
    else:
      print >> out, "  Fitting TLS matrices..."
      dic = {}
      target_best = 1.e+9
      if (params.nproc is Auto) or (params.nproc > 1) :
        process_perms = analyze_permutations(
          groups=groups,
          sites_cart=sites_cart,
          u_cart=u_cart,
          u_iso=u_iso)
        from libtbx import easy_mp
        stdout_and_targets = easy_mp.pool_map(
          processes=params.nproc,
          fixed_func=process_perms,
          args=perms,
          chunksize=100,
          buffer_stdout_stderr=True)
        targets = [ t for so, t in stdout_and_targets ]
        for (perm, target) in zip(perms, targets) :
          dic.setdefault(len(perm), []).append([target,perm])
      else :
        for i_perm, perm in enumerate(perms):
          if i_perm%500==0:
            print >> out, "    ...perm %d of %d"%(i_perm, len(perms))
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
      print >> out, "    Best fits:"
      print >> out, "      No. of         Targets"
      print >> out, "      groups   best   rand.pick diff.  score permutation"
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
        print >> out, "         %3d   %6.1f   %6.1f %6.1f %6.1f"%(
          k,t_best, r, r-t_best, score), perm_best
        if(score > score_best):
          score_best = score
          perm_choice = perm_best[:]
      #
      chains_and_permutations.append([crs[0],perm_choice])
      chains_and_atom_selection_strings.append([crs[0],
        permutations_as_atom_selection_string(groups, perm_choice)])
      #
  if (pdb_inp is not None) :
    external_tls_selections = external_tls(
      pdb_inp       = pdb_inp,
      pdb_hierarchy = pdb_hierarchy,
      sites_cart    = sites_cart,
      u_iso         = u_iso,
      out           = out)
  print_statistics.make_header("SUMMARY", out=out)
  #print "Optimal TLS groups:"
  #for chain_and_permutation in chains_and_permutations:
  #  print chain_and_permutation
  #print
  print >> out, "TLS atom selections for phenix.refine:"
  groups_out = cStringIO.StringIO()
  selection_strings = []
  print >> groups_out, "refinement.refine.adp {"
  for r in chains_and_atom_selection_strings:
    prefix = "chain '%s'"%r[0]
    if(len(r[1])>0 and len(r[1:])>0):
      prefix += " and "
      for r_ in r[1:]:
        for r__ in r_:
          if(len(r__)>0):
            group_selection = prefix+"(%s)"%r__
            print >> groups_out, "  tls = \"%s\"" % group_selection
            selection_strings.append("%s" % group_selection)
    else:
      print >> groups_out, "  tls = \"%s\"" % prefix
      selection_strings.append("%s" % prefix)
  print >> groups_out, "}"
  print >> out, groups_out.getvalue()
  print >> out
  #XXX
  if 0:
    merge_groups_by_connectivity(
      pdb_hierarchy     = pdb_hierarchy,
      xray_structure    = xray_structure,
      selection_strings = selection_strings)
  #XXX
  if(len(selection_strings)>0):
    total_target = total_score(
      pdb_hierarchy     = pdb_hierarchy,
      sites_cart        = sites_cart,
      u_iso             = u_iso,
      selection_strings = selection_strings)
    print >> out, "Overall best total target for automatically found groups: %10.1f"%total_target
    print >> out
  if (return_as_list) :
    return selection_strings
  else :
    return groups_out.getvalue()

# XXX wrapper for running in Phenix GUI
class _run_find_tls (object) :
  def __init__ (self, params, pdb_hierarchy, xray_structure) :
    self.params = params
    self.pdb_hierarchy = pdb_hierarchy
    self.xray_structure = xray_structure

  def __call__ (self, *args, **kwds) :
    return find_tls(
      params=self.params,
      pdb_inp=None,
      pdb_hierarchy=self.pdb_hierarchy,
      xray_structure=self.xray_structure)

if (__name__ == "__main__"):
  t0 = time.time()
  run(args=sys.argv[1:])
  print "Time: %10.3f"%(time.time()-t0)
