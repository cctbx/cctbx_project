from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
import math, sys, os
from libtbx import group_args
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
import libtbx.load_env
from libtbx import easy_pickle
from cctbx import uctbx
from mmtbx.utils import run_reduce_with_timeout

import numpy as np # XXX See if I can avoid it!

from mmtbx.secondary_structure import manager as ss_manager
from mmtbx.secondary_structure import sec_str_master_phil_str

def pymol_water_bonds(model):
  def pymol_atom_selection(atom):
    ag = atom.parent()
    rg = ag.parent()
    chain = rg.parent()
    one = "chain %s and resi %s and name %s and alt '%s'" % (
      chain.id, rg.resseq, atom.name, ag.altloc)
    return one
  #
  ph=model.get_hierarchy()
  outl = ''
  for ag in ph.atom_groups():
    if ag.resname!='HOH': continue
    oxygen = None
    for atom in ag.atoms():
      if atom.element.strip()=='O': oxygen=atom
      elif atom.element_is_hydrogen:
        outl += 'bond %s, %s\n' % (pymol_atom_selection(oxygen),
                                   pymol_atom_selection(atom))
  return outl

mcss = " or ".join(
    ["name %s"%i.strip() for i in iotbx.pdb.protein_atom_names_backbone])

def get_pair_generator(crystal_symmetry, buffer_thickness, sites_cart):
  sst = crystal_symmetry.special_position_settings().site_symmetry_table(
    sites_cart = sites_cart)
  from cctbx import crystal
  conn_asu_mappings = crystal_symmetry.special_position_settings().\
    asu_mappings(buffer_thickness = buffer_thickness)
  conn_asu_mappings.process_sites_cart(
    original_sites      = sites_cart,
    site_symmetry_table = sst)
  conn_pair_asu_table = crystal.pair_asu_table(asu_mappings = conn_asu_mappings)
  conn_pair_asu_table.add_all_pairs(distance_cutoff = buffer_thickness)
  pair_generator = crystal.neighbors_fast_pair_generator(
    conn_asu_mappings, distance_cutoff = buffer_thickness)
  return group_args(
    pair_generator    = pair_generator,
    conn_asu_mappings = conn_asu_mappings)

def apply_symop_to_copy(atom, rt_mx_ji, fm, om):
  atom = atom.detached_copy()
  t1 = fm*flex.vec3_double([atom.xyz])
  t2 = rt_mx_ji*t1[0]
  t3 = om*flex.vec3_double([t2])
  atom.set_xyz(t3[0])
  return atom

def make_atom_id(atom, index):
  return group_args(
    id_str = atom.id_str().replace("pdb=",""),
    index  = index,
    name   = atom.name,
    b      = atom.b,
    occ    = atom.occ,
    chain  = atom.parent().parent().parent().id,
    resseq = atom.parent().parent().resseq,
    altloc = atom.parent().altloc)

def get_stats(data, min_data_size=10):
  if(data.size()<min_data_size): return None
  mean=data.min_max_mean().mean
  sd=data.standard_deviation_of_the_sample()
  assert data.size(), 'no data - may mean no Hydrogen atoms'
  x=data-mean
  skew=kurtosis=None
  if sd:
    skew=(x**3).min_max_mean().mean/sd**3
    kurtosis=(x**4).min_max_mean().mean/sd**4
  return group_args(mean=mean, sd=sd, skew=skew, kurtosis=kurtosis)

master_phil_str = '''
hbond {
  show_hbonds = True
    .type = bool
  output_pymol_file = False
    .type = bool
    .short_caption = Output PyMOL files for visualisation
  output_restraint_file = False
    .type = bool
    .short_caption = Output geometry restraints edits file
  output_stats_pdf = False
    .type = bool
    .short_caption = Output stats in PDF
  add_hydrogens_if_absent = False
    .type = bool
    .short_caption = Add hydrogens if they are not present in the model(s)
  min_data_size = 10
    .type = int
    .style = hidden
  dot_size = 100
    .type = int
    .help = Dot size
}
'''

def show_histogram(data, n_slots, data_min, data_max, log=sys.stdout):
  from cctbx.array_family import flex
  h_data = flex.double()
  hm = flex.histogram(
    data=data, n_slots=n_slots, data_min=data_min, data_max=data_max)
  lc_1 = hm.data_min()
  s_1 = enumerate(hm.slots())
  for (i_1,n_1) in s_1:
    hc_1 = hm.data_min() + hm.slot_width() * (i_1+1)
    #print >> log, "%10.5f - %-10.5f : %d" % (lc_1, hc_1, n_1)
    #print >> log, "%10.2f : %d" % ((lc_1+hc_1)/2, n_1)
    print ("%10.2f : %10.4f" % ((lc_1+hc_1)/2, n_1*100./data.size()), file=log)
    lc_1 = hc_1
  return h_data

# XXX FIND A BETTER PLACE
def get_ss_selections(hierarchy, filter_short=True):
  def get_counts(hierarchy, h_sel, s_sel):
    sh_sel = h_sel | s_sel
    n   = hierarchy.overall_counts().n_residues
    nh  = hierarchy.select(h_sel ).overall_counts().n_residues
    ns  = hierarchy.select(s_sel ).overall_counts().n_residues
    nhs = hierarchy.select(sh_sel).overall_counts().n_residues
    return group_args(
      h  = int(round(nh *100./n,0)),
      s  = int(round(ns *100./n,0)),
      hs = int(round(nhs*100./n,0)))
  def one(hierarchy, method):
    sec_str_master_phil = iotbx.phil.parse(sec_str_master_phil_str)
    params = sec_str_master_phil.fetch().extract()
    params.secondary_structure.protein.search_method = method
    asc = hierarchy.atom_selection_cache()
    ssm = ss_manager(
      hierarchy,
      atom_selection_cache=asc,
      geometry_restraints_manager=None,
      sec_str_from_pdb_file=None,
      params=params.secondary_structure,
      was_initialized=False,
      verbose=-1,
      log=null_out())
    filtered_ann = ssm.actual_sec_str.deep_copy()
    if(filter_short):
      filtered_ann.remove_short_annotations(
        helix_min_len=4, sheet_min_len=4, keep_one_stranded_sheets=True)
    mc_sel = asc.selection(mcss)
    h_sel  = asc.selection(filtered_ann.overall_helices_selection())
    s_sel  = asc.selection(filtered_ann.overall_sheets_selection())
    h_sel  = h_sel & mc_sel
    s_sel  = s_sel & mc_sel
    ss_counts = get_counts(hierarchy=hierarchy, h_sel=h_sel, s_sel=s_sel)
    return group_args(h_sel = h_sel, s_sel = s_sel, counts = ss_counts)
  ksdssp, from_ca, both = None, None, None
  # from_ca
  from_ca = one(hierarchy=hierarchy, method="from_ca")
  # ksdssp
  try:
    ksdssp  = one(hierarchy=hierarchy, method="ksdssp")
  except KeyboardInterrupt: raise
  except: pass # intentional
               # really don't know what else to do here!
  # both
  if([ksdssp, from_ca].count(None)==0):
    h_sel = ksdssp.h_sel & from_ca.h_sel
    s_sel = ksdssp.s_sel & from_ca.s_sel
    ss_counts = get_counts(hierarchy=hierarchy, h_sel=h_sel, s_sel=s_sel)
    both = group_args(h_sel = h_sel, s_sel = s_sel, counts = ss_counts)
  #
  return group_args(ksdssp=ksdssp, from_ca=from_ca, both=both)

def stats(model, prefix, output_stats_pdf, no_ticks=True):
  # Get rid of H, multi-model, no-protein and single-atom residue models
  if(model.percent_of_single_atom_residues()>20):
    return None
  sel = model.selection(string = "protein")
  if(sel.count(True)==0):
    return None
  ssr = "protein and not (element H or element D or resname UNX or resname UNK or resname UNL)"
  sel = model.selection(string = ssr)
  model = model.select(sel)
  if(len(model.get_hierarchy().models())>1):
    return None
  # Add H; this looses CRYST1 !
  rr = run_reduce_with_timeout(
    stdin_lines = model.get_hierarchy().as_pdb_string().splitlines(),
    file_name   = None,
    parameters  = "-oh -his -flip -keep -allalt -pen9999 -",
    override_auto_timeout_with=None)
  # Create model; this is a single-model pure protein with new H added
  pdb_inp = iotbx.pdb.input(source_info = None, lines = rr.stdout_lines)
  model = mmtbx.model.manager(
    model_input      = pdb_inp,
    log              = null_out())
  if(model.crystal_symmetry() is None):
    box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
      sites_cart   = model.get_sites_cart(),
      buffer_layer = 5)
    model.set_sites_cart(box.sites_cart)
    model._crystal_symmetry = box.crystal_symmetry()
  model.process(make_restraints = True)
  #
  N = 10
  #
  import cctbx.geometry_restraints.process_nonbonded_proxies as pnp
  h_bond_params = pnp.h_bond()
  h_bond_params.a_DHA_cutoff=90
  #
  SS = get_ss_selections(hierarchy=model.get_hierarchy())
  HB_all = find(model = model.select(flex.bool(model.size(), True)),
    h_bond_params=h_bond_params)
  HB_alpha = find(model = model.select(SS.both.h_sel), h_bond_params=h_bond_params)
  HB_beta = find(model = model.select(SS.both.s_sel), h_bond_params=h_bond_params)
  #
  if(output_stats_pdf):
    result_dict = {}
    result_dict["all"]   = HB_all.get_params_as_arrays(replace_with_empty_threshold=N)
    result_dict["alpha"] = HB_alpha.get_params_as_arrays(replace_with_empty_threshold=N)
    result_dict["beta"]  = HB_beta.get_params_as_arrays(replace_with_empty_threshold=N)
    # Load histograms for reference high-resolution d_HA and a_DHA
    pkl_fn = libtbx.env.find_in_repositories(
      relative_path="mmtbx")+"/nci/d_HA_and_a_DHA_high_res_p3.pkl"
    assert os.path.isfile(pkl_fn)
    ref = easy_pickle.load(pkl_fn)
    #
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(10,10))
    kwargs = dict(histtype='bar', bins=20, range=[1.6,3.0], alpha=.8)
    for j, it in enumerate([["alpha",1], ["beta",3], ["all",5]]):
      key, i = it
      ax = plt.subplot(int("32%d"%i))
      if(no_ticks):
        #ax.set_xticks([])
        ax.set_yticks([])
      if(j in [0,1]):
        ax.tick_params(bottom=False)
        ax.set_xticklabels([])
      ax.tick_params(axis="x", labelsize=12)
      ax.tick_params(axis="y", labelsize=12, left=False, pad=-2)
      ax.text(0.98,0.92,key, size=12, horizontalalignment='right',
        transform=ax.transAxes)
      HB = result_dict[key]
      if HB is None: continue
      w1 = np.ones_like(HB.d_HA)/HB.d_HA.size()
      ax.hist(HB.d_HA, color="orangered", weights=w1, rwidth=0.3, **kwargs)
      #
      start, end1, end2 = 0, max(ref.distances[key].vals), \
        round(max(ref.distances[key].vals),2)
      if(not no_ticks):
        plt.yticks([0.01,end1], ["0", end2], visible=True, rotation="horizontal")

      if  (key=="alpha"): plt.ylim(0, end2+0.02)
      elif(key=="beta"):  plt.ylim(0, end2+0.02)
      elif(key=="all"):  plt.ylim(0, end2+0.02)
      else: assert 0
      #
      if(j==0): ax.set_title("Distance", size=15)
      bins = list(flex.double(ref.distances[key].bins))
      ax.bar(bins, ref.distances[key].vals, alpha=.3, width=0.07)
    #
    kwargs = dict(histtype='bar', bins=20, range=[90,180], alpha=.8)
    for j, it in enumerate([["alpha",2], ["beta",4], ["all",6]]):
      key, i = it
      ax = plt.subplot(int("32%d"%i))
      if(j in [0,1]):
        ax.tick_params(bottom=False)
        ax.set_xticklabels([])
      if(no_ticks):
        #ax.set_xticks([])
        ax.set_yticks([])
      ax.tick_params(axis="x", labelsize=12)
      ax.tick_params(axis="y", labelsize=12, left=False, pad=-2)
      ax.text(0.98,0.92,key, size=12, horizontalalignment='right',
        transform=ax.transAxes)

      ax.text(0.98,0.92,key, size=12, horizontalalignment='right',
        transform=ax.transAxes)
      #if(j in [0,1]): ax.plot_params(bottom=False)
      HB = result_dict[key]
      if HB is None: continue
      w1 = np.ones_like(HB.a_DHA)/HB.a_DHA.size()
      ax.hist(HB.a_DHA, color="orangered", weights=w1, rwidth=0.3, **kwargs)
      #
      start, end1, end2 = 0, max(ref.angles[key].vals), \
        round(max(ref.angles[key].vals),2)
      if(not no_ticks):
        plt.yticks([0.01,end1], ["0", end2], visible=True, rotation="horizontal")

      if  (key=="alpha"): plt.ylim(0, end2+0.02)
      elif(key=="beta"):  plt.ylim(0, end2+0.02)
      elif(key=="all"):  plt.ylim(0, end2+0.02)
      else: assert 0
      #
      if(j==0): ax.set_title("Angle", size=15)
      ax.bar(ref.angles[key].bins, ref.angles[key].vals, width=4.5, alpha=.3)
    plt.subplots_adjust(wspace=0.12, hspace=0.025)
    if(no_ticks):
      plt.subplots_adjust(wspace=0.025, hspace=0.025)
    #fig.savefig("%s.png"%prefix, dpi=1000)
    fig.savefig("%s.pdf"%prefix)
    # end save as PDF part
  return group_args(all = HB_all, alpha = HB_alpha, beta = HB_beta)


def precheck(atoms, i, j, Hs, As, Ds, fsc0, tolerate_altloc=False):
  """
  Check if two atoms are potential H bond partners, based on element and altloc
  """
  ei, ej = atoms[i].element, atoms[j].element
  altloc_i = atoms[i].parent().altloc
  altloc_j = atoms[j].parent().altloc
  resseq_i = atoms[i].parent().parent().resseq
  resseq_j = atoms[j].parent().parent().resseq
  one_is_Hs = ei in Hs or ej in Hs
  other_is_acceptor = ei in As or ej in As
  if tolerate_altloc:
    is_candidate = one_is_Hs and other_is_acceptor and resseq_i != resseq_j
  else:
    is_candidate = one_is_Hs and other_is_acceptor and \
      altloc_i == altloc_j and resseq_i != resseq_j
  if(ei in Hs):
    bound_to_h = fsc0[i]
    if(not bound_to_h): # exclude 'lone' H
      is_candidate = False
    elif(atoms[bound_to_h[0]].element not in Ds): # Use only first atom bound to H
      is_candidate = False
  if(ej in Hs):
    bound_to_h = fsc0[j]
    if(not bound_to_h):
      is_candidate = False
    elif(atoms[bound_to_h[0]].element not in Ds):
      is_candidate = False
  return is_candidate

def get_D_H_A_Y(i, j, Hs, fsc0, rt_mx_ji, fm, om, atoms):
  """
  Get atom objects for donor and acceptor atoms
  Apply symmetry op if necessary, so that correct geometry can be calculated
  """
  Y = []
  if(atoms[i].element in Hs):
    H = atoms[i]
    D = atoms[fsc0[i][0]]
    A = atoms[j]
    Y_iseqs = fsc0[j]
    if(len(Y_iseqs)>0):
      Y = [atoms[k] for k in fsc0[j]]
    atom_H = make_atom_id(atom = H, index = i)
    atom_A = make_atom_id(atom = A, index = j)
    atom_D = make_atom_id(atom = D, index = D.i_seq)
    if(rt_mx_ji is not None and str(rt_mx_ji) != "x,y,z"):
      A = apply_symop_to_copy(A, rt_mx_ji, fm, om)
      if(len(Y_iseqs)>0):
        Y = [apply_symop_to_copy(y, rt_mx_ji, fm, om) for y in Y]
  if(atoms[j].element in Hs):
    H = atoms[j]
    D = atoms[fsc0[j][0]]
    A = atoms[i]
    Y_iseqs = fsc0[i]
    if(len(Y_iseqs)>0):
      Y = [atoms[k] for k in fsc0[i]]
    atom_A = make_atom_id(atom = A, index = i)
    atom_H = make_atom_id(atom = H, index = j)
    atom_D = make_atom_id(atom = D, index = D.i_seq)
    if(rt_mx_ji is not None and str(rt_mx_ji) != "x,y,z"):
      H = apply_symop_to_copy(H, rt_mx_ji, fm, om)
      D = apply_symop_to_copy(D, rt_mx_ji, fm, om)
  return D, H, A, Y, atom_A, atom_H, atom_D

class find(object):
  def __init__(self,
        model,
        h_bond_params  = None,
        protein_only   = False,
        pair_proxies   = None):
    assert not protein_only
    if(h_bond_params is None):
      import cctbx.geometry_restraints.process_nonbonded_proxies as pnp
      h_bond_params = pnp.h_bond()
    Hs           = h_bond_params.Hs
    As           = h_bond_params.As
    Ds           = h_bond_params.Ds
    d_HA_cutoff  = h_bond_params.d_HA_cutoff
    d_DA_cutoff  = h_bond_params.d_DA_cutoff
    a_DHA_cutoff = h_bond_params.a_DHA_cutoff
    a_YAH_cutoff = h_bond_params.a_YAH_cutoff
    #
    self.result = []
    self.model = model
    self.pair_proxies = pair_proxies
    self.external_proxies = False
    if(self.pair_proxies is not None):
      self.external_proxies = True
    atoms = self.model.get_hierarchy().atoms()
    geometry = self.model.get_restraints_manager()
    fsc0 = geometry.geometry.shell_sym_tables[0].full_simple_connectivity()
    bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
      sites_cart = self.model.get_sites_cart())
    sites_cart = self.model.get_sites_cart()
    crystal_symmetry = self.model.crystal_symmetry()
    fm = crystal_symmetry.unit_cell().fractionalization_matrix()
    om = crystal_symmetry.unit_cell().orthogonalization_matrix()
    pg = get_pair_generator(
      crystal_symmetry = crystal_symmetry,
      buffer_thickness = d_HA_cutoff[1],
      sites_cart       = sites_cart)
    get_class = iotbx.pdb.common_residue_names_get_class
    # find proxies if not provided
    if(self.pair_proxies is None):
      pp = []
      self.pair_proxies = []
      pp = [p for p in pg.pair_generator]
    else:
      pp = self.pair_proxies
    # now loop over proxies
    for p in pp:
      i, j = p.i_seq, p.j_seq
      if(self.external_proxies): # making sure proxies point to same atoms
        a_i = make_atom_id(atom = atoms[i], index = i).id_str
        a_j = make_atom_id(atom = atoms[j], index = j).id_str
        assert a_i == p.atom_A.id_str, [a_i, p.atom_A.id_str]
        assert a_j == p.atom_H.id_str, [a_j, p.atom_H.id_str]
      # presecreen candidates
      ei, ej = atoms[i].element, atoms[j].element
      is_candidate = precheck(
        atoms = atoms,
        i = i,
        j = j,
        Hs = Hs,
        As = As,
        Ds = Ds,
        fsc0 = fsc0)
      if(protein_only):
        for it in [i,j]:
          resname = atoms[it].parent().resname
          is_candidate &= get_class(name=resname) == "common_amino_acid"
      if(not is_candidate): continue
      # pre-screen candidates end
      # symop tp map onto symmetry related
      rt_mx_ji = None
      if(not self.external_proxies):
        rt_mx_i = pg.conn_asu_mappings.get_rt_mx_i(p)
        rt_mx_j = pg.conn_asu_mappings.get_rt_mx_j(p)
        rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)
      else:
        rt_mx_ji = p.rt_mx_ji
      #
      D, H, A, Y, atom_A, atom_H, atom_D = get_D_H_A_Y(
        i        = i,
        j        = j,
        Hs       = Hs,
        fsc0     = fsc0,
        rt_mx_ji = rt_mx_ji,
        fm       = fm,
        om       = om,
        atoms    = atoms)
      if(len(Y) == 0): continue # don't use 'lone' acceptors
      #
      d_DA = D.distance(A)
      if(not self.external_proxies):
        if(d_DA < d_DA_cutoff[0] or d_DA > d_DA_cutoff[1]):
          continue
      #
      d_HA = A.distance(H)
      if(not self.external_proxies):
        assert d_HA <= d_HA_cutoff[1]
        assert approx_equal(math.sqrt(p.dist_sq), d_HA, 1.e-3)
        if(d_HA < d_HA_cutoff[0]): continue
#      assert H.distance(D) < 1.15, [H.distance(D), H.name, D.name]
      # filter by a_DHA
      a_DHA = H.angle(A, D, deg=True)
      if(not self.external_proxies):
        if(a_DHA < a_DHA_cutoff): continue
      # filter by a_YAH
      a_YAH = []
      if(len(Y)>0):
        for Y_ in Y:
          a_YAH_ = A.angle(Y_, H, deg=True)
          a_YAH.append(a_YAH_)
      if(not self.external_proxies):
        flags = []
        for a_YAH_ in a_YAH:
          flags.append(
            not (a_YAH_ >= a_YAH_cutoff[0] and a_YAH_ <= a_YAH_cutoff[1]))
        flags = list(set(flags))
        if(len(flags)>1 or (len(flags)==1 and flags[0])): continue
      #
      assert approx_equal(d_HA, H.distance(A), 1.e-3)
      #a_YAD = []
      #if(len(Y)>0):
      #  for Y_ in Y:
      #    a_YAD_ = A.angle(Y_, D, deg=True)
      #    a_YAD.append(a_YAD_)
      self.result.append(group_args(
        i       = i,
        j       = j,
        atom_H  = atom_H,
        atom_A  = atom_A,
        atom_D  = atom_D,
        symop   = rt_mx_ji,
        d_HA    = d_HA,
        a_DHA   = a_DHA,
        a_YAH   = a_YAH,
        #a_YAD   = a_YAD,
        d_AD    = A.distance(D)
      ))
      if(not self.external_proxies):
        proxy_custom = group_args(i_seq = i, j_seq = j, rt_mx_ji = rt_mx_ji,
          atom_H = atom_H, atom_A = atom_A)
        self.pair_proxies.append(proxy_custom)

  def get_params_as_arrays(self, b=None, occ=None,
                           replace_with_empty_threshold=None):
    d_HA  = flex.double()
    a_DHA = flex.double()
    a_YAH = flex.double()
    for r in self.result:
      if(b   is not None and r.atom_H.b>b): continue
      if(b   is not None and r.atom_A.b>b): continue
      if(occ is not None and r.atom_H.occ<occ): continue
      if(occ is not None and r.atom_A.occ<occ): continue
      d_HA .append(r.d_HA )
      a_DHA.append(r.a_DHA)
      if(len(r.a_YAH)>0):
        a_YAH.extend(flex.double(r.a_YAH))
    if(replace_with_empty_threshold is not None and
       d_HA.size()<replace_with_empty_threshold):
      d_HA  = flex.double()
      a_DHA = flex.double()
      a_YAH = flex.double()
    return group_args(d_HA=d_HA, a_DHA=a_DHA, a_YAH=a_YAH)

  def get_counts(self, b=None, occ=None, filter_id_str=None, min_data_size=10):
    theta_1 = flex.double()
    theta_2 = flex.double()
    d_HA    = flex.double()
    n_sym = 0
    for r in self.result:
      if(str(r.symop) != "x,y,z"):
        n_sym += 1
      if(b   is not None and r.atom_H.b>b): continue
      if(b   is not None and r.atom_A.b>b): continue
      if(occ is not None and r.atom_H.occ<occ): continue
      if(occ is not None and r.atom_A.occ<occ): continue
      if(filter_id_str is not None and
         (r.atom_A.id_str.find(filter_id_str)==-1 and r.atom_H.id_str.find(filter_id_str)==-1)
         ):
         continue
      theta_1.append(r.a_DHA)
      theta_2.extend(flex.double(r.a_YAH))
      d_HA   .append(r.d_HA)
    n_filter=len(theta_1)
    bpr=float(len(self.result))/\
      len(list(self.model.get_hierarchy().residue_groups()))
    if len(theta_1)==0: return None
    theta_1 = get_stats(theta_1, min_data_size=min_data_size)
    theta_2 = get_stats(theta_2, min_data_size=min_data_size)
    d_HA    = get_stats(d_HA, min_data_size=min_data_size)
    if([theta_1, theta_2, d_HA].count(None)>0): return None
    return group_args(
      theta_1 = theta_1,
      theta_2 = theta_2,
      d_HA    = d_HA,
      n       = len(self.result),
      n_filter= n_filter,
      n_sym   = n_sym,
      bpr     = bpr)

  def show_summary(self, log = sys.stdout):
    def printit(o, prefix):
      fmt="%8s %7.3f %7.3f %7.3f %7.3f"
      print(fmt%(prefix, o.mean, o.sd, o.skew, o.kurtosis), file=log)
    c = self.get_counts()
    if(c is None): return
    if(c.theta_1 is None): return
    print("Total:       %d"%c.n,     file=log)
    print("Symmetry:    %d"%c.n_sym, file=log)
    print("Per residue: %7.4f"%c.bpr,   file=log)
    min_stat_limit=50
    if c.n<min_stat_limit:
      print('Statistics not displayed due to low number of H bonds found.', file=log)
    else:
      print("            Mean      SD    Skew   Kurtosis",   file=log)
      printit(c.theta_1, "theta_1:")
      printit(c.theta_2, "theta_2:")
      printit(c.d_HA,    "d_HA:")

  def show(self, log = sys.stdout, sym_only=False):
    for r in self.result:
      ids_i = r.atom_H.id_str
      ids_j = r.atom_A.id_str
      if(sym_only):
        if(str(r.symop)=="x,y,z"): continue
      print("%4d %4d"%(r.i,r.j), "%s<>%s"%(ids_i, ids_j), \
        "d_HA=%5.3f"%r.d_HA, "d_AD=%5.3f"%r.d_AD, "a_DHA=%7.3f"%r.a_DHA, \
        "symop: %s"%str(r.symop), " ".join(["a_YAH=%d"%i for i in r.a_YAH]),
        file=log)

  def as_pymol(self, prefix="hbonds_pymol"):
    pdb_file_name = "%s.pdb"%prefix
    with open(pdb_file_name, "w") as of:
      print(self.model.model_as_pdb(), file=of)
    water_bonds = pymol_water_bonds(self.model)
    with open("%s.pml"%prefix, "w") as of:
      print("load", "/".join([os.getcwd(), pdb_file_name]), file=of)
      for r in self.result:
        if(str(r.symop) != "x,y,z"): continue
        ai = r.atom_H
        aj = r.atom_A
        one = "chain %s and resi %s and name %s and alt '%s'"%(
          ai.chain, ai.resseq, ai.name, ai.altloc)
        two = "chain %s and resi %s and name %s and alt '%s'"%(
          aj.chain, aj.resseq, aj.name, aj.altloc)
        print("dist %s, %s"%(one, two), file=of)
      if water_bonds:
        print(water_bonds, file=of)

  def as_restraints(self, file_name="hbond.eff", distance_ideal=None, sigma_dist=0.1,
       angle_ideal = None, sigma_angle=2, use_actual=True):
    f = "chain %s and resseq %s and name %s"
    with open(file_name, "w") as of:
      print("geometry_restraints.edits {", file=of)
      for r in self.result:
        h = f%(r.atom_H.chain, r.atom_H.resseq, r.atom_H.name)
        a = f%(r.atom_A.chain, r.atom_A.resseq, r.atom_A.name)
        d = f%(r.atom_D.chain, r.atom_D.resseq, r.atom_D.name)
        if(not use_actual):
          if(r.d_HA<2.5): dt = 2.05
          else:           dt = 2.8
          if(r.a_DHA<130): at = 115
          else:            at = 160
        else:
          dt = r.d_HA
          at = r.a_DHA
        dis = """    bond {
          atom_selection_1 = %s
          atom_selection_2 = %s
          symmetry_operation = %s
          distance_ideal = %f
          sigma = 0.05
         }"""%(h,a,str(r.symop),dt)
        if(str(r.symop)!="x,y,z"): continue
        ang = """    angle {
          atom_selection_1 = %s
          atom_selection_2 = %s
          atom_selection_3 = %s
          angle_ideal = %f
          sigma = 5
          }"""%(a,h,d,at)
        print(dis, file=of)
        print(ang, file=of)
      print("}", file=of)
