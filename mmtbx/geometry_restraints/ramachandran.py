from __future__ import division
import libtbx.load_env
import iotbx.phil
from libtbx.utils import Sorry
from libtbx import adopt_init_args
import sys
import os
import boost.python
from scitbx.array_family import flex

ext = boost.python.import_ext("mmtbx_ramachandran_restraints_ext")
from mmtbx_ramachandran_restraints_ext import lookup_table, \
    ramachandran_residual_sum

# rosetta_adaptbx is not considered as working anymore.
# Probably should be removed at some point.
if libtbx.env.has_module("rosetta_adaptbx") :
  potential_phil = """\
   rama_potential = *oldfield emsley rosetta
    .type = choice(multi=False)
    .short_caption = Ramachandran potential
    .caption = Oldfield Coot Rosetta"""
else :
  potential_phil = """\
   rama_potential = *oldfield emsley
    .type = choice(multi=False)
    .short_caption = Ramachandran potential
    .caption = Oldfield Coot"""

master_phil = iotbx.phil.parse("""
  rama_weight = 1.0
    .type = float
    .short_caption = Ramachandran gradients weight
    .expert_level = 1
  scale_allowed = 1.0
    .type = float
    .short_caption = Rescale allowed region pseudo-energy by
  %s
  oldfield
    .short_caption = Oldfield potential parameters
    .style = box auto_align
  {
    esd = 10.0
      .type = float
      .expert_level = 2
      .short_caption = E.S.D.
    weight_scale = 1.0
      .type = float
      .expert_level = 2
    dist_weight_max = 10.0
      .type = float
      .expert_level = 2
    weight = None
      .type = float
      .expert_level = 2
  }
  rama_selection = None
    .type = atom_selection
    .short_caption = Atom selection for Ramachandran restraints
    .expert_level = 1
  rama_exclude_sec_str = False
    .type = bool
    .expert_level = 1
    .short_caption = Exclude secondary structure from Ramachandran restraints
""" % potential_phil)

def is_proxy_present(proxies, n_seq, proxy):
  p_iseqs = list(proxy.get_i_seqs())
  ps = proxies.proxy_select(n_seq=n_seq,
      iselection=flex.size_t(p_iseqs))
  return ps.size() > 0

class ramachandran_manager(object):
  def __init__ (self, pdb_hierarchy, atom_selection=None, params=None,
      log=sys.stdout, proxies=None, tables=None, initialize=True):
    assert pdb_hierarchy is not None
    assert not pdb_hierarchy.atoms().extract_i_seq().all_eq(0), ""+\
        "Probably all atoms have i_seq = 0 which is wrong"
    adopt_init_args(self, locals())
    self.bool_atom_selection = None
    if self.atom_selection is None:
      self.bool_atom_selection = flex.bool(pdb_hierarchy.atoms().size(), True)
    else:
      cache = pdb_hierarchy.atom_selection_cache()
      self.bool_atom_selection = cache.selection(atom_selection)
    if params is None:
      params = master_phil.fetch().extract()
    self.params = params
    if initialize:
      if(self.params.rama_potential == "oldfield"):
        self.tables = ramachandran_plot_data()
      elif (self.params.rama_potential == "rosetta") :
        import rosetta_adaptbx
        rosetta_adaptbx.init()
        self.tables = rosetta_adaptbx.ext.ramachandran()
      else :
        self.tables = load_tables(params)
      # get proxies
      self.extract_proxies()
    else:
      assert proxies is not None

  def proxy_select(self, n_seq, iselection):
    result_proxies = self.proxies.proxy_select(n_seq, iselection)
    return ramachandran_manager(
        pdb_hierarchy=self.pdb_hierarchy,
        atom_selection=self.atom_selection,
        params=self.params,
        log=self.log,
        proxies=result_proxies,
        tables=self.tables,
        initialize=False)

  def extract_proxies(self):
    self.proxies = ext.shared_phi_psi_proxy()
    from mmtbx.conformation_dependent_library import generate_protein_threes
    selected_h = self.pdb_hierarchy.select(self.bool_atom_selection)
    n_seq = flex.max(selected_h.atoms().extract_i_seq())
    for three in generate_protein_threes(
        hierarchy=selected_h,
        geometry=None):
      phi_atoms, psi_atoms = three.get_phi_psi_atoms()
      rama_key = three.get_ramalyse_key()
      i_seqs = [atom.i_seq for atom in phi_atoms] + [psi_atoms[-1].i_seq]
      resnames = three.get_resnames()
      r_name = resnames[1]
      assert rama_key in ["general", "glycine", "cis-proline", "trans-proline",
                 "pre-proline", "isoleucine or valine"]
      proxy = ext.phi_psi_proxy(
          residue_name=r_name,
          residue_type=rama_key,
          i_seqs=i_seqs)
      if not is_proxy_present(self.proxies, n_seq, proxy):
        self.proxies.append(proxy)
    print >> self.log, ""
    print >> self.log, "  %d Ramachandran restraints generated." % (
        self.get_n_proxies())

  def target_and_gradients (self,
      unit_cell,
      sites_cart,
      gradient_array=None,
      residuals_array=None) :
    assert self.proxies is not None
    if(gradient_array is None) :
      gradient_array = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
    if residuals_array is None:
      residuals_array = flex.double(self.proxies.size())
    target = 0
    if(self.params.rama_potential == "oldfield"):
      op = self.params.oldfield
      w = op.weight
      if w is None:
        w = -1
      res = ramachandran_residual_sum(
          sites_cart=sites_cart,
          proxies=self.proxies,
          gradient_array=gradient_array,
          general_table=self.tables.general,
          gly_table=self.tables.gly,
          cispro_table=self.tables.cispro,
          transpro_table=self.tables.transpro,
          prepro_table=self.tables.prepro,
          ileval_table=self.tables.ileval,
          weights=(w, op.esd, op.dist_weight_max, op.weight_scale),
          residuals_array=residuals_array)
      return res
    elif (self.params.rama_potential == "rosetta") :
      # Moving these cycles to C++ part would speed them up only up to 10%
      for i, proxy in enumerate(self.proxies):
        residuals_array[i] = self.tables.residue_target_and_gradients(
          gradient_array=gradient_array,
          sites_cart=sites_cart,
          proxy=proxy,
          weight=self.params.rama_weight)
    else: # emsley
      assert (self.params.rama_weight >= 0.0)
      # bad hack to keep emsley potential in working(?) condition after
      # changing from rama500 to rama8000
      new_to_old_conversion = {"general":"ala", "glycine":"gly",
          "cis-proline":"pro", "trans-proline":"pro", "pre-proline":"prepro",
          "isoleucine or valine":"ala"}
      # Moving these cycles to C++ part would speed them up only up to 10%
      for i, proxy in enumerate(self.proxies):
        rama_table = self.tables[new_to_old_conversion[proxy.residue_type]]
        residuals_array[i] = rama_table.compute_gradients(
          gradient_array=gradient_array,
          sites_cart=sites_cart,
          proxy=proxy,
          weight=self.params.rama_weight,
          epsilon=0.001)
    return flex.sum(residuals_array)

  def get_n_proxies(self):
    if self.proxies is not None:
      return self.proxies.size()
    else:
      return 0

  def _get_sorted_proxies_for_show(self,
      by_value,
      sites_cart,
      site_labels=None):
    class rama_for_show:
      def __init__(self,
          labels,
          residual):
        adopt_init_args(self, locals())
    assert by_value in ["residual", "delta"]
    assert site_labels is None or len(site_labels) == sites_cart.size()
    if self.get_n_proxies() == 0:
      return
    residuals_array = flex.double(self.proxies.size())
    self.target_and_gradients(
        unit_cell=None,
        sites_cart=sites_cart,
        residuals_array=residuals_array)
    result = []
    labels = site_labels if site_labels is not None \
        else [str(i) for i in range(sites_cart.size())]
    for i, pr in enumerate(self.proxies):
      i_seqs = pr.get_i_seqs()
      result.append(rama_for_show(
          [labels[i_seqs[0]],
           labels[i_seqs[1]],
           labels[i_seqs[2]],
           labels[i_seqs[3]],
           labels[i_seqs[4]]],
          residuals_array[i]))
    if by_value == "residual":
      result.sort(key=lambda x: x.residual, reverse=True)
    return result

  def show_sorted(self,
      by_value,
      sites_cart,
      site_labels=None,
      proxy_label=None, # not used yet
      f=None,
      prefix="",
      max_items=None):
    if self.get_n_proxies() == 0:
      return
    if (f is None): f = sys.stdout
    if by_value != "residual":
      by_value = "residual"
    sorted_proxies_for_show = self._get_sorted_proxies_for_show(
      by_value=by_value,
      sites_cart=sites_cart,
      site_labels=site_labels)
    print >> f, "Ramachandran plot restraints: %d"%len(sorted_proxies_for_show)
    print >> f, "Sorted by %s:" % by_value
    for p in sorted_proxies_for_show:
      print >> f, "phi-psi angles formed by             residual"
      print >> f, "    %s            %5.2e" % (p.labels[0], p.residual)
      for i in range(1,5):
        print >> f, "    %s" % p.labels[i]
    print >> f, ""


def load_tables (params=None) :
  import boost.python
  if (params is None) :
    params = master_phil.fetch().extract()
  if (params.scale_allowed <= 0.0) :
    raise Sorry("Ramachandran restraint parameter scale_allowed must be "+
      "a positive number (current value: %g)." % params.scale_allowed)
  tables = {}
  for residue_type in ["ala", "gly", "prepro", "pro"] :
    file_name = libtbx.env.find_in_repositories(
      relative_path="chem_data/rotarama_data/%s.rama.combined.data" %
        residue_type,
      test=os.path.isfile)
    f = open(file_name, "r")
    data = flex.double()
    for line in f.readlines() :
      val, phi, psi = line.split()
      assert ((int(phi) % 2 == 1) and (int(psi) % 2 == 1))
      data.append(float(val))
    t = lookup_table(data, 180, params.scale_allowed)
    tables[residue_type] = t
  return tables

class ramachandran_plot_data(object):
  def __init__(self):
    self.general = None
    self.gly = None
    self.cispro = None
    self.transpro = None
    self.prepro = None
    self.ileval = None
    stuff = [None]*6
    data = [("general", "rama8000-general-noGPIVpreP.data", 0),
        ("glycine", "rama8000-gly-sym.data", 1),
        ("cis-proline", "rama8000-cispro.data", 2),
        ("trans-proline", "rama8000-transpro.data", 3),
        ("pre-proline", "rama8000-prepro-noGP.data", 4),
        ("isoleucine or valine", "rama8000-ileval-nopreP.data", 5)]

    for (rama_key, file_name, selfstore) in data:
      file_name = libtbx.env.find_in_repositories(
        relative_path="chem_data/rotarama_data/%s" % file_name,
        test = os.path.isfile)
      stuff[selfstore] = flex.vec3_double()
      fo = open(file_name, "r")
      for line in fo.readlines():
        line = line.split()
        if(len(line)==3):
          phi_, psi_, val = float(line[0]),float(line[1]),float(line[2])
          triplet = [phi_, psi_, val]
          stuff[selfstore].append(triplet)

    # This strange cutting of tables at arbitrary thresholds is necessary
    # to make underlying 'algorithm' for gradient calculation
    # work in reasonable time (~2*n where n is number of points in the table,
    # for every gradient for each proxy).
    # Essentially, for each region on Ramachandran
    # plot it selects points with highest values (>0.5-0.7 of maximum region
    # value) and works only with them. This selection has nothing to do with
    # actual limits of allowed/favored regions. Funny thing is as soon as
    # phi-psi angles arrive to selected area, corresponding restraint is
    # effectively disabled by weight - calculations are still performed.
    self.general = self.normalize_general(data=stuff[0])
    self.gly     = self.normalize_gly(data=stuff[1])
    self.cispro  = self.normalize_pro(data=stuff[2])
    self.transpro= self.normalize_pro(data=stuff[3])
    self.prepro  = self.normalize_prepro(data=stuff[4])
    self.ileval  = self.normalize_general(data=stuff[5])

  def norm_to_max(self, data, val, sel, threshold=0.5):
    vmax = flex.max(val.select(sel))
    sel1 = val > vmax*threshold
    return data.select((sel1&sel))

  def thin_data(self, x, step = 1):
    result = x.select(flex.size_t(range(0,x.size(),step)))
    return result

  def split_array(self, data):
    phi = flex.double()
    psi = flex.double()
    val = flex.double()
    for x,y,z in data:
      phi.append(x)
      psi.append(y)
      val.append(z)
    return phi, psi, val

  def normalize_general(self, data):
    phi, psi, val = self.split_array(data=data)
    s0=(phi>0)&(phi< 180) & (psi<  -5)&(psi>-180)
    s1=(phi>0)&(phi< 180) & (psi>  -5)&(psi< 65)
    s2=(phi<0)&(phi>-180) & (psi<-100)&(psi>-180)
    s3=(phi<0)&(phi>-180) & (psi> -65)&(psi<  50)
    s4=(phi<0)&(phi>-180) & (psi>  50)&(psi< 180)
    s5=(phi<0)&(phi>-180) & (psi<-65)&(psi>-100)
    s6=(phi>0)&(phi< 180) & (psi>65)&(psi< 180)
    d0 = self.norm_to_max(data=data, val=val, sel=s0, threshold=0.7)
    d1 = self.norm_to_max(data=data, val=val, sel=s1, threshold=0.5)
    d2 = self.norm_to_max(data=data, val=val, sel=s2, threshold=0.5)
    d3 = self.norm_to_max(data=data, val=val, sel=s3, threshold=0.5)
    d4 = self.norm_to_max(data=data, val=val, sel=s4, threshold=0.5)
    d5 = self.norm_to_max(data=data, val=val, sel=s5, threshold=0.7)
    d6 = self.norm_to_max(data=data, val=val, sel=s6, threshold=0.7)
    d1.extend(d5)
    d1.extend(d6)
    d1.extend(d0)
    d1.extend(d2)
    d1.extend(d3)
    d1.extend(d4)
    return self.thin_data(d1)

  def normalize_prepro(self, data):
    phi, psi, val = self.split_array(data=data)
    s1 =(phi<0)&(phi>-180) & (psi<10)&(psi>-110)
    s2 =(phi<0)&(phi>-180) & (psi>10)&(psi<180)
    s3 =(phi>0)&(phi<90) & (psi> 0)&(psi< 110)
    s6 =(phi<0)&(phi>-180)& (psi<-150)&(psi>-180)
    s5 =(phi>120)&(phi<180)& (psi>-180)&(psi<180)
    s4 =(phi>0)&(phi<120)& (psi>120)&(psi<180)
    d1 = self.norm_to_max(data=data, val=val, sel=s1)
    d2 = self.norm_to_max(data=data, val=val, sel=s2)
    d3 = self.norm_to_max(data=data, val=val, sel=s3)
    d6 = self.norm_to_max(data=data, val=val, sel=s6)
    d4 = self.norm_to_max(data=data, val=val, sel=s4)
    d5 = self.norm_to_max(data=data, val=val, sel=s5)
    d1.extend(d4)
    d1.extend(d5)
    d1.extend(d2)
    d1.extend(d3)
    d1.extend(d6)
    return self.thin_data(d1)

  def normalize_pro(self, data):
    phi, psi, val = self.split_array(data=data)
    s1=(phi<0)&(phi>-130)& (psi<-100)&(psi>-180)
    s2=(phi<0)&(phi>-130)& (psi>-100)&(psi<  30)
    s3=(phi<0)&(phi>-130)& (psi>  30)&(psi<  90)
    s4=(phi<0)&(phi>-130)& (psi>  90)&(psi< 180)
    d1 = self.norm_to_max(data=data, val=val, sel=s1)
    d2 = self.norm_to_max(data=data, val=val, sel=s2)
    d3 = self.norm_to_max(data=data, val=val, sel=s3)
    d4 = self.norm_to_max(data=data, val=val, sel=s4)
    d1.extend(d2)
    d1.extend(d3)
    d1.extend(d4)
    return self.thin_data(d1)

  def normalize_gly(self, data):
    phi, psi, val = self.split_array(data=data)
    s1=(phi<0)&(phi>-180)& (psi<-90)&(psi>-180)
    s2=(phi>0)&(phi<180) & (psi<-90)&(psi>-180)
    s3=(phi<0)&(phi>-180)& (psi>-90)&(psi<70)
    s4=(phi>0)&(phi<180) & (psi>-90)&(psi<70)
    s5=(phi<0)&(phi>-180)& (psi>70)&(psi<180)
    s6=(phi>0)&(phi<180) & (psi>70)&(psi<180)
    d1 = self.norm_to_max(data=data, val=val, sel=s1)
    d2 = self.norm_to_max(data=data, val=val, sel=s2)
    d3 = self.norm_to_max(data=data, val=val, sel=s3)
    d4 = self.norm_to_max(data=data, val=val, sel=s4)
    d5 = self.norm_to_max(data=data, val=val, sel=s5)
    d6 = self.norm_to_max(data=data, val=val, sel=s6)
    d1.extend(d2)
    d1.extend(d3)
    d1.extend(d4)
    d1.extend(d5)
    d1.extend(d6)
    return self.thin_data(d1)


# def process_refinement_settings (
#     params,
#     pdb_hierarchy,
#     secondary_structure_manager,
#     d_min=None,
#     log=sys.stdout,
#     scope_name="refinement.pdb_interpretation.peptide_link") :
#   # seems that it is never used except in the test (which is not run...)
#   make_header("Extracting atoms for Ramachandran restraints", out=log)
#   from scitbx.array_family import flex
#   if (params.rama_selection is None) :
#     atom_selection = flex.bool(pdb_hierarchy.atoms().size(), True)
#   else :
#     cache = pdb_hierarchy.atom_selection_cache()
#     try :
#       sele = cache.selection(params.rama_selection)
#     except Exception, e :
#       raise Sorry("""Atom selection error:
#   %s
# Selection string resulting in error:
#   %s""" % (str(e), params.rama_selection))
#     else :
#       if (sele.count(True) == 0) :
#         raise Sorry("""Empty atom selection for %s.rama_selection.
# Current selection string:
#   %s""" % (scope_name, params.rama_selection))
#       atom_selection = sele
#   if params.rama_exclude_sec_str :
#     secondary_structure_manager.initialize(log=log)
#     alpha_sele = secondary_structure_manager.alpha_selection()
#     beta_sele = secondary_structure_manager.beta_selection()
#     ss_sele = alpha_sele | beta_sele
#     atom_selection &= ~ss_sele
#   return atom_selection
