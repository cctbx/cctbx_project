
from __future__ import division
import libtbx.load_env
import iotbx.phil
from libtbx.str_utils import make_header
from libtbx.utils import Sorry
from libtbx import adopt_init_args
import sys
import os
import mmtbx.rotamer
import boost.python
from scitbx.array_family import flex


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

class ramachandran_manager(object):
  def __init__ (self, pdb_hierarchy, atom_selection, params=None, log=sys.stdout,
      proxies=None, tables=None, initialize=True):
    assert pdb_hierarchy is not None
    adopt_init_args(self, locals())
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

  def select(self, n_seq, iselection):
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
    assert [a.i_seq for a in self.pdb_hierarchy.atoms()].count(0) == 1 ,\
        "Probably all atoms have i_seq = 0 which is wrong"
    ext = boost.python.import_ext("mmtbx_ramachandran_restraints_ext")
    angles = mmtbx.rotamer.extract_phi_psi(
      pdb_hierarchy=self.pdb_hierarchy,
      atom_selection=self.atom_selection)
    self.proxies = ext.shared_phi_psi_proxy()
    for angle in angles :
      residue_name = angle.residue_name
      if (residue_name == "MSE") :
        residue_name = "MET"
      proxy = ext.phi_psi_proxy(
        residue_name=residue_name,
        residue_type=angle.residue_type,
        i_seqs=angle.i_seqs)
      self.proxies.append(proxy)
    print >> self.log, ""
    print >> self.log, "  %d Ramachandran restraints generated." % self.get_n_proxies()

  def restraints_residual_sum (self,
                               sites_cart,
                               gradient_array=None,
                               residuals_array=None) :
    from scitbx.array_family import flex
    assert self.proxies is not None
    ext = boost.python.import_ext("mmtbx_ramachandran_restraints_ext")
    from mmtbx_ramachandran_restraints_ext import rama_target_and_gradients, \
      target_phi_psi
    if(gradient_array is None) :
      gradient_array = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
    if residuals_array is None:
      residuals_array = flex.double()
    residuals_array.clear()
    target = 0
    if(self.params.rama_potential == "oldfield"):
      op = self.params.oldfield
      for proxy in self.proxies:
        if(proxy.residue_type=="gly"):    rama_table = self.tables.gly
        if(proxy.residue_type=="pro"):    rama_table = self.tables.pro
        if(proxy.residue_type=="prepro"): rama_table = self.tables.prepro
        if(proxy.residue_type=="ala"):    rama_table = self.tables.general
        r = target_phi_psi(
          rama_table     = rama_table,
          sites_cart     = sites_cart,
          proxy          = proxy)
        if(op.weight is None):
          weight = 1./(op.esd**2)*min(r[2],op.dist_weight_max)*op.weight_scale
        else: weight = op.weight
        tg = rama_target_and_gradients(
           gradient_array = gradient_array,
           phi_target     = r[0],
           psi_target     = r[1],
           weight         = weight,
           rama_table     = rama_table,
           sites_cart     = sites_cart,
           proxy          = proxy)
        residuals_array.append(tg.target())
    elif (self.params.rama_potential == "rosetta") :
      for proxy in self.proxies :
        residuals_array.append(self.tables.residue_target_and_gradients(
          gradient_array=gradient_array,
          sites_cart=sites_cart,
          proxy=proxy,
          weight=self.params.rama_weight))
    else:
      assert (self.params.rama_weight >= 0.0)
      for proxy in self.proxies :
        rama_table = self.tables[proxy.residue_type]
        residuals_array.append(rama_table.compute_gradients(
          gradient_array=gradient_array,
          sites_cart=sites_cart,
          proxy=proxy,
          weight=self.params.rama_weight,
          epsilon=0.001))
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
    residuals_array = flex.double()
    self.restraints_residual_sum(
        sites_cart,
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
    print >> f, "Ramachandran plot restraints: %d" % len(sorted_proxies_for_show)
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
  ext = boost.python.import_ext("mmtbx_ramachandran_restraints_ext")
  from mmtbx_ramachandran_restraints_ext import lookup_table
  from scitbx.array_family import flex
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
    from scitbx.array_family import flex
    self.gly = None
    self.pro = None
    self.prepro = None
    self.general = None
    files = [
      ["gly",     "rama500-gly-sym.data"],
      ["pro",     "rama500-pro.data"],
      ["prepro",  "rama500-prepro.data"],
      ["general", "rama500-general.data"]]
    for rtf in files:
      rt,f = rtf
      file_name = libtbx.env.find_in_repositories(
        relative_path="chem_data/rotarama_data/%s"%f,
        test = os.path.isfile)
      if(rt=="gly"):     self.gly     = flex.vec3_double()
      if(rt=="pro"):     self.pro     = flex.vec3_double()
      if(rt=="prepro"):  self.prepro  = flex.vec3_double()
      if(rt=="general"): self.general = flex.vec3_double()
      fo = open(file_name, "r")
      for line in fo.readlines():
        line = line.split()
        if(len(line)==3):
          phi_, psi_, val = float(line[0]),float(line[1]),float(line[2])
          triplet = [phi_, psi_, val]
          if(rt=="gly"):     self.gly    .append(triplet)
          if(rt=="pro"):     self.pro    .append(triplet)
          if(rt=="prepro"):  self.prepro .append(triplet)
          if(rt=="general"): self.general.append(triplet)
    #
    self.gly = self.normalize_gly(data=self.gly)
    self.pro = self.normalize_pro(data=self.pro)
    self.prepro = self.normalize_prepro(data=self.prepro)
    self.general = self.normalize_general(data=self.general)

  def norm_to_max(self, data, val, sel, threshold=0.5):
    from scitbx.array_family import flex
    vmax = flex.max(val.select(sel))
    sel1 = val > vmax*threshold
    return data.select((sel1&sel))

  def thin_data(self, x, step = 1):
    from scitbx.array_family import flex
    return x.select(flex.size_t(range(0,x.size(),step)))

  def split_array(self, data):
    from scitbx.array_family import flex
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
