
from __future__ import division
import libtbx.load_env
import libtbx.phil
from libtbx.utils import Sorry
from libtbx import adopt_init_args
import sys
import os
from scitbx.array_family import flex

import boost.python
ext = boost.python.import_ext("mmtbx_ramachandran_ext")
from mmtbx_ramachandran_ext import *

ext = boost.python.import_ext("mmtbx_ramachandran_restraints_ext")
from mmtbx_ramachandran_restraints_ext import target_and_gradients, target_phi_psi

master_phil = libtbx.phil.parse("""
  rama_weight = 1.0
    .type = float
    .short_caption = Ramachandran gradients weight
    .expert_level = 1
  scale_allowed = 1.0
    .type = float
    .short_caption = Rescale allowed region pseudo-energy by
  use_finite_differences = False
    .type = bool
    .help = Used for testing - not suitable for real structures.
    .short_caption = Use finite differences (DEVELOPERS ONLY)
    .expert_level = 3
#  type = *oldfield emsley rosetta
   type = *oldfield emsley
    .type = choice(multi=False)
  oldfield
  {
    esd = 10.0
      .type = float
      .expert_level = 2
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
#  rosetta_scale = 5.0
#    .type = float
""")

refine_opt_params = libtbx.phil.parse("""
#  min_allowed_d_min = 3.0
#    .type = float
#    .short_caption = Resolution cutoff for Ramachandran restraints
#    .expert_level = 2
  rama_selection = None
    .type = str
    .short_caption = Atom selection for Ramachandran restraints
    .style = selection
    .expert_level = 1
  exclude_secondary_structure = False
    .type = str
    .expert_level = 1
""")

def load_tables (params=None) :
  if (params is None) :
    params = master_phil.fetch().extract()
  if (params.scale_allowed <= 0.0) :
    raise Sorry("Ramachandran restraint parameter scale_allowed must be "+
      "a positive number (current value: %g)." % params.scale_allowed)
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

class generic_proxy (object) :
  restraint_type = None

class proxy (generic_proxy) :
  restraint_type = "ramachandran"
  def __init__ (self, i_seqs, residue_type, residue_name) :
    assert (len(i_seqs) == 5)
    self.i_seqs = i_seqs
    self.residue_type = residue_type
    # XXX the Rosetta AA lookup function will crash if passed a non-standard
    # residue name, so it is validated first if Rosetta is being used
    self.residue_name = residue_name
    self.ignore_flag = None

def show_histogram(data, n_slots, log):
  hm = flex.histogram(data = data, n_slots = n_slots)
  lc_1 = hm.data_min()
  s_1 = enumerate(hm.slots())
  print >> log, "      Map Values          Number of points"
  for (i_1,n_1) in s_1:
    hc_1 = hm.data_min() + hm.slot_width() * (i_1+1)
    print >> log, "%10.3f - %-10.3f : %d" % (lc_1, hc_1, n_1)
    lc_1 = hc_1

class ramachandran_plot_data(object):
  def __init__(self):
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
    vmax = flex.max(val.select(sel))
    sel1 = val > vmax*threshold
    return data.select((sel1&sel))

  def thin_data(self, x, step = 1):
    return x.select(flex.size_t(range(0,x.size(),step)))

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

class generic_restraints_helper (object) :
  def __init__ (self, params) :
    adopt_init_args(self, locals())
    if(self.params.type == "oldfield"):
      self.tables = ramachandran_plot_data()
    elif (self.params.type == "rosetta") :
      try :
        from rosetta_adaptbx import scoring
      except ImportError, e :
        print e
        raise Sorry("Rosetta not included or not configured in this build.")
      else :
        self.tables = scoring.ramachandran() # not a list
    else :
      self.tables = load_tables(params)

  def restraints_residual_sum (self,
                               sites_cart,
                               proxies,
                               gradient_array=None,
                               unit_cell=None) :
    from scitbx.array_family import flex
    ramachandran_proxies = []
    for proxy in proxies :
      if (proxy.restraint_type == "ramachandran") :
        ramachandran_proxies.append(proxy)
    if(self.params.type == "oldfield"):
      op = self.params.oldfield
      if(gradient_array is None) :
        from scitbx.array_family import flex
        gradient_array = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
      target = 0
      for proxy in ramachandran_proxies:
        if(proxy.residue_type=="gly"):    rama_table = self.tables.gly
        if(proxy.residue_type=="pro"):    rama_table = self.tables.pro
        if(proxy.residue_type=="prepro"): rama_table = self.tables.prepro
        if(proxy.residue_type=="ala"):    rama_table = self.tables.general
        r = target_phi_psi(
          rama_table     = rama_table,
          sites_cart     = sites_cart,
          i_seqs         = proxy.i_seqs)
        if(op.weight is None):
          weight = 1./(op.esd**2)*min(r[2],op.dist_weight_max)*op.weight_scale
        else: weight = op.weight
        tg = target_and_gradients(
           gradient_array = gradient_array,
           phi_target     = r[0],
           psi_target     = r[1],
           weight         = weight,
           rama_table     = rama_table,
           sites_cart     = sites_cart,
           i_seqs         = proxy.i_seqs)
        target += tg.target()
      return target
    elif (self.params.type == "rosetta") :
      return self.target_and_gradients_rosetta(
        sites_cart=sites_cart,
        proxies=ramachandran_proxies,
        gradient_array=gradient_array)
    else:
      return self._phi_psi_restraints_residual_sum(
        sites_cart=sites_cart,
        proxies=ramachandran_proxies,
        gradient_array=gradient_array)

  def _phi_psi_restraints_residual_sum (self,
                                        proxies,
                                        sites_cart,
                                        gradient_array=None) :
    if (gradient_array is None) :
      from scitbx.array_family import flex
      gradient_array = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
    sum = 0
    assert (self.params.rama_weight >= 0.0)
    for proxy in proxies :
      rama_table = self.tables[proxy.residue_type]
      if self.params.use_finite_differences :
        sum += rama_table.compute_gradients_finite_differences(
          gradient_array=gradient_array,
          sites_cart=sites_cart,
          i_seqs=proxy.i_seqs,
          weight=self.params.rama_weight,
          epsilon=0.001)
      else :
        sum += rama_table.compute_gradients(
          gradient_array=gradient_array,
          sites_cart=sites_cart,
          i_seqs=proxy.i_seqs,
          weight=self.params.rama_weight,
          epsilon=0.001)
    return sum

  def target_and_gradients_rosetta (self,
                                    proxies,
                                    sites_cart,
                                    gradient_array=None) :
    from iotbx import pdb
    if (gradient_array is None) :
      from scitbx.array_family import flex
      gradient_array = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
    sum = 0
    assert (self.params.rama_weight >= 0.0)
    weight_scale = self.params.rosetta_scale
    if (weight_scale is None) :
      weight_scale = 2.0
    for proxy in proxies :
      if proxy.ignore_flag : # non-standard AA
        continue
      elif (proxy.ignore_flag is None) :
        if (not proxy.residue_name in pdb.common_residue_names_amino_acid) :
          proxy.ignore_flag = True
          continue
        else :
          proxy.ignore_flag = False
      sum += self.tables.target_and_gradients(
        gradient_array=gradient_array,
        sites_cart=sites_cart,
        residue_name=proxy.residue_name,
        i_seqs=proxy.i_seqs,
        weight=self.params.rama_weight * weight_scale)
    return sum

def extract_proxies (pdb_hierarchy,
                     atom_selection=None,
                     log=sys.stdout) :
  from iotbx.pdb.amino_acid_codes import one_letter_given_three_letter
  from cctbx import geometry_restraints
  from scitbx.array_family import flex
  if (atom_selection is None) :
    atom_selection = flex.bool(pdb_hierarchy.atoms().size(), True)
  proxies = []
  for model in pdb_hierarchy.models() :
    for chain in model.chains() :
      for conformer in chain.conformers() :
        residues = conformer.residues()
        for i, residue in enumerate(residues) :
          if (not residue.resname in one_letter_given_three_letter) :
            continue
          next_res, prev_res = None, None
          resseq2 = residue.resseq_as_int()
          resseq1, resseq3 = None, None
          if (i > 0):
            resseq1 = residues[i-1].resseq_as_int()
            if (resseq2 != (resseq1 + 1)) :
              continue
            prev_res = residues[i-1]
          if (i < (len(residues) - 1)) :
            resseq3 = residues[i+1].resseq_as_int()
            if (resseq2 != (resseq3 - 1)) :
              continue
            next_res = residues[i+1]
          if (next_res is not None) and (prev_res is not None) :
            c1, n2, ca2, c2, n3 = None, None, None, None, None
            for atom in prev_res.atoms() :
              if (atom.name == " C  ") :
                c1 = atom
                break
            for atom in residue.atoms() :
              if (atom.name == " N  ") :
                n2 = atom
              elif (atom.name == " CA ") :
                ca2 = atom
              elif (atom.name == " C  ") :
                c2 = atom
            for atom in next_res.atoms() :
              if (atom.name == " N  ") :
                n3 = atom
            if (None in [c1, n2, ca2, c2, n3]) :
              #print >> log, "  incomplete backbone for %s %d-%d, skipping." % \
              #  (chain.id, resseq1, resseq3)
              continue
            i_seqs = [c1.i_seq,n2.i_seq,ca2.i_seq,c2.i_seq,n3.i_seq]
            for i_seq in i_seqs :
              if (not atom_selection[i_seq]) :
                continue
            pep1 = geometry_restraints.bond(
              sites=[c1.xyz,n2.xyz],
              distance_ideal=1,
              weight=1)
            pep2 = geometry_restraints.bond(
              sites=[c2.xyz,n3.xyz],
              distance_ideal=1,
              weight=1)
            if (pep1.distance_model > 4) or (pep2.distance_model > 4) :
              continue
            residue_name = residue.resname
            if (residue_name == "PRO") :
              residue_type = "pro"
            elif (residue_name == "GLY") :
              residue_type = "gly"
            elif (residues[i+1].resname == "PRO") :
              residue_type = "prepro"
            else :
              residue_type = "ala"
            if (residue_name == "MSE") :
              residue_name = "MET"
            phi_psi = proxy(i_seqs, residue_type, residue_name)
            proxies.append(phi_psi)
  print >> log, "%d Ramachandran restraints generated." % len(proxies)
  return proxies

def process_refinement_settings (
    params,
    pdb_hierarchy,
    secondary_structure_manager,
    d_min=None,
    log=sys.stdout,
    scope_name="refinement.ramachandran_restraints") :
  if (params.rama_selection is None) :
    atom_selection = flex.bool(pdb_hierarchy.atoms().size(), True)
  else :
    cache = pdb_hierarchy.atom_selection_cache()
    try :
      sele = cache.selection(params.rama_selection)
    except Exception, e :
      raise Sorry("""Atom selection error:
  %s
Selection string resulting in error:
  %s""" % (str(e), params.rama_selection))
    else :
      if (sele.count(True) == 0) :
        raise Sorry("""Empty atom selection for %s.rama_selection.
Current selection string:
  %s""" % (scope_name, params.rama_selection))
      atom_selection = sele
  if params.exclude_secondary_structure :
    alpha_sele = secondary_structure_manager.alpha_selection()
    beta_sele = secondary_structure_manager.beta_selection()
    ss_sele = alpha_sele | beta_sele
    atom_selection &= ~ss_sele
  return atom_selection
