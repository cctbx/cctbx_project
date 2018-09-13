from __future__ import division
import libtbx.load_env
import iotbx.phil
from libtbx.utils import Sorry
from libtbx import adopt_init_args
import sys
import os
import boost.python
from scitbx.array_family import flex
from mmtbx.validation import ramalyze
from mmtbx.conformation_dependent_library import generate_protein_threes

ext = boost.python.import_ext("mmtbx_ramachandran_restraints_ext")
from mmtbx_ramachandran_restraints_ext import lookup_table, \
    ramachandran_residual_sum, phi_psi_targets
ext2 = boost.python.import_ext("mmtbx_validation_ramachandran_ext")
from mmtbx_validation_ramachandran_ext import rama_eval

master_phil = iotbx.phil.parse("""
  rama_weight = 1.0
    .type = float
    .short_caption = Ramachandran gradients weight
    .expert_level = 1
  scale_allowed = 1.0
    .type = float
    .short_caption = Rescale allowed region pseudo-energy by
  rama_potential = *oldfield emsley
    .type = choice(multi=False)
    .short_caption = Ramachandran potential
    .caption = Oldfield Coot
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
    plot_cutoff = 0.027
      .type = float
      .expert_level = 2
  }
  rama_selection = None
    .type = atom_selection
    .short_caption = Atom selection for Ramachandran restraints
    .help = Selection of part of the model for which \
        Ramachandran restraints will be set up.
    .expert_level = 1
  restrain_rama_outliers = True
    .type = bool
    .help = Apply restraints to Ramachandran outliers
  restrain_rama_allowed = True
    .type = bool
    .help = Apply restraints to residues in allowed region on Ramachandran plot
""")

def is_proxy_present(proxies, n_seq, proxy):
  p_iseqs = list(proxy.get_i_seqs())
  ps = proxies.proxy_select(n_seq=n_seq,
      iselection=flex.size_t(p_iseqs))
  return ps.size() > 0

class ramachandran_manager(object):
  def __init__ (self, pdb_hierarchy, params=None,
      log=sys.stdout, proxies=None, tables=None, initialize=True):
    assert pdb_hierarchy is not None
    assert not pdb_hierarchy.atoms().extract_i_seq().all_eq(0), ""+\
        "Probably all atoms have i_seq = 0 which is wrong"
    adopt_init_args(self, locals(), exclude=["log"])
    if self.params is None:
      self.params = master_phil.fetch().extract()
    self.bool_atom_selection = None
    if self.params.rama_selection is None:
      self.bool_atom_selection = flex.bool(pdb_hierarchy.atoms_size(), True)
    else:
      cache = pdb_hierarchy.atom_selection_cache()
      self.bool_atom_selection = cache.selection(self.params.rama_selection)
    if initialize:
      if(self.params.rama_potential == "oldfield"):
        self.tables = ramachandran_plot_data(
            plot_cutoff=self.params.oldfield.plot_cutoff)
      else :
        self.tables = load_tables(params)
      # get proxies
      self.extract_proxies(log=log)
    else:
      assert proxies is not None
    if(self.params.rama_potential == "oldfield"):
      self.target_phi_psi = self.update_phi_psi_targets(
        sites_cart = self.pdb_hierarchy.atoms().extract_xyz())

  def proxy_select(self, n_seq, iselection, log=sys.stdout):
    result_proxies = self.proxies.proxy_select(n_seq, iselection)
    return ramachandran_manager(
        pdb_hierarchy=self.pdb_hierarchy,
        params=self.params,
        log=log,
        proxies=result_proxies,
        tables=self.tables,
        initialize=False)

  def extract_proxies(self, log):
    self.proxies = ext.shared_phi_psi_proxy()
    selected_h = self.pdb_hierarchy.select(self.bool_atom_selection)
    n_seq = flex.max(selected_h.atoms().extract_i_seq())
    # it would be great to save rama_eval, but the fact that this is called in
    # pdb_interpretation, not in mmtbx.model makes it impossible
    need_filtering = not (self.params.restrain_rama_outliers and
                          self.params.restrain_rama_allowed)
    if need_filtering:
      self.rama_eval = rama_eval()
    for three in generate_protein_threes(
        hierarchy=selected_h,
        geometry=None):
      rc = three.get_phi_psi_atoms()
      if rc is None: continue
      rama_key = three.get_ramalyze_key()
      if need_filtering:
        angles = three.get_phi_psi_angles()
        rama_score = self.rama_eval.get_score(rama_key, angles[0], angles[1])
        r_evaluation = self.rama_eval.evaluate_score(rama_key, rama_score)
        if (r_evaluation == ramalyze.RAMALYZE_OUTLIER and
            not self.params.restrain_rama_outliers):
          continue
        if (r_evaluation == ramalyze.RAMALYZE_ALLOWED and
            not self.params.restrain_rama_allowed):
          continue
      phi_atoms, psi_atoms = rc
      i_seqs = [atom.i_seq for atom in phi_atoms] + [psi_atoms[-1].i_seq]
      resnames = three.get_resnames()
      r_name = resnames[1]
      assert rama_key in range(6)
      text_rama_key = ramalyze.res_types[rama_key]
      assert text_rama_key in ["general", "glycine", "cis-proline", "trans-proline",
                 "pre-proline", "isoleucine or valine"]
      proxy = ext.phi_psi_proxy(
          residue_name=r_name,
          residue_type=text_rama_key,
          i_seqs=i_seqs)
      if not is_proxy_present(self.proxies, n_seq, proxy):
        self.proxies.append(proxy)
    print >> log, ""
    print >> log, "  %d Ramachandran restraints generated." % (
        self.get_n_proxies())

  def update_phi_psi_targets(self, sites_cart):
    if(self.params.rama_potential != "oldfield"): return None
    self.target_phi_psi = phi_psi_targets(
      sites_cart=sites_cart,
      proxies=self.proxies,
      general_table=self.tables.general,
      gly_table=self.tables.gly,
      cispro_table=self.tables.cispro,
      transpro_table=self.tables.transpro,
      prepro_table=self.tables.prepro,
      ileval_table=self.tables.ileval)
    return self.target_phi_psi

  def target_and_gradients(self,
      unit_cell,
      sites_cart,
      gradient_array=None,
      residuals_array=None):
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
        phi_psi_targets = self.target_phi_psi,
        weights=(w, op.esd, op.dist_weight_max, 2.0, op.weight_scale),
        residuals_array=residuals_array)
      return res
    else: # emsley
      #assert (self.params.rama_weight >= 0.0)
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
  def __init__(self, plot_cutoff=0.027):
    self.plot_cutoff = plot_cutoff
    self.general = None
    self.gly = None
    self.cispro = None
    self.transpro = None
    self.prepro = None
    self.ileval = None
    stuff = [None]*6
    data = [
      ("general",              "rama8000-general-noGPIVpreP.data", 0),
      ("glycine",              "rama8000-gly-sym.data", 1),
      ("cis-proline",          "rama8000-cispro.data", 2),
      ("trans-proline",        "rama8000-transpro.data", 3),
      ("pre-proline",          "rama8000-prepro-noGP.data", 4),
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
    self.general  = self.select_good(data=stuff[0], step=1)
    self.gly      = self.select_good(data=stuff[1], step=1)
    self.cispro   = self.select_good(data=stuff[2], step=1)
    self.transpro = self.select_good(data=stuff[3], step=1)
    self.prepro   = self.select_good(data=stuff[4], step=1)
    self.ileval   = self.select_good(data=stuff[5], step=1)

  def select_good(self, data, step):
    phi, psi, val = self.split_array(data=data)
    # 0.02 is border for favored, 0.007 - arbitrary addition to ensure more
    # residues in favored region. Moved to plot_cutoff.
    sel = (val>self.plot_cutoff)
    phi_values = set(list(phi))
    psi_values = set(list(psi))
    phi_list = sorted(phi_values)
    psi_list = sorted(psi_values)
    needed_phi_values = set([phi_list[i] for i in range(0, len(phi_list),step)])
    needed_psi_values = set([psi_list[i] for i in range(0, len(psi_list),step)])
    result = flex.vec3_double()
    selected = data.select(sel)
    for phi, psi, val in selected:
      if phi in needed_phi_values and psi in needed_psi_values:
        result.append([phi, psi, val])
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
