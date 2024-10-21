from __future__ import absolute_import, division, print_function
import libtbx.load_env
import iotbx.phil
from libtbx import adopt_init_args
import sys
import os, math
import boost_adaptbx.boost.python as bp
from scitbx.array_family import flex
from mmtbx.validation import ramalyze
from mmtbx.conformation_dependent_library import generate_protein_threes
from six.moves import range
from collections import OrderedDict
from mmtbx.rotamer import ramachandran_eval

ext = bp.import_ext("mmtbx_ramachandran_restraints_ext")
from mmtbx_ramachandran_restraints_ext import lookup_table, \
    ramachandran_residual_sum, phi_psi_targets
ext2 = bp.import_ext("mmtbx_validation_ramachandran_ext")
from mmtbx_validation_ramachandran_ext import rama_eval

old_master_phil = iotbx.phil.parse("""
  rama_weight = 1.0
    .type = float
    .short_caption = Ramachandran gradients weight
    .expert_level = 1
    .style = hidden
  scale_allowed = 1.0
    .type = float
    .short_caption = Rescale allowed region pseudo-energy by
    .style = hidden
  rama_potential = *oldfield emsley
    .type = choice(multi=False)
    .short_caption = Ramachandran potential
    .caption = Oldfield Coot
    .style = hidden
  oldfield
    .short_caption = Oldfield potential parameters
    .style = box auto_align hidden
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
    .style = hidden
  restrain_rama_outliers = True
    .type = bool
    .help = Apply restraints to Ramachandran outliers
    .style = hidden
  restrain_rama_allowed = True
    .type = bool
    .help = Apply restraints to residues in allowed region on Ramachandran plot
    .style = hidden
  restrain_allowed_outliers_with_emsley = False
    .type = bool
    .help = In case of restrain_rama_outliers=True and/or restrain_rama_allowed=True \
      still restrain these residues with emsley. Make sense only in case of \
      using oldfield potential.
    .style = hidden
""")

master_phil = iotbx.phil.parse("""\
ramachandran_plot_restraints {
  enabled = False
    .type = bool
    .short_caption = Ramachandran restraints

  favored = *oldfield emsley emsley8k phi_psi_2
    .type = choice(multi=False)

  allowed = *oldfield emsley emsley8k phi_psi_2
    .type = choice(multi=False)

  outlier = *oldfield emsley emsley8k phi_psi_2
    .type = choice(multi=False)

  selection = None
    .type = atom_selection
    .short_caption = Atom selection for Ramachandran restraints
    .help = Selection of part of the model for which \
        Ramachandran restraints will be set up.
    .expert_level = 1
  inject_emsley8k_into_oldfield_favored = True
    .type=bool
    .expert_level = 3
    .help = Backdoor to disable temporary dirty hack to use both
  oldfield
    .short_caption = Oldfield settings
    .expert_level = 2
    .style = box
  {
    weight = 0.
      .type = float
      .expert_level = 2
      .help = Direct weight value. If 0 the weight will \
          be calculated as following: \
          (w, op.esd, op.dist_weight_max, 2.0, op.weight_scale) \
                                                                \
           1 / esd^2  *  max(2.0,   min(current_distance_to_allowed, dist_weight_max))    * weight_scale \
                         max(2.0,   current_distance_to_allowed) \
                                                                  \
           1 / esd^2  * weight_scale  *  max(distance_to_allowed_cutoff,   current_distance_to_allowed)   \
                weight_scale(=0.01)  *  max(distance_weight_min(=2.), min(distance_weight_max(=10.), current_distance_to_allowed))
    weight_scale = 0.01
      .type = float
      .expert_level = 2
    distance_weight_min = 2.0
      .type = float
      .expert_level = 2
      .help = minimum coefficient when scaling depending on how far the residue \
          is from allowed region.
    distance_weight_max = 10.0
      .type = float
      .expert_level = 2
      .help = maximum coefficient when scaling depending on how far the residue \
          is from allowed region.
    plot_cutoff = 0.027
      .type = float
      .expert_level = 2
  }
  emsley
    .short_caption = Emsley settings
    .style = box
  {
    weight = 1.0
      .type = float
      .short_caption = Ramachandran plot restraints weight (emsley)
      .expert_level = 1
    scale_allowed = 1.0
      .type = float
      .short_caption = Rescale allowed region pseudo-energy by
  }
  emsley8k
    .short_caption = Emsley8k settings
    .expert_level = 1
    .style = box
  {
    weight_favored = 5.0
      .type = float
      .short_caption = Ramachandran plot restraints weight (emsley8k, favored)
      .expert_level = 1
    weight_allowed = 10.0
      .type = float
      .short_caption = Ramachandran plot restraints weight (emsley8k, allowed)
      .expert_level = 1
    weight_outlier = 10.0
      .type = float
      .short_caption = Ramachandran plot restraints weight (emsley8k, outlier)
      .expert_level = 1
  }
  phi_psi_2
    .style = box
  {
    favored_strategy = *closest highest_probability random weighted_random
      .type = choice(multi=False)
    allowed_strategy = *closest highest_probability random weighted_random
      .type = choice(multi=False)
    outlier_strategy = *closest highest_probability random weighted_random
      .type = choice(multi=False)
  }
}
  """)

# Transformation from old to new parameters:
# weight = old.weight if (old.weight is None or old.weight > 0) else 0
# weight_scale = 1/old.esd^2 * old.weight_scale
# distance_to_allowed_cutoff = 2 if old.dist_weight_max > 2 else old.dist_weight_max
#

def is_proxy_present(proxies, n_seq, proxy):
  p_iseqs = list(proxy.get_i_seqs())
  ps = proxies.proxy_select(n_seq=n_seq,
      iselection=flex.size_t(p_iseqs))
  return ps.size() > 0

class ramachandran_manager(object):
  def __init__(self, pdb_hierarchy, params=None,
      log=sys.stdout,
      proxies=None, tables=None,
      initialize=True):
    assert pdb_hierarchy is not None
    assert not pdb_hierarchy.atoms().extract_i_seq().all_eq(0), ""+\
        "Probably all atoms have i_seq = 0 which is wrong"

    if params is None:
      # print ('init, params is None')
      w_params = master_phil.fetch().extract()
      w_params = w_params.ramachandran_plot_restraints
    elif hasattr(params, 'enabled'):
      # print ("init, hasattr(params, 'enabled')")
      # New params
      w_params = params
    elif hasattr(params, 'ramachandran_plot_restraints'):
      # print ("init, hasattr(params, 'ramachandran_plot_restraints'")
      # print ("init, ", type(params), type(params.ramachandran_plot_restraints), params.ramachandran_plot_restraints)
      w_params = params.ramachandran_plot_restraints
    else:
      # print ("init, else")
      w_params = master_phil.fetch().extract()
      w_params = w_params.ramachandran_plot_restraints
      # old params, make transfer
      w_params.selection = params.rama_selection
      # oldfield
      w_params.enabled = True
      w_params.oldfield.weight = \
          params.oldfield.weight if (params.oldfield.weight is None or params.oldfield.weight > 0) else 0
      w_params.oldfield.weight_scale = \
          1/(params.oldfield.esd**2) * params.oldfield.weight_scale
      w_params.oldfield.distance_weight_min = 2.0
      w_params.oldfield.distance_weight_max = params.oldfield.dist_weight_max

      # emsley
      w_params.emsley.weight = params.rama_weight
      w_params.emsley.scale_allowed = params.scale_allowed
      # strategy
      if params.rama_potential == 'oldfield':
        pass
      elif params.rama_potential == 'emsley':
        w_params.favored = 'emsley'
        w_params.allowed = 'emsley'
        w_params.outlier = 'emsley'
      if params.restrain_rama_outliers:
        w_params.outlier = params.rama_potential
      else:
        w_params.outlier = None
      if params.restrain_rama_allowed:
        w_params.allowed = params.rama_potential
      else:
        w_params.allowed = None
      if params.restrain_allowed_outliers_with_emsley:
        if not params.restrain_rama_allowed:
          w_params.allowed = 'emsley'
        if not params.restrain_rama_outliers:
          w_params.outlier = 'emsley'

    self.params = w_params
    self.rama_eval = rama_eval()
    self.hierarchy = pdb_hierarchy # only for def select()
    self.log = log
    self._oldfield_proxies = ext.shared_phi_psi_proxy()
    self._emsley_proxies   = ext.shared_phi_psi_proxy()
    self._emsley8k_proxies = ext.shared_phi_psi_proxy()
    self._phi_psi_2_proxies = ext.shared_phi_psi_proxy()
    self._oldfield_tables = None
    self._emsley_tables   = None
    self._emsley8k_tables = None
    self._phi_psi_2_tables = None
    if proxies is not None:
      self._oldfield_proxies, \
      self._emsley_proxies, \
      self._emsley8k_proxies, \
      self._phi_psi_2_proxies = proxies
    if tables is not None:
      self._oldfield_tables, \
      self._emsley_tables, \
      self._emsley8k_tables, \
      self._phi_psi_2_tables = tables
    self.initialize = initialize
    # bad hack to keep emsley potential in working(?) condition after
    # changing from rama500 to rama8000
    self.new_to_old_conversion = {"general":"ala", "glycine":"gly",
        "cis-proline":"pro", "trans-proline":"pro", "pre-proline":"prepro",
        "isoleucine or valine":"ala"}
    bool_atom_selection = self._determine_bool_atom_selection(pdb_hierarchy)
    fao = [self.params.favored, self.params.allowed, self.params.outlier]
    if initialize:
      if 'oldfield' in fao:
        self._oldfield_tables = ramachandran_plot_data(
          plot_cutoff=self.params.oldfield.plot_cutoff)
      if 'emsley' in fao:
        self._emsley_tables = load_tables()
      #
      ### THIS IS CRUEL. REMOVE ONCE favored/allowed/outlier are made multiple!
      #
      if 'emsley8k' in fao or self.params.inject_emsley8k_into_oldfield_favored:
        self._emsley8k_tables = load_emsley8k_tables()
      if 'phi_psi_2' in fao:
        self._phi_psi_2_tables = load_phi_psi_2_tables()
      # get proxies
      self.extract_proxies(pdb_hierarchy)
    if 'oldfield' in fao:
      self.target_phi_psi = self.update_phi_psi_targets_on_init(
        hierarchy = pdb_hierarchy)
    self.initialize = False

  def _determine_bool_atom_selection(self, hierarchy):
    bool_atom_selection = None
    if self.params.selection is None:
      bool_atom_selection = flex.bool(hierarchy.atoms_size(), True)
    else:
      cache = hierarchy.atom_selection_cache()
      bool_atom_selection = cache.selection(self.params.selection)
    return bool_atom_selection

  def proxy_select(self, n_seq, iselection):
    new_manager = ramachandran_manager(
      pdb_hierarchy = self.hierarchy,
      params        = self.params,
      log           = self.log,
      proxies = (None if self.get_n_oldfield_proxies() == 0 else self._oldfield_proxies.proxy_select(n_seq, iselection),
                 None if self.get_n_emsley_proxies()   == 0 else self._emsley_proxies.proxy_select(n_seq, iselection),
                 None if self.get_n_emsley8k_proxies() == 0 else self._emsley8k_proxies.proxy_select(n_seq, iselection),
                 None if self.get_n_phi_psi_2_proxies() == 0 else self._phi_psi_2_proxies.proxy_select(n_seq, iselection),
                 ),
      tables = (self._oldfield_tables, self._emsley_tables, self._emsley8k_tables, self._phi_psi_2_tables),
      initialize=False)
    return new_manager

  def extract_proxies(self, hierarchy):
    def _get_motifs():
      from phenix.programs.phi_psi_2 import results_manager as pp2
      pp2_manager = pp2(model=None, log=self.log)
      phi_psi_2_motifs = pp2_manager.get_overall_motif_count_and_output(
        None,
        self.hierarchy,
        return_rama_restraints=True,
        )
      return phi_psi_2_motifs
    phi_psi_2_motifs = None
    favored = ramalyze.RAMALYZE_FAVORED
    allowed = ramalyze.RAMALYZE_ALLOWED
    outlier = ramalyze.RAMALYZE_OUTLIER
    self.hierarchy = hierarchy
    bool_atom_selection = self._determine_bool_atom_selection(hierarchy)
    selected_h = hierarchy.select(bool_atom_selection)
    n_seq = flex.max(selected_h.atoms().extract_i_seq())
    # Drop all previous proxies
    self._oldfield_proxies = ext.shared_phi_psi_proxy()
    self._emsley_proxies   = ext.shared_phi_psi_proxy()
    self._emsley8k_proxies = ext.shared_phi_psi_proxy()
    self._phi_psi_2_proxies = ext.shared_phi_psi_proxy()
    # it would be great to save rama_eval, but the fact that this is called in
    # pdb_interpretation, not in mmtbx.model makes it impossible
    self.rama_eval = rama_eval()
    outl = []
    for three in generate_protein_threes(hierarchy=selected_h, geometry=None):
      rc = three.get_phi_psi_atoms()
      if rc is None: continue
      rama_key = three.get_ramalyze_key()
      angles = three.get_phi_psi_angles()
      rama_score = self.rama_eval.get_score(rama_key, angles[0], angles[1])
      r_eval = self.rama_eval.evaluate_score(rama_key, rama_score)
      phi_atoms, psi_atoms = rc
      i_seqs = [atom.i_seq for atom in phi_atoms] + [psi_atoms[-1].i_seq]
      resnames = three.get_resnames()
      r_name = resnames[1]
      assert rama_key in range(6)
      text_rama_key = ramalyze.res_types[rama_key]
      assert text_rama_key in ["general", "glycine", "cis-proline",
        "trans-proline", "pre-proline", "isoleucine or valine"]
      # pick where to put...
      ev_match_dict = {
        favored: self.params.favored,
        allowed: self.params.allowed,
        outlier: self.params.outlier}
      r_type = ev_match_dict[r_eval]
      if r_type == 'oldfield':
        proxy = ext.phi_psi_proxy(
          residue_type = text_rama_key,
          i_seqs       = i_seqs,
          weight       = 1) # XXX Not used in oldfield
        self.append_oldfield_proxies(proxy, n_seq)

        ### THIS IS CRUEL. REMOVE ONCE favored/allowed/outlier are made multiple!
        if(self.params.inject_emsley8k_into_oldfield_favored):
          proxy = ext.phi_psi_proxy(
            residue_type = text_rama_key,
            i_seqs       = i_seqs,
            weight       = 5)
          self.append_emsley8k_proxies(proxy, n_seq)
        ###

      elif r_type == 'emsley':
        weight = self.params.emsley.weight
        proxy = ext.phi_psi_proxy(
          residue_type = text_rama_key,
          i_seqs       = i_seqs,
          weight       = weight)
        self.append_emsley_proxies(proxy, n_seq)
      elif r_type == 'emsley8k':
        if(  r_eval is favored): weight=self.params.emsley8k.weight_favored
        elif(r_eval is allowed): weight=self.params.emsley8k.weight_allowed
        elif(r_eval is outlier): weight=self.params.emsley8k.weight_outlier
        else:                    raise RuntimeError("Rama eveluation failed.")
        proxy = ext.phi_psi_proxy(
          residue_type = text_rama_key,
          i_seqs       = i_seqs,
          weight       = weight)
        self.append_emsley8k_proxies(proxy, n_seq)
      elif r_type == 'phi_psi_2':
        from mmtbx.validation.phi_psi_2_data import get_phi_psi_key_for_rama_proxy
        if phi_psi_2_motifs is None: phi_psi_2_motifs = _get_motifs()
        if(  r_eval is favored): strategy=self.params.phi_psi_2.favored_strategy
        elif(r_eval is allowed): strategy=self.params.phi_psi_2.allowed_strategy
        elif(r_eval is outlier): strategy=self.params.phi_psi_2.outlier_strategy
        else:                    raise RuntimeError("Rama eveluation failed.")
        if strategy=='closest':
          strategy+='_%0.1f_%0.1f' % tuple(three.get_phi_psi_angles())
        pp2_key = get_phi_psi_key_for_rama_proxy(phi_psi_2_motifs,
                                                 three,
                                                 strategy=strategy,
                                                 )
        if pp2_key is None: continue
        weight=1
        proxy = ext.phi_psi_proxy(
          residue_type = pp2_key,
          i_seqs       = i_seqs,
          weight       = weight)
        outl.append([proxy.residue_type, three])
        self.append_phi_psi_2_proxies(proxy, n_seq)
      elif(r_type is None): pass
      else:
        raise RuntimeError("Not an option: %s"%str(r_type))

    print("", file=self.log)
    print("  %d Ramachandran restraints generated." % (
        self.get_n_proxies()), file=self.log)
    print("    %d Oldfield, %d Emsley, %d emsley8k and %d Phi/Psi/2." % (
      self.get_n_oldfield_proxies(),
      self.get_n_emsley_proxies(),
      self.get_n_emsley8k_proxies(),
      self.get_n_phi_psi_2_proxies()), file=self.log)
    if outl:
      print('    Rama restraints by Phi/Psi/2')
      for pp2, three in outl:
        print('      %s : %s' % (three[1].id_str(), pp2.split('|')[0]), file=self.log)

  @staticmethod
  def _append_proxies(proxies, proxy, n_seq):
    if not is_proxy_present(proxies, n_seq, proxy):
      proxies.append(proxy)

  def append_oldfield_proxies(self, proxy, n_seq):
    ramachandran_manager._append_proxies(self._oldfield_proxies, proxy, n_seq)

  def append_emsley_proxies(self, proxy, n_seq):
    ramachandran_manager._append_proxies(self._emsley_proxies, proxy, n_seq)

  def append_emsley8k_proxies(self, proxy, n_seq):
    ramachandran_manager._append_proxies(self._emsley8k_proxies, proxy, n_seq)

  def append_phi_psi_2_proxies(self, proxy, n_seq):
    ramachandran_manager._append_proxies(self._phi_psi_2_proxies, proxy, n_seq)

  def update_phi_psi_targets_on_init(self, hierarchy):
    if 'oldfield' in [self.params.favored, self.params.allowed, self.params.outlier]:
      sites_cart = hierarchy.atoms().extract_xyz()
      self.target_phi_psi = phi_psi_targets(
        sites_cart     = hierarchy.atoms().extract_xyz(),
        proxies        = self._oldfield_proxies,
        general_table  = self._oldfield_tables.general,
        gly_table      = self._oldfield_tables.gly,
        cispro_table   = self._oldfield_tables.cispro,
        transpro_table = self._oldfield_tables.transpro,
        prepro_table   = self._oldfield_tables.prepro,
        ileval_table   = self._oldfield_tables.ileval)
      return self.target_phi_psi
    return None

  def update_phi_psi_targets(self, hierarchy):
    self.hierarchy = hierarchy
    if not self.initialize:
      self.extract_proxies(hierarchy)
    self.update_phi_psi_targets_on_init(hierarchy)

  def target_and_gradients(self,
      unit_cell,
      sites_cart,
      gradient_array=None):
    if(gradient_array is None):
      gradient_array = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
    overall_residual_sum = 0
    # Oldfield
    self.residuals_array_oldfield = None #residuals_array_oldfield
    n_oldfield_proxies = self.get_n_oldfield_proxies()
    oldfield_residual_sum = 0
    if n_oldfield_proxies > 0:
      if self.residuals_array_oldfield is None:
        self.residuals_array_oldfield = flex.double(n_oldfield_proxies, 0.)
      op = self.params.oldfield
      w = op.weight
      if w is None:
        w = 0.
      oldfield_residual_sum = ramachandran_residual_sum(
          sites_cart=sites_cart,
          proxies=self._oldfield_proxies,
          gradient_array=gradient_array,
          phi_psi_targets = self.target_phi_psi,
          weights=(w, op.weight_scale, op.distance_weight_min, op.distance_weight_max),
          residuals_array=self.residuals_array_oldfield)
      overall_residual_sum += oldfield_residual_sum
    # Emsley
    self.residuals_array_emsley = None #residuals_array_emsley
    n_emsley_proxies = self.get_n_emsley_proxies()
    if n_emsley_proxies > 0:
      if self.residuals_array_emsley is None:
        self.residuals_array_emsley = flex.double(n_emsley_proxies, 0.)
      #assert (self.params.rama_weight >= 0.0)
      # Moving these cycles to C++ part would speed them up only up to 10%
      for i, proxy in enumerate(self._emsley_proxies):
        rama_table = self._emsley_tables[self.new_to_old_conversion[proxy.residue_type]]
        self.residuals_array_emsley[i] = rama_table.compute_gradients(
            gradient_array=gradient_array,
            sites_cart=sites_cart,
            proxy=proxy,
            epsilon=0.001)
      overall_residual_sum += flex.sum(self.residuals_array_emsley)
    # emsley8k
    self.residuals_array_emsley8k = None
    n_emsley8k_proxies = self.get_n_emsley8k_proxies()
    if n_emsley8k_proxies > 0:
      if self.residuals_array_emsley8k is None:
        self.residuals_array_emsley8k = flex.double(n_emsley8k_proxies, 0.)
      for i, proxy in enumerate(self._emsley8k_proxies):
        rama_table = self._emsley8k_tables[proxy.residue_type]
        self.residuals_array_emsley8k[i] = rama_table.compute_gradients(
          gradient_array = gradient_array,
          sites_cart     = sites_cart,
          proxy          = proxy,
          epsilon        = 1.0) # XXX
      overall_residual_sum += flex.sum(self.residuals_array_emsley8k)
    # phi/psi/2
    from mmtbx.validation.phi_psi_2_data import get_rama_table
    self.residuals_array_phi_psi_2 = None
    n_phi_psi_2_proxies = self.get_n_phi_psi_2_proxies()
    if n_phi_psi_2_proxies:
      if self.residuals_array_phi_psi_2 is None:
        self.residuals_array_phi_psi_2 = flex.double(n_phi_psi_2_proxies, 0.)
      for i, proxy in enumerate(self._phi_psi_2_proxies):
        # assert self.params.phi_psi_2.outlier_strategy=='highest_probability'
        rama_table = get_rama_table(
          proxy,
          self._phi_psi_2_tables,
          )
        if rama_table is None: continue
        self.residuals_array_phi_psi_2[i] = rama_table.compute_gradients(
          gradient_array = gradient_array,
          sites_cart     = sites_cart,
          proxy          = proxy,
          epsilon        = 1.0) # XXX
      overall_residual_sum += flex.sum(self.residuals_array_phi_psi_2)
    return overall_residual_sum

  def get_n_oldfield_proxies(self):
    if self._oldfield_proxies is not None:
      return self._oldfield_proxies.size()
    return 0

  def get_n_emsley_proxies(self):
    if self._emsley_proxies is not None:
      return self._emsley_proxies.size()
    return 0

  def get_n_emsley8k_proxies(self):
    if self._emsley8k_proxies is not None:
      return self._emsley8k_proxies.size()
    return 0

  def get_n_phi_psi_2_proxies(self):
    if self._phi_psi_2_proxies is not None:
      return self._phi_psi_2_proxies.size()
    return 0

  def get_n_proxies(self):
    return self.get_n_emsley_proxies() + \
           self.get_n_oldfield_proxies() + \
           self.get_n_emsley8k_proxies() + self.get_n_phi_psi_2_proxies()

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
    self.target_and_gradients(
        unit_cell=None,
        sites_cart=sites_cart)
    result_oldfield = []
    result_emsley   = []
    result_emsley8k  = []
    result_phi_psi_2 = []
    labels = site_labels if site_labels is not None \
        else [str(i) for i in range(sites_cart.size())]
    for proxies, residual_array, result in [
        (self._oldfield_proxies, self.residuals_array_oldfield, result_oldfield),
        (self._emsley_proxies, self.residuals_array_emsley, result_emsley),
        (self._emsley8k_proxies, self.residuals_array_emsley8k, result_emsley8k),
        (self._phi_psi_2_proxies, self.residuals_array_phi_psi_2, result_phi_psi_2),
        ]:
      if proxies is not None and proxies.size() > 0:
        for i, pr in enumerate(proxies):
          i_seqs = pr.get_i_seqs()
          result.append(rama_for_show(
              [labels[i_seqs[0]],
               labels[i_seqs[1]],
               labels[i_seqs[2]],
               labels[i_seqs[3]],
               labels[i_seqs[4]]],
              residual_array[i]))
        if by_value == "residual":
          result.sort(key=lambda x: x.residual, reverse=True)
    return result_oldfield, result_emsley, result_emsley8k, result_phi_psi_2

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
    sorted_oldfield_proxies_for_show, \
    sorted_emsley_proxies_for_show, \
    sorted_emsley8k_proxies_for_show, \
    sorted_phi_psi_2_proxies_for_show = \
      self._get_sorted_proxies_for_show(
        by_value=by_value,
        sites_cart=sites_cart,
        site_labels=site_labels)
    for proxies, label in [
        (sorted_oldfield_proxies_for_show, "Oldfield"),
        (sorted_emsley_proxies_for_show, "Emsley"),
        (sorted_emsley8k_proxies_for_show, "emsley8k"),
        (sorted_phi_psi_2_proxies_for_show, "phi/psi/2"),
        ]:
      print("Ramachandran plot restraints (%s): %d" % (label, len(proxies)), file=f)
      print("Sorted by %s:" % by_value, file=f)
      for p in proxies:
        print("phi-psi angles formed by             residual", file=f)
        print("    %s            %5.2e" % (p.labels[0], p.residual), file=f)
        for i in range(1,5):
          print("    %s" % p.labels[i], file=f)
      print("", file=f)

def load_tables():
  tables = {}
  for residue_type in ["ala", "gly", "prepro", "pro"] :
    file_name = libtbx.env.find_in_repositories(
      relative_path="chem_data/rotarama_data/%s.rama.combined.data" %
        residue_type,
      test=os.path.isfile)
    data = flex.double()
    with open(file_name, "r") as f:
      lines = f.readlines()
    for line in lines:
      val, phi, psi = line.split()
      assert ((int(phi) % 2 == 1) and (int(psi) % 2 == 1))
      data.append(float(val))
    t = lookup_table(data, 180)
    tables[residue_type] = t
  return tables

def load_emsley8k_tables():
  tables = {}
  name_to_file = [
    ("general",              "rama8000-general-noGPIVpreP.data", 0),
    ("glycine",              "rama8000-gly-sym.data", 1),
    ("cis-proline",          "rama8000-cispro.data", 2),
    ("trans-proline",        "rama8000-transpro.data", 3),
    ("pre-proline",          "rama8000-prepro-noGP.data", 4),
    ("isoleucine or valine", "rama8000-ileval-nopreP.data", 5)]
  tmp = OrderedDict()
  rr = [i for i in range(-179,180,2)]
  for i in rr:
    for j in rr:
      tmp[(i,j)]=0
  R = ramachandran_eval.RamachandranEval()
  outlier = ramalyze.RAMALYZE_OUTLIER
  favored = ramalyze.RAMALYZE_FAVORED
  allowed = ramalyze.RAMALYZE_ALLOWED
  for (rama_key, file_name, selfstore) in name_to_file:
    file_name = libtbx.env.find_in_repositories(
      relative_path="chem_data/rotarama_data/%s"%(file_name),
      test=os.path.isfile)
    di        = {}
    outlier_vals = flex.double()
    favored_vals = flex.double()
    allowed_vals = flex.double()
    status       = {}
    with open(file_name, "r") as f:
      lines = f.readlines()
    for line in lines:
      if line[0]=="#": continue
      phi, psi, val = line.split()
      phi=int(float(phi))
      psi=int(float(psi))
      val=float(val)
      di[(phi,psi)]=val
      rama_score = R.rama_eval.get_score(selfstore, phi, psi)
      evaluation = R.rama_eval.evaluate_score(selfstore, rama_score)
      if  (evaluation==outlier):
        outlier_vals.append(val)
        status[(phi,psi)]=outlier
      elif(evaluation==favored):
        favored_vals.append(val)
        status[(phi,psi)]=favored
      elif(evaluation==allowed):
        allowed_vals.append(val)
        status[(phi,psi)]=allowed
      else: raise RuntimeError("Not supposed to be here.")
    data = flex.double()
    max_outlier = flex.max(outlier_vals)
    max_favored = flex.max(favored_vals)
    min_favored = flex.min(favored_vals)
    max_allowed = flex.max(allowed_vals)
    for k, v in zip(tmp.keys(), tmp.values()):
      try:
        val = di[k]
        if  (status[k]==outlier): val = -1 + val/max_outlier
        elif(status[k]==favored): val = val #math.exp(val)**0.5/2.71828182846**0.5 #math.exp(val)**3
        elif(status[k]==allowed): val = val #math.exp(val)**0.5/2.71828182846**0.5 #math.exp(val/max_allowed)**3
      except KeyError:
        val = -1
      data.append(val)
    t = lookup_table(data, 180)
    tables[rama_key] = t
  return tables

def load_phi_psi_2_tables():
  from mmtbx.validation.phi_psi_2_data import load_phi_psi_2_rama_restraints_tables
  return load_phi_psi_2_rama_restraints_tables()

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
      with open(file_name, "r") as f:
        lines = f.readlines()
      for line in lines:
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
