
from __future__ import absolute_import, division, print_function
from mmtbx.building import alternate_conformations as alt_confs
import mmtbx.building
from libtbx.str_utils import make_sub_header
from libtbx import Auto, adopt_init_args
from libtbx.utils import null_out
from libtbx import easy_mp
import random
import time
import os
import sys

master_params_str = """
n_trials = 4
  .type = int
filter_rotamer_outliers = True
  .type = bool
cc_min = 0.7
  .type = float
min_fofc = 2.0
  .type = float
rsr_after_anneal = False
  .type = bool
adjust_b_factors_if_poor_density = False
  .type = bool
negate_surrounding = False
  .type = bool
partial_occupancy = 0.5
  .type = floats
shake_sites = None
  .type = float
simulated_annealing {
  include scope mmtbx.dynamics.simulated_annealing.master_params_str
}
"""

class refine_into_difference_density(object):
  """
  Driver for the annealing into mFo-DFc density with selected atoms at partial
  occupancy.  Can be run in parallel with different random seeds or occupancy.
  """
  def __init__(self,
      fmodel,
      pdb_hierarchy,
      processed_pdb_file,
      params,
      selection,
      selection_score=None,
      nproc=Auto,
      out=None):
    adopt_init_args(self, locals())
    from scitbx.array_family import flex
    if (self.out is None) : self.out = sys.stdout
    self.sites_start = fmodel.xray_structure.sites_cart().deep_copy()
    if (type(selection).__name__ == 'bool'):
      assert (self.selection.count(True) > 0)
      self.iselection = self.selection.iselection()
    else : # actually an iselection
      assert (len(self.selection) > 0)
      self.iselection = self.selection
      self.selection = flex.bool(self.sites_start.size(), False).set_selected(
        self.iselection, True)
    if (self.selection_score is None):
      self.selection_score = self.selection
    use_mp = (self.params.n_trials > 1)
    assert (params.partial_occupancy is not None)
    if (len(params.partial_occupancy) > 1):
      use_mp = True
    if (not use_mp):
      self._trials = [ self.run_trial(params.partial_occupancy[0]) ]
    else :
      self.out = null_out()
      args = []
      for occ in params.partial_occupancy :
        assert (occ > 0) and (occ <= 1)
        for k in range(self.params.n_trials):
          args.append(occ)
      self._trials = easy_mp.pool_map(
        fixed_func=self.run_trial,
        args=args,
        processes=self.nproc)

  def get_trials(self):
    """
    Return a list of results, unfiltered.
    """
    return self._trials

  def run_trial(self, occ):
    """
    Actually run the refinement - called in parallel via easy_mp.pool_map
    """
    from scitbx.array_family import flex
    seed = int(time.time() / os.getpid())
    random.seed(seed)
    flex.set_random_seed(seed)
    fmodel = self.fmodel
    selection = self.selection.deep_copy()
    iselection = self.iselection.deep_copy()
    d_min = fmodel.f_obs().d_min()
    xrs = fmodel.xray_structure.deep_copy_scatterers()
    hd_sel = xrs.hd_selection()
    sites_start = xrs.sites_cart().deep_copy()
    assert (len(selection) == len(sites_start))
    fmodel.update_xray_structure(xrs, update_f_calc=True)
    sites_start_selected = sites_start.select(iselection)
    occ_start = xrs.scatterers().extract_occupancies().deep_copy()
    u_start = xrs.extract_u_cart_plus_u_iso()
    if (self.params.adjust_b_factors_if_poor_density):
      fofc_map = fmodel.map_coefficients(
        map_type="mFo-DFc",
        exclude_free_r_reflections=True).fft_map(
          resolution_factor=1/4.).apply_sigma_scaling().real_map_unpadded()
      sites_frac = xrs.sites_frac()
      b_scale_isel = flex.size_t()
      for i_seq in iselection :
        map_value = fofc_map.tricubic_interpolation(sites_frac[i_seq])
        if (map_value < -2.5):
          b_scale_isel.append(i_seq)
      if (len(b_scale_isel) > 0):
        b_scale_sel = flex.bool(sites_frac.size(), False).set_selected(
          b_scale_isel, True)
        xrs.scale_adp(factor=0.75, selection=b_scale_sel)
    nearby_water_selection = mmtbx.building.get_nearby_water_selection(
      pdb_hierarchy=self.pdb_hierarchy,
      xray_structure=xrs,
      selection=selection)
    two_fofc_map, fofc_map = alt_confs.get_partial_omit_map(
      fmodel=fmodel,
      selection=iselection,
      selection_delete=None,#nearby_water_selection,
      negate_surrounding=self.params.negate_surrounding,
      partial_occupancy=occ)
    #xrs.set_occupancies(occ_start)
    make_sub_header("Simulated annealing into mFo-DFc map", out=self.out)
    if (self.params.shake_sites is not None):
      xrs.shake_sites_in_place(self.params.shake_sites, selection=selection)
    sites_new = mmtbx.building.run_real_space_annealing(
      xray_structure=xrs,
      pdb_hierarchy=self.pdb_hierarchy,
      processed_pdb_file=self.processed_pdb_file,
      selection=selection,
      target_map=fofc_map,
      d_min=d_min,
      params=self.params.simulated_annealing,
      #wc=5, # FIXME why does this need to be scaled?
      target_map_rsr=two_fofc_map,
      rsr_after_anneal=self.params.rsr_after_anneal,
      out=self.out,
      debug=True)
    # now compute CC of refined sites to difference map
    fmodel.update_xray_structure(xrs, update_f_calc=True)
    fc_coeffs = fmodel.f_model()
    fc_fft_map = fc_coeffs.fft_map(resolution_factor=1/4.)
    fc_map = fc_fft_map.apply_sigma_scaling().real_map_unpadded()
    pdb_atoms = self.pdb_hierarchy.atoms()
    # XXX should this only calculate statistics for the central atoms?
    map_stats = mmtbx.building.get_model_map_stats(
      selection=self.selection_score,
      target_map=fofc_map,
      model_map=fc_map,
      unit_cell=xrs.unit_cell(),
      sites_cart=sites_new,
      pdb_atoms=pdb_atoms,
      local_sampling=False)
    # reset xray structure
    xrs.set_sites_cart(sites_start)
    xrs.set_u_cart(u_start)
    fmodel.update_xray_structure(xrs, update_f_calc=True)
    # we may only want the rmsd and max. dev. from a subset of sites, e.g.
    # the central residue of a sliding window (minus hydrogens)
    selection_score = self.selection_score.deep_copy()
    if (type(selection_score).__name__ != 'bool'):
      selection_score = flex.bool(hd_sel.size(), False).set_selected(
        self.selection_score, True)
    selection_score &= ~hd_sel
    site_stats = alt_confs.coord_stats_with_flips(
      sites1=sites_start.select(selection_score),
      sites2=sites_new.select(selection_score),
      atoms=self.pdb_hierarchy.atoms().select(selection_score))
    return alt_confs.trial_result(
      sites_cart=sites_new.select(self.iselection),
      min_fofc=map_stats.min,
      mean_fofc=map_stats.mean,
      rmsd=site_stats.rmsd,
      max_dev=site_stats.max_dev,
      cc=map_stats.cc)

  def get_filtered_trials(self, include_pdb_hierarchies=False, log=None):
    """
    Filter the results based on map statistics and rotamer scoring.

    :returns: a list of acceptable trial objects, or if include_pdb_hierarchies
              is True, a list of tuples of each trial and the modified pdb
              hierarchy.
    """
    if (log is None) : log = null_out()
    trials = sorted(self.get_trials(), lambda a,b: cmp(b.cc, a.cc))
    filtered = []
    for k, trial in enumerate(trials):
      hierarchy = self.pdb_hierarchy.deep_copy()
      if (trial.min_fofc < self.params.min_fofc) or (trial.cc < self.params.cc_min):
        print("  discarding trial %d [poor map quality]:" % (k+1), file=log)
        trial.show_summary(out=log, prefix="    ")
        continue
      sites = self.sites_start.deep_copy()
      sites.set_selected(self.iselection, trial.sites_cart)
      hierarchy.atoms().set_xyz(sites)
      if (self.params.filter_rotamer_outliers):
        n_outliers = alt_confs.score_rotamers(hierarchy=hierarchy,
          selection=self.iselection)
        if (n_outliers > 0):
          print("  discarding trial %d [%d rotamer outlier(s)]:" % \
            (k+1, n_outliers), file=log)
          trial.show_summary(out=log, prefix="    ")
          continue
      if (include_pdb_hierarchies):
        filtered.append((trial, hierarchy))
      else :
        filtered.append(trial)
    return filtered

  def as_pdb_ensemble(self, log=None):
    """
    Create a multi-MselfDEL PDB hierarchy with the trial results after filtering
    by density and rotamer score.

    :returns: the new multi-model PDB hierarchy, or if no models passed the
              filtering stage, None
    """
    if (log is None) : log = null_out()
    import iotbx.pdb.hierarchy
    root = iotbx.pdb.hierarchy.root()
    trials_and_models = self.get_filtered_trials(include_pdb_hierarchies=True,
      log=log)
    n_kept = 0
    for k, (trial, hierarchy) in enumerate(trials_and_models):
      n_kept += 1
      new_model = hierarchy.only_model().detached_copy()
      new_model.id = str(n_kept)
      print("MODEL %d:" % (k+1), file=log)
      trial.show_summary(prefix="  ", out=log)
      root.append_model(new_model)
    if (n_kept == 0):
      return None
    else :
      return root
