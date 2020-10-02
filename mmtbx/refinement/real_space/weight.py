from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import mmtbx.refinement.real_space.individual_sites
import random
from six.moves import range

class run(object):
  def __init__(
      self,
      map_data,
      xray_structure,
      pdb_hierarchy,
      geometry_restraints_manager,
      gradients_method="fd",
      ncs_groups=None,
      rms_bonds_limit=0.015,
      rms_angles_limit=2.0,
      real_space_gradients_delta=1./4,
      max_iterations = 100,
      range_size=10,
      n_ranges=10,
      default_weight=50):
    """
Fast determination of optimal data/restraints weight for real-space refinement
of individual sites.
    """
    self.msg_strings = []
    # split chains into chunks
    result = []
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        if(chain.is_protein() or chain.is_na()):
          residue_range_sel = flex.size_t()
          cntr = 0
          for rg in chain.residue_groups():
            i_seqs = rg.atoms().extract_i_seq()
            cntr += 1
            if(cntr<10):
              residue_range_sel.extend(i_seqs)
            else:
              result.append(residue_range_sel)
              residue_range_sel = flex.size_t()
              residue_range_sel.extend(i_seqs)
              cntr = 0
          if(len(result)==0):
            assert residue_range_sel.size()>0
            result.append(residue_range_sel)
    self.msg_strings.append("number of chunks: %d"%len(result))
    # randomly pick chunks
    random_chunks = []
    if(len(result)>0):
      if len(result) <= n_ranges:
        # just try them all, no need to randomize
        random_chunks = list(range(len(result)))
      else:
        while len(random_chunks) <= n_ranges:
          # Put only unique choices until got enough lenght.
          # Could be slightly slow when len(random_chunks) slightly > n_ranges
          rc = random.choice(range(len(result)))
          if rc not in random_chunks:
            random_chunks.append(rc)
    self.msg_strings.append("random chunks:"%random_chunks)
    # setup refinery
    xrs_dc = xray_structure.deep_copy_scatterers()
    sel_all = flex.bool(xrs_dc.scatterers().size(), True)
    grm_dc = geometry_restraints_manager.select(sel_all)
    ro = mmtbx.refinement.real_space.individual_sites.box_refinement_manager(
      xray_structure              = xrs_dc,
      target_map                  = map_data,
      geometry_restraints_manager = grm_dc.geometry,
      real_space_gradients_delta  = real_space_gradients_delta,
      max_iterations              = max_iterations,
      ncs_groups                  = ncs_groups,
      gradients_method            = gradients_method)
    optimal_weights = flex.double()
    # loop over chunks: determine best weight for each chunk
    if(len(result)==0):
      random_chunks = [None]
    for chunk in random_chunks:
      if(chunk is None): sel = flex.bool(xrs_dc.scatterers().size(), True)
      else:
        sel = result[chunk]
        sel = flex.bool(xrs_dc.scatterers().size(), sel)
      ro.refine(
        selection        = sel,
        rms_bonds_limit  = rms_bonds_limit,
        rms_angles_limit = rms_angles_limit)
      self.msg_strings.append("chunk %s optimal weight: %9.4f"%(
          str(chunk), ro.weight_optimal))
      if(ro.weight_optimal is not None):
        optimal_weights.append(ro.weight_optimal)
    # select overall best weight
    sel = flex.sort_permutation(optimal_weights)
    optimal_weights = optimal_weights.select(sel)
    self.weight = flex.mean_default(
      optimal_weights[:optimal_weights.size()//2], default_weight)
    #mean = flex.mean(optimal_weights)
    #sel  = optimal_weights < mean*3
    #sel &= optimal_weights > mean/3
    #if(sel.count(True)>0):
    #  optimal_weights = optimal_weights.select(sel)
    #self.weight = flex.mean_default(optimal_weights, default_weight)
    self.msg_strings.append("overall best weight: %9.4f"%self.weight)

  def show(self, log, prefix=""):
    for m in self.msg_strings:
      print("%s %s"%(prefix, m), file=log)
