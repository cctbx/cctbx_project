from __future__ import division
from cctbx.array_family import flex
import mmtbx.refinement.real_space.individual_sites
import random

def run(
      map_data,
      xray_structure,
      pdb_hierarchy,
      geometry_restraints_manager,
      rms_bonds_limit=0.015,
      rms_angles_limit=2.0,
      real_space_gradients_delta=1./4,
      max_iterations = 100,
      range_size=10,
      n_ranges=10,
      log=None,
      prefix=""):
  """
Fast determination of optimal data/restraints weight for real-space refinement
of individual sites.
  """
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
  if(log is not None):
    print >> log, "%s number of chunks: %d"%(prefix, len(result))
  # randomly pick chunks
  random_chunks = []
  for i in xrange(n_ranges):
    random_chunks.append(random.choice(xrange(len(result))))
  if(log is not None):
    print >> log, "%s random chunks:"%prefix, random_chunks
  # setup refinery
  xrs_dc = xray_structure.deep_copy_scatterers()
  sel_all = flex.bool(xrs_dc.scatterers().size(), True)
  grm_dc = geometry_restraints_manager.select(sel_all)
  ro = mmtbx.refinement.real_space.individual_sites.box_refinement_manager(
    xray_structure              = xrs_dc,
    target_map                  = map_data,
    geometry_restraints_manager = grm_dc.geometry,
    real_space_gradients_delta  = real_space_gradients_delta,
    max_iterations              = max_iterations)
  optimal_weights = flex.double()
  # loop over chunks: determine best weight for each chunk
  for chunk in random_chunks:
    sel = result[chunk]
    sel = flex.bool(xrs_dc.scatterers().size(), sel)
    ro.refine(
      selection        = sel,
      rms_bonds_limit  = rms_bonds_limit,
      rms_angles_limit = rms_angles_limit)
    if(log is not None):
      print >> log,"%s chunk %3d optimal weight: %9.4f"%(
        prefix,chunk,ro.weight_optimal)
    if(ro.weight_optimal is not None):
      optimal_weights.append(ro.weight_optimal)
  # select overall best weight
  mean = flex.mean(optimal_weights)
  sel  = optimal_weights < mean*3
  sel &= optimal_weights > mean/3
  if(sel.count(True)>0):
    optimal_weights = optimal_weights.select(sel)
  weight = flex.mean(optimal_weights)
  if(log is not None):
    print >> log, "%s overall best weight: %9.4f"%(prefix, weight)
  return weight
