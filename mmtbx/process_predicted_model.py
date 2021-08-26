from __future__ import division, print_function
import sys

################################################################################
#################### TOOLS FOR PROCESSING A PREDICTED MODEL ####################
################################################################################

"""
   process_predicted_model: tools to update B-values used in some
    model predictions as an error estimate indicator and to split up model
    into domains that contain the more confident predictions
"""

from scitbx.array_family import flex
from libtbx.utils import Sorry
from scitbx.matrix import col
from libtbx import group_args
from cctbx.maptbx.segment_and_split_map import get_co

################################################################################
####################   process_predicted_model  ################################
################################################################################

def process_predicted_model(model,
  b_value_field_is = 'lddt',
  remove_low_confidence_residues = None,
  minimum_lddt = None,
  maximum_rmsd = 1.5,
  input_lddt_is_fractional = None,
  split_model_by_compact_regions = None,
  chain_id = None,
  log = sys.stdout):

  """
  process_predicted_model:
  Purpose:  Convert values in B-value field to pseudo-B-values, remove
    low_confidence residues, optionally split into compact regions.
  Rationale: predicted models may have regions of low and high confidence.
    This routine uses values in the B-value field to identify confidence,
    removes low-confidence regions, and then examines the remaining model to
    find regions that are compact (residues have high contact with neighbors)
    and that are separate from other regions (low contact with neigbors).

  Inputs:
    model:  iotbx.model.model object containing model information.
           Normally contains a single chain.   If multiple chains, process
           each separately.

    b_value_field_is:  'lddt' or 'rmsd'.  For AlphaFold models
                        the b-value field is a value of LDDT (confidence)
                        on scale of 0-1 or 0-100
                        For RoseTTAFold, the B-value field is rmsd (A)

    input_lddt_is_fractional:  if True, input lddt is scale of 0 to 1,
        otherwise 0 - 100
       If None, set to True if all lddt are from 0 to 1
    remove_low_confidence_residues: remove residues with low confidence
        (lddt or rmsd as set below)
    minimum_lddt: minimum lddt to keep residues (on same scale as b_value_field,
      if not set, calculated from maximum_rmsd).
    maximum_rmsd: alternative specification of minimum confidence based on rmsd.
        If not set, calculated from minimum_lddt. Default is 1.5 A
        maximum_rmsd or minimum_lddt required.
    split_model_by_compact_regions: split resulting model into compact regions
    chain_id: if model contains more than one chain, split this chain only.
              NOTE: only one chain can be processed at a time.

  Output:
    processed_model_info: group_args object containing:
      processed_model:  single model with regions identified in SEGID field
      processed_model_list:  list of iotbx.model.model objects, one for
                              each compact region

  """

  # Make sure we have what we expect:
  import mmtbx.model
  assert isinstance(model, mmtbx.model.manager)
  assert b_value_field_is in ('lddt','rmsd')

  # Decide what to do

  # Determine if input lddt is fractional and get b values

  b_value_field = model.get_hierarchy().atoms().extract_b()
  if b_value_field_is == 'lddt':
    if input_lddt_is_fractional is None:
      sel = (b_value_field < 0) | (b_value_field > 1)
      input_lddt_is_fractional = (sel.count(True) == 0)

    b_values = get_b_values_from_lddt(b_value_field,
       input_lddt_is_fractional = input_lddt_is_fractional)

    if input_lddt_is_fractional:
      print("B-value field interpreted as LDDT %s" %("(0 - 1)"), file = log)
    else:
      print("B-value field interpreted as LDDT %s" %("(0 - 100)"), file = log)

  elif b_value_field_is == 'rmsd':
    b_values = get_b_values_rmsd(b_value_field)
    print("B-value field interpreted as rmsd %s" %("(0 - 1)"), file = log)
  else:
    raise AssertionError("Please set b_value_field_is to either lddt or rmsd")

  if input_lddt_is_fractional:
    if minimum_lddt is not None: # convert to fractional
      minimum_lddt = minimum_lddt * 0.01
      print("Minimum LDDT converted to %.2f" %(minimum_lddt), file = log)

  # From here on we work only with fractional lddt

  # Get confidence cutoff if needed
  if remove_low_confidence_residues:
    maximum_b_value = get_cutoff_b_value(
      maximum_rmsd,
      minimum_lddt,
      log = log)
  else:
    maximum_b_value = None


  # Make a new model with new B-values

  ph  = model.get_hierarchy().deep_copy()
  ph.atoms().set_b(b_values)

  # Remove low_confidence regions if desired
  if remove_low_confidence_residues:
    n_before = ph.overall_counts().n_residues
    selection_string = " (bfactor < %s)" %maximum_b_value
    asc1 = ph.atom_selection_cache()
    sel = asc1.selection(selection_string)
    ph = ph.select(sel)
    n_after = ph.overall_counts().n_residues
    print("Total of %s of %s residues kept after B-factor filtering" %(
       n_after, n_before), file = log)

  # Get a new model
  new_model = model.as_map_model_manager().model_from_hierarchy(
     ph, return_as_model = True)

  # Get high-confidence regions as domains if desired:
  if split_model_by_compact_regions:
    # Make sure we have just 1 chain or a chain ID supplied
    chain_id = get_chain_id(model, chain_id, log = log)

    info = split_model_into_compact_units(new_model, log = log)
    if info is None:
      print("No compact regions identified", file = log)
      segid_list = []
      model_list = []
    else:
      new_model = info.model
      segid_list = info.segid_list
      print("Total of %s regions identified" %(
        len(segid_list)), file = log)
      model_list = split_model_by_segid(new_model, segid_list)
  else:
    model_list = []
    segid_list = []

  return group_args(
    group_args_type = 'processed predicted model',
    model = new_model,
    model_list = model_list,
    segid_list = segid_list
    )


def split_model_by_segid(m, segid_list):
  """
   Split a model into pieces based on SEGID
  """
  split_model_list = []
  for segid in segid_list:
    selection_string = "segid %s" %(segid)
    ph = m.get_hierarchy()
    asc1 = ph.atom_selection_cache()
    sel = asc1.selection(selection_string)
    m1 = m.select(sel)
    split_model_list.append(m1)
  return split_model_list


def get_chain_id(model, chain_id, log = sys.stdout):
  from mmtbx.secondary_structure.find_ss_from_ca import get_chain_ids
  chain_ids = get_chain_ids(model.get_hierarchy())
  if len(chain_ids) == 1:
    chain_id = chain_ids[0]
  elif chain_id is not None:
    if chain_id in chain_ids:
      print("Working on chain '%s' only" %(chain_id), file = log)
    else:
      raise Sorry("Chain %s is not in the chains present (%s)" %(
        chain_id, " ".join(chain_ids)))
  else:
    chain_id = chain_ids[0]
    print("Working on chain '%s' only (skipping other chains: %s)" %(
      chain_id, " ".join(chain_ids[1:])), file = log)

def get_cutoff_b_value(
    maximum_rmsd,
    minimum_lddt,
    log = sys.stdout):

  # Get B-value cutoff

  if maximum_rmsd is not None:
    print("Maximum rmsd of %.2f A used" %(maximum_rmsd), file = log)
  elif minimum_lddt:
    print("Minimum confidence level is %.2f" %(
      minimum_lddt), file = log)
    if minimum_lddt< 0 or \
      minimum_lddt> 1:
      raise Sorry("minimum_lddt must "+
         "be between 0 and 1")
    maximum_rmsd = get_rmsd_from_lddt(
           flex.double(1,minimum_lddt),
           is_fractional = True)[0]
    print("Maximum rmsd set to %.2f A based on confidence cutoff of %.2f" %(
       maximum_rmsd, minimum_lddt), file = log)
  else:
     raise Sorry( "Need to set either maximum_rmsd or " +
          "minimum_lddt")

  maximum_b_value = get_b_values_rmsd(
     flex.double(1,maximum_rmsd))[0]

  print("Maximum B-value to be included: %.2f A**2" %(maximum_b_value),
    file = log)
  return maximum_b_value




################################################################################
####################   get_b_values_from_lddt  ################################
################################################################################

def get_b_values_from_lddt(lddt_values,
    input_lddt_is_fractional = True):
  """
  get_b_values_from_lddt:
  Purpose:  AlphaFold models are supplied with values of LDDT (predicted
   local-distance difference test) in the B-value field.  This routine
   uses the formula from:

     Hiranuma, N., Park, H., Baek, M. et al. Improved protein structure
       refinement guided by deep learning based accuracy estimation.
       Nat Commun 12, 1340 (2021).
       https://doi.org/10.1038/s41467-021-21511-x

   to convert these values to error estimates,
   and then uses the relation between B-values and coordinate rms to
   generate pseudo-B-factors

  NOTE: formulas taken from phaser_voyager implementation by
  Claudia Millan and Massimo Sammito at
  phaser_voyager/src/Voyager/MDSLibraries/pdb_structure.py


  Inputs:
    lddt_values: flex array of lddt values
    input_lddt_is_fractional: if False, convert by multiplying * 0.01
  Outputs:
    flex array of B-values
  """

  if input_lddt_is_fractional:
    rmsd = get_rmsd_from_lddt(lddt_values) # usual
  else:
    rmsd = get_rmsd_from_lddt(0.01 * lddt_values)
  b_values = get_b_values_rmsd(rmsd)

  return b_values

################################################################################
####################   get_rmsd_from_lddt  ################################
################################################################################

def get_rmsd_from_lddt(lddt_values, is_fractional = None):
  """
  get_rmsd_from_lddt:
  Purpose:  AlphaFold models are supplied with values of LDDT (confidence)
   in the B-value field.  This routine uses the formula
   from the alphafold paper to convert these values to error estimates.
  NOTE: lddt_values can be fractional (0 to 1) or percentage (0 to 100)
  If is_fractional is not set, assume fractional if all between 0 and 1
  rmsd_est = 1.5 * flex.exp(4*(0.7-fractional_values))

  NOTE: formulas taken from phaser_voyager implementation by
  Claudia Millan and Massimo Sammito at
  phaser_voyager/src/Voyager/MDSLibraries/pdb_structure.py


  Inputs:
    lddt_values: flex array of lddt values
  Outputs:
    flex array of error estimates (A)
  """

  if is_fractional is None:
    is_fractional = ( lddt_values.min_max_mean().min >= 0  and
      lddt_values.min_max_mean().max <= 1 )


  if is_fractional:
    fractional_values = lddt_values.deep_copy()
  else:
    fractional_values = lddt_values * 0.01

  fractional_values.set_selected((fractional_values < 0), 0)
  fractional_values.set_selected((fractional_values > 1), 1)

  rmsd_est = 1.5 * flex.exp(4*(0.7-fractional_values))
  return rmsd_est

################################################################################
####################   get_b_values_rmsd #######################
################################################################################

def get_b_values_rmsd(rmsd):
  """
  get_b_values_rmsd:
  Purpose:  TTAFold models are supplied with values of rmsd (A)
   in the B-value field.  This routine converts error estimates into
   b-values

  NOTE: formulas taken from phaser_voyager implementation by
  Claudia Millan and Massimo Sammito at
  phaser_voyager/src/Voyager/MDSLibraries/pdb_structure.py


  Inputs:
    rmsd: flex array of error estimates (A)
  Outputs:
    flex array of B-values (A**2)
  """

  rmsd = rmsd.deep_copy() # do not change original

  # Make sure error estimates are in reasonable range
  rmsd.set_selected((rmsd < 0), 0)
  rmsd.set_selected((rmsd > 20), 20)

  b_values = flex.pow2(rmsd) * ((8 * (3.14159 ** 2)) / 3.0)
  return b_values


################################################################################
####################   split_model_into_compact_units   ########################
################################################################################

def split_model_into_compact_units(
     m,
     d_min = 10,
     grid_resolution = 6,
     close_distance = 15,
     minimum_region_size = 10,
     log = sys.stdout):

  """
   Function to identify groups of atoms in a model that form compact units
   (domains).  Normally used after trimming low-confidence regions in
   predicted models to isolate domains that are likely to have indeterminate
   relationships.

   Method: calculate a low-resolution map based on the input model; identify
    large blobs corresponding to domains.  Assign all atoms in structure to
    a domain.  Then regroup in order to have few cases where small parts of
    model are part of one domain but neighboring parts are part of another.

   Inputs:
   m:  cctbx.model.model object containing information about the input model
   d_min:  resolution used for low-res map
   grid_resolution:  resolution of map used to define the gridding
   close_distance:  distance between two CA (or P) atoms considered close
                    NOTE: may be useful to double default for P compared to CA
   minimum_region_size: typical size (CA or P) of the smallest segments to keep
   bfactor_min: smallest bfactor for atoms to include in calculations
   bfactor_max: largest bfactor for atoms to include in calculations


   Output:
   group_args object with members:
    m:  new model with SEGID values from 0 to N where there are N domains
      SEGID value of 0 is all atoms not included in the calculations
      SEGID 1 to N are the N domains, roughly in order along the chain
    segid_list:  list of all the SEGID values (0 is last)

   On failure:  returns None
  """


  # Make sure the model has crystal_symmetry.  Just put a box around it if nec
  m.add_crystal_symmetry_if_necessary()

  # Select CA and P atoms with B-values in range
  selection_string = '(name ca or name p)'
  new_m = m.apply_selection_string(selection_string)

  # Put the model inside a box and get a map_model_manager
  put_model_inside_cell(new_m, grid_resolution)
  mmm = new_m.as_map_model_manager()

  # Generate map at medium_res for this model and use it to get domains

  # First generate a map to get the gridding (not the best way but quick)
  mmm.set_resolution(grid_resolution)
  mmm.generate_map()

  # Now redo it at low res for the real map to use
  mmm.generate_map(d_min, map_id = 'map_manager')

  # Box the map and set SD to 1 mean to 0
  box_mmm = mmm.extract_all_maps_around_model()
  box_mmm.map_manager().set_mean_zero_sd_one()

  # Now get regions where there is model
  map_data = box_mmm.map_manager().map_data()

  #  Get a connectivity analysis of this map data
  co_info = get_best_co(map_data)
  if not co_info:
    return None # failed

  # Assign all points in box to a grouping
  co_info = assign_all_points(co_info, map_data, log = log)

  #  Assign all CA in model to a region
  regions_list = assign_ca_to_region(co_info, new_m, minimum_region_size,
     close_distance,  log = log)

  # Set segid based on regions_list

  atoms = new_m.get_hierarchy().atoms()  # new
  region_name_dict = {}
  used_regions = []
  i = 0
  segid_list = []
  for region in regions_list:
    if not region in used_regions:
      used_regions.append(region)
      i += 1
      region_name_dict[region] = i
      segid_list.append("%s" %(i))

  region_dict = {}
  for at, region in zip(atoms, regions_list):
    resseq_int = at.parent().parent().resseq_as_int()
    region_dict[resseq_int] = region_name_dict[region]

  # And apply to full model
  for at in m.get_hierarchy().atoms():
    resseq_int = at.parent().parent().resseq_as_int()
    region = region_dict.get(resseq_int,0)
    at.set_segid("%s"  %(region))

  # All done
  return group_args(
    group_args_type = 'model_info',
    model = m,
    segid_list = segid_list)

def assign_ca_to_region(co_info, m, minimum_region_size, close_distance,
     log = sys.stdout):
  region_id_map = co_info.region_id_map
  id_list = co_info.id_list
  regions_list = flex.int()
  sites_frac = m.crystal_symmetry().unit_cell().fractionalize(m.get_sites_cart())
  for sf in sites_frac:
    regions_list.append(int(region_id_map.value_at_closest_grid_point(sf)))
  # Now remove occasional ones out of place
  for cycle in range(10):
    regions_list = replace_lone_sites(regions_list)
    regions_list = replace_short_segments(regions_list, minimum_region_size)
    for i in range(len(get_unique_values(regions_list))):
     new_regions_list = merge_close_regions(
        m.get_sites_cart(), regions_list, minimum_region_size, close_distance)
     if new_regions_list:
       regions_list = new_regions_list
     else:
       break
  return regions_list

def get_unique_values(regions_list):
  unique_values = []
  for x in regions_list:
    if not x in unique_values:
      unique_values.append(x)
  return unique_values

def merge_close_regions(sites_cart, regions_list, minimum_region_size,
    close_distance = None):

  # Count number of residues in each pair that are close to the other
  # Split a group if some residues are close to other and not to self


  sites_dict = {}
  index_dict = {}
  id_list = get_unique_values(regions_list)

  for co_id in id_list:
    sel = (regions_list == co_id)
    sites_dict[co_id] = sites_cart.select(sel)
    index_dict[co_id] = sel.iselection()

  n_close_list = []
  typical_n_close = 0
  typical_n_close_n = 0
  close_to_other_list = []
  found_something = None
  for i in id_list:
    for j in id_list:
      if i==j: continue
      n_close = 0 # number in i close to j

      for k in range(sites_dict[i].size()):
         index = index_dict[i][k]
         distances = (sites_dict[j] - col(sites_dict[i][k])).norms()
         local_n_close = (distances < close_distance).count(True)
         if local_n_close > 0:
           n_close += 1
         distances_self = (sites_dict[i] - col(sites_dict[i][k])).norms()
         self_local_n_close = (distances_self < close_distance).count(True) - 1
         if local_n_close > self_local_n_close: #ZZ + minimum_region_size/2:
           close_to_other_list.append(
             group_args(group_args_type = 'closer to other',
             excess = local_n_close - self_local_n_close,
             index = index,
             i = i,
             k = k,
             j = j))


      n_close_list.append(group_args(  # how many in i close to j
        group_args_type = 'n close',
        n_close = n_close,
        i = i,
        j = j,))

  closer_to_other_swaps = get_closer_to_other(close_to_other_list,
      minimum_region_size)
  # Apply close swaps
  for s in closer_to_other_swaps:
    c = s.k_list[0]
    for k in range(c.start, c.end+1):
      regions_list[k] = s.j
      found_something = True

  update_regions_list(regions_list)

  if found_something:
    return regions_list

def get_closer_to_other(close_to_other_list, minimum_region_size):
  close_dict = {}
  for c in close_to_other_list:
    i,j = c.i,c.j
    if not i in close_dict.keys():
       close_dict[i] = {}
    if not j in close_dict[i].keys():
      close_dict[i][j] = 0
    close_dict[i][j] += 1
  for i in close_dict.keys():
    for j in close_dict[i].keys():
      if close_dict[i][j] >= 0: # minimum_region_size//2:
        pass
      else:
        del close_dict[i][j]
        if not close_dict[i]:
          del close_dict[i]
  all_k_list = []
  for i in close_dict.keys():
    for j in close_dict[i].keys():
      k_list = get_k_list(i,j,close_to_other_list)
      k_list = merge_k_list(k_list, minimum_region_size)
      if not k_list:continue
      all_k_list.append(group_args(
        group_args_type = 'k list',
        i = i,
        j = j,
        k_list = k_list,
       ))
  return all_k_list


def merge_k_list(k_list, minimum_region_size):
  n = len(k_list)
  for i in range(n):
    last_n = len(k_list)
    k_list = merge_k_list_once(k_list)
    if len(k_list) == last_n:
      break

  new_k_list = []
  for k1 in k_list:
    n1 = k1.end - k1.start + 1
    if n1 >= minimum_region_size//3:
      new_k_list.append(k1)
  return new_k_list

def merge_k_list_once(k_list):
  new_k_list = []
  for k1,k2 in zip(k_list,k_list[1:]):
    n1 = k1.end - k1.start + 1
    n2 = k2.end - k2.start + 1
    n_between = k2.start - k1.end - 1
    if n_between < min(n1,n2):
      k1.end = k2.end
      k2.start  = None
      k2.end = None
      break
  for k1 in k_list:
    if k1.start is not None:
      new_k_list.append(k1)
  return new_k_list


def get_k_list(i,j,close_to_other_list):
  k_list = []
  for c in close_to_other_list:
    if c.i==i and c.j==j:
      k_list.append(c.index)
  k_list.sort()
  k_list_as_groups = get_indices_as_ranges(k_list)
  return k_list_as_groups

def replace_short_segments(regions_list, minimum_region_size):
  id_list = get_unique_values(regions_list)
  new_regions_list = regions_list.deep_copy()
  for co_id in id_list:
    indices = (regions_list == co_id).iselection()
    indices_as_ranges = get_indices_as_ranges(indices)
    for r in indices_as_ranges:
      if r.end - r.start + 1 < minimum_region_size:
        value = regions_list[r.start - 1] if r.start > 0 else \
           regions_list[min(regions_list.size() - 1, r.end + 1)]
        for i in range(r.start,r.end+1):
          new_regions_list[i] = value

  regions_list = new_regions_list
  update_regions_list(regions_list)
  return regions_list

def update_regions_list(regions_list):
  id_list = get_unique_values(regions_list)
  id_list.sort()
  new_id_dict = {}
  i = 0
  for id_value in id_list:
    i += 1
    new_id_dict[id_value] = i
  new_id_list = list(new_id_dict.keys())
  new_id_list.sort()

  for i in range(regions_list.size()):
    regions_list[i] = new_id_dict[regions_list[i]]

def get_indices_as_ranges(indices):
  ranges = []
  grouping = None
  for index in indices:
    if not grouping or index != grouping.end + 1: # new grouping
      grouping = group_args(
        group_args_type = 'grouping',
        start = index,
        end = index)
      ranges.append(grouping)
    else:
      grouping.end = index
  return ranges

def replace_lone_sites(regions_list):
  regions_list[0] = regions_list[1]
  regions_list[-1] = regions_list[-2]
  o0 = regions_list[:-2]
  o1 = regions_list[1:-1]
  o2 = regions_list[2:]
  # find o1 is different than 0 or 2 and 0 and 2 are the same
  same_02 = (o2 == o0)
  different_01 = (o1 != 00)
  lone = (same_02 & different_01)
  for i in lone.iselection():
    index = i + 1
    regions_list[i + 1] =  regions_list[i]
  update_regions_list(regions_list)
  return regions_list

def put_model_inside_cell(m, grid_resolution):
  # Put model inside cell
  sc = m.get_sites_cart()
  sc -= col(sc.min())
  sc += col((grid_resolution,grid_resolution,grid_resolution))  # inside box
  m.set_sites_cart(sc)
  return m

def assign_all_points(co_info, map_data, log = sys.stdout):
  # add shells around all co until everything is covered
  co = co_info.co
  id_list = []
  for i in range(1,len(co_info.sorted_by_volume)):
    id_list.append(co_info.sorted_by_volume[i][1])


  # Set starting points
  region_id_map = co.result()

  done = False
  for i in range(1,map_data.all()[0]):  # max possible
    if done: continue
    for co_id in id_list:
      if done: continue
      available = (region_id_map == 0)
      if available.count(True) == 0:
        done = True
        break

      bool_region_mask = co.expand_mask(
        id_to_expand = co_id, expand_size = i)
      new = (bool_region_mask & available)
      region_id_map.set_selected(new, co_id)

  co_info.region_id_map = region_id_map
  co_info.id_list = id_list
  return co_info

def get_best_co(map_data, min_cutoff = 0.5):
  max_value = map_data.as_1d().min_max_mean().max

  # Find max number of clusters in range of 0.5 to 1.0 * max
  n = 100
  max_clusters = None
  cutoff= None
  for t in range(int(min_cutoff*n),n+1):
    threshold = t * max_value/n
    co, sorted_by_volume, min_b, max_b  = get_co(
      map_data, threshold = threshold, wrapping = False)
    if ((not max_clusters) or (len(sorted_by_volume) > max_clusters)) and (
        len(sorted_by_volume) > 1 ):
     max_clusters = len(sorted_by_volume)
     cutoff = threshold
  if max_clusters is None:
    return None
  print("Clusters: %s   Threshold: %.2f " %(max_clusters, cutoff))
  co, sorted_by_volume, min_b, max_b  = get_co(
      map_data, threshold = cutoff , wrapping = False)
  return group_args(
    group_args_type = 'co info',
    co = co,
    sorted_by_volume = sorted_by_volume)

################################################################################
####################   end of split_model_into_compact_units   #################
################################################################################
