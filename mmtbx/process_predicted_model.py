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
import iotbx.phil

################################################################################
####################   process_predicted_model  ################################
################################################################################


master_phil_str = """
  process_predicted_model{

    remove_low_confidence_residues = True
      .type = bool
      .help = Remove low-confidence residues (based on minimum lddt or \
             maximum_rmsd, whichever is specified)
      .short_caption = Remove low-confidence residues
      .expert_level = 3

    split_model_by_compact_regions = True
      .type = bool
      .help = Split model into compact regions after removing \
           low-confidence residues.
      .short_caption = Split model into compact regions
      .expert_level = 3

    maximum_domains = 3
      .type = int
      .help = Maximum domains to obtain.  You can use this to merge \
                the closest domains at the end of splitting the model. Make\
                it bigger (and optionally make domain_size smaller) to \
                get more domains.
      .short_caption = Maximum domains

    domain_size = 15
      .type = float
      .help = Approximate size of domains to be found (A units).  This is the \
               resolution that \
              will be used to make a domain map.  If you are getting too many \
              domains, try making domain_size bigger (maximum is 70 A).
      .short_caption = Domain size (A)

    minimum_domain_length = 10
      .type = float
      .help = Minimum length of a domain to keep (reject at end if smaller).
      .short_caption = Minimum domain length (residues)

    maximum_fraction_close = 0.3
      .type = float
      .help = Maximum fraction of CA in one domain close to one in another \
              before merging them
      .short_caption = Maximum fraction close

    minimum_sequential_residues = 5
      .type = int
      .help = Minimum length of a short segment to keep (reject at end ).
      .short_caption = Minimum sequential_residues

    minimum_remainder_sequence_length = 15
      .type = int
      .help = used to choose whether the sequence of a removed \
               segment is written to the remainder sequence file.
      .short_caption = Minimum remainder sequence length

    b_value_field_is = *lddt rmsd b_value
      .type = choice
      .help = The B-factor field in predicted models can be LDDT \
             (confidence, 0-1 or 0-100) or rmsd (A) or a B-factor
      .short_caption = Contents of B-value field (required)
      .expert_level = 3

    input_lddt_is_fractional = None
      .type = bool
      .help = You can specify if the input lddt values (in B-factor field) \
                are fractional (0-1) or not (0-100). By default if all  \
               values are between 0 and 1 it is fractional.
      .short_caption = Input lddt is fractional

    minimum_lddt = None
      .type = float
      .help = If low-confidence residues are removed, the cutoff is defined by \
          minimum_lddt or maximum_rmsd, whichever is defined (you cannot \
          define both).  A minimum lddt of 0.70 corresponds to a maximum rmsd \
          of 1.5.  Minimum lddt values are fractional or not depending on \
          the value of input_lddt_is_fractional.
      .short_caption = Minimum lddt


    maximum_rmsd = 1.5
      .type = float
      .help = If low-confidence residues are removed, the cutoff is defined by \
          minimum_lddt or maximum_rmsd, whichever is defined (you cannot \
          define both).  A minimum lddt of 0.70 corresponds to a maximum rmsd \
          of 1.5.  Minimum lddt values are fractional or not depending on \
          the value of input_lddt_is_fractional.
      .short_caption = Maximum rmsd

    default_maximum_rmsd = 1.5
      .type = float
      .help = Default value of maximum_rmsd, used if maximum_rmsd is not set
      .short_caption = default_maximum_rmsd

    subtract_minimum_b = False
      .type = bool
      .help = If set, subtract the lowest B-value from all B-values \
          just before writing \
          out the final files.  Does not affect the cutoff for removing low-\
           confidence residues.

     pae_power = 1
       .type = float
       .help = If PAE matrix (predicted alignment error matrix) is supplied,\
            each edge in the graph will be weighted proportional to \
              (1/pae**pae_power). Use this to try and get the number of domains\
              that you want (try 1, 0.5, 1.5, 2)
       .short_caption = PAE power (if PAE matrix supplied)

     pae_cutoff = 5
       .type = float
       .help = If PAE matrix (predicted alignment error matrix) is supplied,\
            graph edges will only be created for residue pairs with \
            pae<pae_cutoff
       .short_caption = PAE cutoff (if PAE matrix supplied)

     pae_graph_resolution = 0.5
       .type = float
       .help = If PAE matrix (predicted alignment error matrix) is supplied,\
            pae_graph_resolution regulates how aggressively the clustering \
            algorithm is. Smaller values lead to larger clusters. Value \
            should be larger than zero, and values larger than 5 are \
            unlikely to be useful
       .short_caption = PAE graph resolution (if PAE matrix supplied)

     weight_by_ca_ca_distance = False
       .type = bool
       .help = Adjust the edge weighting for each residue pair according  \
             to the distance between CA residues. If this is True, \
             then distance_model can be provided, otherwise supplied model \
             will be used. See also distance_power
       .short_caption = Weight by CA-CA distance

     distance_power = 1
       .type = float
       .help = If weight_by_ca_ca_distance is True, then edge weights will \
          be multiplied by 1/distance**distance_power.
       .short_caption = Distance power (for weighting by CA-CA distance)

     stop_if_no_residues_obtained = True
      .type = bool
      .help = Raise Sorry and stop if processing yields no residues
      .short_caption = Stop if no result

     keep_all_if_no_residues_obtained = False
      .type = bool
      .help = Keep everything if processing yields no residues
      .short_caption = Keep all if no result

     vrms_from_rmsd_intercept = 0.25
       .type = float
       .help = Estimate of vrms (error in model) from pLDDT will be based on\
           vrms_from_rmsd_intercept + vrms_from_rmsd_slope * pLDDT \
           where mean pLDDT of non-low_confidence_residues is used.
       .short_caption = vRMS intercept

     vrms_from_rmsd_slope = 1.0
       .type = float
       .help = Estimate of vrms (error in model) from pLDDT will be based on\
           vrms_from_rmsd_intercept + vrms_from_rmsd_slope * pLDDT \
           where mean pLDDT of non-low_confidence_residues is used.
       .short_caption = vRMS slope

    }

    """


def process_predicted_model(
    model,
    params,
    pae_matrix = None,
    distance_model = None,
    mark_atoms_to_keep_with_occ_one = False,
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

  Inputs (supplied as model and a params object):
    model:  iotbx.model.model object containing model information.
           Normally contains a single chain.   If multiple chains, process
           each separately.

    b_value_field_is:  'lddt' or 'rmsd' or 'b_value'.  For AlphaFold models
                        the b-value field is a value of LDDT (confidence)
                        on scale of 0-1 or 0-100
                        For RoseTTAFold, the B-value field is rmsd (A)
                        If b_value... it is left as is.

    input_lddt_is_fractional:  if True, input lddt is scale of 0 to 1,
        otherwise 0 - 100
       If None, set to True if all lddt are from 0 to 1
    remove_low_confidence_residues: remove residues with low confidence
        (lddt or rmsd as set below)
    minimum_lddt: minimum lddt to keep residues (on same scale as b_value_field,
      if not set, calculated from maximum_rmsd).
    maximum_rmsd: alternative specification of minimum confidence based on rmsd.
        If not set, calculated from minimum_lddt.
    default_maximum_rmsd:  used as default if nothing specified for
         maximum_rmsd or minimum_lddt .Default is 1.5 A,
    split_model_by_compact_regions: split resulting model into compact regions
      and return a list of models in the group_arg return object
    pae_matrix:  matrix of predicted aligned errors (e.g., from AlphaFold2), NxN
      matrix of RMSD values, N = number of residues in model.
      Alternative to splitting by compact regions. Split to minimize predicted
          aligned errors in each grouping.
        pae_power (default=1): each edge in the graph will be weighted
           proportional to (1/pae**pae_power)
        pae_cutoff (optional, default=5): graph edges will only be created for
         residue pairs with pae<pae_cutoff
    weight_by_ca_ca_distance: (optional, default=False): adjust the edge
        weighting for each residue pair according to the distance between
        CA residues. If this is True, then distance_model must be provided.
    distance_power (optional, default=1): If weight_by_ca_ca_distance` is True,
         then edge weights will be multiplied by 1/distance**distance_power.
    distance_model ((optional, default=None): A model corresponding to the
        PAE matrix. Only needed if weight_by_ca_ca_distances is True.

    domain_size: typical size of domains (resolution used for filtering is
       the domain size)
    minimum_domain_length: minimum length (residues) of a domain to keep
    maximum_fraction_close: Merge domains with more than this fraction of close
                           CA atoms
    maximum_domains: if more than this many domains, merge close ones to reduce
       number
    chain_id: if model contains more than one chain, split this chain only.
              NOTE: only one chain can be processed at a time.
    if subtract_minimum_b is set, subtract minimum(B values) from all B values
       after applying any B value cutoffs

    If mark_atoms_to_keep_with_occ_one is set, return list of models, each
      of which is complete, but in which occupancy = 1 marks atoms to include
      and occupancy=0 marks those to exclude

    If stop_if_no_residues_obtained (default), stop with Sorry if no residues
      are obtained after processing, except if
        keep_all_if_no_residues_obtained (not default), then take everything.

  Output:
    processed_model_info: group_args object containing:
      processed_model:  single model with regions identified in chainid field
      model_list:  list of models representing domains
      lddt_list: one lddt on scale of 0 to 1 for each residue in input model.
      vrms_list: one vrms estimate (rms model error in A for each model)

  How to get the parameters object set up:

    You can set up a parameters object like this (see example at end of this
    file as well:

    master_phil = iotbx.phil.parse(master_phil_str)
    params = master_phil.extract()

    The default values are set in the master_phil_str string above.
    You can then set values of params:

    params.process_predicted_model.split_model_by_compact_regions = True


  """

  # Make sure we have what we expect:
  import mmtbx.model
  assert isinstance(model, mmtbx.model.manager)

  # Make sure we have just 1 chain or a chain ID supplied
  chain_ids = model.chain_ids()
  if len(chain_ids) != 1:
    chain_id = model.first_chain_id()
    model.add_crystal_symmetry_if_necessary()
    model = model.apply_selection_string('chain %s' %(chain_id))

  # Decide what to do
  p = params.process_predicted_model

  # Determine if input lddt is fractional and get b values

  b_value_field = model.get_hierarchy().atoms().extract_b()
  if p.b_value_field_is == 'lddt':
    if p.input_lddt_is_fractional is None:
      sel = (b_value_field < 0) | (b_value_field > 1)
      p.input_lddt_is_fractional = (sel.count(True) == 0)

    b_values = get_b_values_from_lddt(b_value_field,
       input_lddt_is_fractional = p.input_lddt_is_fractional)

    if p.input_lddt_is_fractional:
      print("B-value field interpreted as LDDT %s" %("(0 - 1)"), file = log)
    else:
      print("B-value field interpreted as LDDT %s" %("(0 - 100)"), file = log)

  elif p.b_value_field_is == 'rmsd':
    b_values = get_b_values_rmsd(b_value_field)
    print("B-value field interpreted as rmsd %s" %("(0 - 1)"), file = log)

  elif p.b_value_field_is == 'b_value':
    b_values = b_value_field
    print("B-value field interpreted as b_values", file = log)
  else:
    raise AssertionError("Please set b_value_field_is to either lddt or rmsd")

  if (not p.input_lddt_is_fractional):
    if p.minimum_lddt is not None: # convert to fractional
      p.minimum_lddt = p.minimum_lddt * 0.01
      print("Minimum LDDT converted to %.2f" %(p.minimum_lddt), file = log)

  # From here on we work only with fractional lddt

  # Get confidence cutoff if needed
  if p.remove_low_confidence_residues:
    maximum_b_value = get_cutoff_b_value(
      p.maximum_rmsd,
      p.minimum_lddt,
      default_maximum_rmsd = p.default_maximum_rmsd,
      log = log)
  else:
    maximum_b_value = None


  # Offset b-values and cutoff if requested
  if p.subtract_minimum_b:
    minimum_b = b_values.min_max_mean().min
    b_values -= minimum_b
    assert b_values.min_max_mean().min == 0
    if maximum_b_value is not None:
      maximum_b_value -= minimum_b  # offset this too
    print("Subtracting minimum B of " +
      "%.2f from values and from cutoff (now %s)" %(
      minimum_b, " %.2f" %maximum_b_value if maximum_b_value is not None else "None"), file = log)

  # Make a new model with new B-values

  ph  = model.get_hierarchy().deep_copy()
  ph.atoms().set_b(b_values)

  # Remove low_confidence regions if desired
  if p.remove_low_confidence_residues:
    n_before = ph.overall_counts().n_residues
    selection_string = " (bfactor < %s)" %maximum_b_value
    asc1 = ph.atom_selection_cache()
    sel = asc1.selection(selection_string)
    working_ph = ph.select(sel)
    if p.minimum_sequential_residues:  #
      # Remove any very short segments
      asc1 = working_ph.atom_selection_cache()
      sel1 = asc1.selection('name ca')
      ca_ph = working_ph.select(sel1)
      selection_to_remove = get_selection_for_short_segments(ca_ph,
         p.minimum_sequential_residues)
      if selection_to_remove:
        print("Removing short segments: %s" %(selection_to_remove), file = log)
        asc1 = ph.atom_selection_cache() # original ph
        sel2 = asc1.selection(selection_to_remove)
        sel = ~ (~sel | sel2)

    new_ph = ph.select(sel)
    n_after = new_ph.overall_counts().n_residues
    print("Total of %s of %s residues kept after B-factor filtering" %(
       n_after, n_before), file = log)
    keep_all = False
    remainder_sequence_str = None
    if n_after == 0:
      if p.stop_if_no_residues_obtained:
        raise Sorry("No residues remaining after filtering...please check if "+
         "B-value field is really '%s'. Adjust maximum_rmsd if necessary." %(
           p.b_value_field_is))
      elif p.keep_all_if_no_residues_obtained:
        keep_all = True
        print("Keeping everything as no residues obtained after filtering",
           file = log)
      else:
        return group_args(
         group_args_type = 'processed predicted model',
         model = None,
         model_list = [],
         chainid_list = [],
         remainder_sequence_str = "",
         b_values = [],
         )

    if not keep_all:
      removed_ph = ph.select(~sel)
      from mmtbx.secondary_structure.find_ss_from_ca import model_info, \
         split_model
      remainder_sequence_str = ""
      for m in split_model(model_info(removed_ph)):
        seq = m.hierarchy.as_sequence(as_string = True)
        if len(seq) >= p.minimum_remainder_sequence_length:
          remainder_sequence_str += "\n> fragment sequence "
          remainder_sequence_str += "\n%s\n" %(
            m.hierarchy.as_sequence(as_string = True))
      ph = new_ph
  else:
    remainder_sequence_str = None

  # Get a new model
  new_model = model.as_map_model_manager().model_from_hierarchy(
     ph, return_as_model = True)

  # Get high-confidence regions as domains if desired:
  if p.split_model_by_compact_regions:

    if pae_matrix is not None: # use pae matrix method
      info = split_model_with_pae(model, new_model, pae_matrix,
        maximum_domains = p.maximum_domains,
        pae_power = p.pae_power,
        pae_cutoff = p.pae_cutoff,
        pae_graph_resolution = p.pae_graph_resolution,
        minimum_domain_length = p.minimum_domain_length,
        weight_by_ca_ca_distance = p.weight_by_ca_ca_distance,
        distance_power = p.distance_power,
        distance_model = distance_model,
        log = log)
    else: # usual
      info = split_model_into_compact_units(new_model,
        d_min = p.domain_size,
        maximum_domains = p.maximum_domains,
        minimum_domain_length = p.minimum_domain_length,
        maximum_fraction_close = p.maximum_fraction_close,
        log = log)
    if info is None:
      print("No compact regions identified", file = log)
      chainid_list = []
      model_list = []
    else:
      new_model = info.model
      chainid_list = info.chainid_list
      print("Total of %s regions identified" %(
        len(chainid_list)), file = log)
      model_list = split_model_by_chainid(new_model, chainid_list,
        mark_atoms_to_keep_with_occ_one = mark_atoms_to_keep_with_occ_one)
  else:
    model_list = []
    chainid_list = []

  # Estimate vrms (model error) for each domain
  vrms_list = get_vrms_list(p, model_list)
  return group_args(
    group_args_type = 'processed predicted model',
    model = new_model,
    model_list = model_list,
    chainid_list = chainid_list,
    remainder_sequence_str = remainder_sequence_str,
    b_values = b_values,
    vrms_list = vrms_list,
    )

def get_vrms_list(p, model_list):
  vrms_list = []
  for m in model_list:
    b_values = m.apply_selection_string(
      'name ca and not element ca').get_b_iso()
    rmsd = get_rmsd_from_lddt(get_lddt_from_b(b_values)).min_max_mean().mean
    vrms = rmsd * p.vrms_from_rmsd_slope + p.vrms_from_rmsd_intercept
    vrms_list.append(vrms)
  return vrms_list

def get_selection_for_short_segments(ph, minimum_sequential_residues):
  chain_dict = {}
  for model in ph.models():
    for chain in model.chains():
      residue_list = []
      for rg in chain.residue_groups():
        resseq_int = rg.resseq_as_int()
        residue_list.append(resseq_int)
      residue_list = sorted(residue_list)
      chain_dict[chain.id] = residue_list
  selections = []
  for chain_id in chain_dict.keys():
    residue_list = chain_dict[chain_id]
    for r in get_indices_as_ranges(residue_list):
      if r.end - r.start + 1 < minimum_sequential_residues:
        selections.append("(chain %s and resseq %s:%s)" %(
          chain_id, r.start, r.end))
  selection_string = " or ".join(selections)
  return selection_string




def split_model_by_chainid(m, chainid_list,
    mark_atoms_to_keep_with_occ_one = False):
  """
   Split a model into pieces based on chainid
   Optionally write out everything for each model, using
      occupancy=0 to mark everything that is not select3ed
  """
  split_model_list = []
  for chainid in chainid_list:
    selection_string = "chain %s" %(chainid)
    ph = m.get_hierarchy()
    asc1 = ph.atom_selection_cache()
    sel = asc1.selection(selection_string)
    if (not mark_atoms_to_keep_with_occ_one): # usual
      m1 = m.select(sel)
    else:  # for Voyager, mark unused with zero occupancies
      m1 = m.deep_copy()
      ph1 = m1.get_hierarchy()
      atoms = ph1.atoms()
      occupancies = atoms.extract_occ()
      occupancies.set_selected(sel, 1)
      occupancies.set_selected(~sel, 0)
      atoms.set_occ(occupancies)
    split_model_list.append(m1)
  return split_model_list

def get_cutoff_b_value(
    maximum_rmsd,
    minimum_lddt,
    default_maximum_rmsd = None,
    log = sys.stdout):

  # Get B-value cutoff

  if maximum_rmsd is None and minimum_lddt is None:
    maximum_rmsd = default_maximum_rmsd
    assert maximum_rmsd is not None

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
####################   get_lddt_from_b          ################################
################################################################################

def get_lddt_from_b(b_values, input_lddt_is_fractional = True):
  """  Inverse of get_b_values_from_lddt
  Inputs:
    flex array of B-values
    input_lddt_is_fractional: if False, convert by multiplying * 100 at end
  Outputs:
    lddt_values: flex array of lddt values
  """
  if not b_values:
    return None

  # b_values = flex.pow2(rmsd) * ((8 * (3.14159 ** 2)) / 3.0)
  rmsd = flex.sqrt( b_values/ ((8 * (3.14159 ** 2)) / 3.0))

  # rmsd  = 1.5 * flex.exp(4*(0.7-lddt))
  lddt = 0.7 - 0.25 * flex.log(rmsd/1.5)


  if not input_lddt_is_fractional:
    lddt = lddt * 100

  return lddt


################################################################################
####################   get_rmsd_from_lddt  ################################
################################################################################

def get_rmsd_from_lddt(lddt_values, is_fractional = None):
  """
  get_rmsd_from_lddt:
  Purpose:  AlphaFold models come with predicted LDDT values in the B-value
   field to indicate confidence.  This routine uses a formula provided in the
   supplementary material of the RoseTTAFold paper to convert these values
   to error estimates.
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
####################   split_model_with_pae  ###################################
################################################################################

def split_model_with_pae(
     model,
     m,
     pae_matrix,
     maximum_domains = None,
     pae_power = 1.,
     pae_cutoff = 5.,
     pae_graph_resolution = 0.5,
     minimum_domain_length = 10,
     weight_by_ca_ca_distance = False,
     distance_power = 1,
     distance_model = None,
     log = sys.stdout):

  """
   Function to identify groups of atoms in a model that form compact units
   using a predicted alignment error matrix (pae_matrix).
   Normally used after trimming low-confidence regions in
   predicted models to isolate domains that are likely to have indeterminate
   relationships.

   m:  cctbx.model.model object containing information about the input model
     after trimming
   model: model before trimming
   pae_matrix:  matrix of predicted aligned errors (e.g., from AlphaFold2), NxN
       matrix of RMSD values, N = number of residues in model.
   maximum_domains:  If more than this many domains, merge closest ones until
     reaching this number
   pae_power (default=1): each edge in the graph will be weighted
       proportional to (1/pae**pae_power)
   pae_cutoff (optional, default=5): graph edges will only be created for
       residue pairs with pae<pae_cutoff
   pae_graph_resolution (optional, default = 0.5): regulates how aggressively
       the clustering algorithm is. Smaller values lead to larger clusters.
       Value should be larger than zero, and values larger than 5 are
        unlikely to be useful
   weight_by_ca_ca_distance: (optional, default=False): adjust the edge
        weighting for each residue pair according to the distance between
        CA residues. If this is True, then distance_model must be provided.
   distance_power (optional, default=1): If weight_by_ca_ca_distance` is True,
         then edge weights will be multiplied by 1/distance**distance_power.
   distance_model ((optional, default=None): A model corresponding to the
        PAE matrix. Only needed if weight_by_ca_ca_distances is True.
   minimum_domain_length:  if a region is smaller than this, skip completely

   Output:
   group_args object with members:
    m:  new model with chainid values from 0 to N where there are N domains
      chainid 1 to N are the N domains, roughly in order along the chain.
    chainid_list:  list of all the chainid values

   On failure:  returns None
  """

  print("\nSelecting domains with predicted alignment uncertainty estimates",
     file = log)
  # Select CA and P atoms with B-values in range
  selection_string = '(name ca or name p)'
  m_ca= m.apply_selection_string(selection_string)
  n = model.apply_selection_string(selection_string
       ).get_hierarchy().overall_counts().n_residues

  # Make sure matrix matches
  if tuple(pae_matrix.shape) != (n, n):
     raise Sorry("The pae matrix has a size of (%s,%s) " %(
      tuple(pae_matrix.shape)) +
      "but the number of residues in the model is %s" %(n))
  first_resno = model.first_resseq_as_int()

  #  Assign all CA in model to a region
  from mmtbx.domains_from_pae import get_domain_selections_from_pae_matrix
  selection_list = get_domain_selections_from_pae_matrix(
    pae_matrix = pae_matrix,
     pae_power = pae_power,
     pae_cutoff = pae_cutoff,
     graph_resolution = pae_graph_resolution,
     first_resno = first_resno,
     weight_by_ca_ca_distance = weight_by_ca_ca_distance,
     distance_power = distance_power,
     distance_model = distance_model,
    )

  # And apply to full model

  unique_regions = list(range(len(selection_list)))

  keep_list = []
  good_selections = []
  ph = m.get_hierarchy()
  for selection_string, region_number in zip(selection_list,unique_regions):
    asc1 = ph.atom_selection_cache()
    sel = asc1.selection(selection_string)
    if sel.count(True) >= minimum_domain_length:
      keep_list.append(True)
      good_selections.append(selection_string)
    else:
      keep_list.append(False)
      print("Skipping region with selection '%s' that contains %s residues" %(
         selection_string,sel.count(True)),
        file = log)

  region_name_dict, chainid_list = get_region_name_dict(m, unique_regions,
    keep_list = keep_list)
  print("\nSelection list based on PAE values:", file =log)

  # Now create new model with chains based on region list
  full_new_model = None
  for keep, selection_string, region_number in zip(
     keep_list, selection_list,unique_regions):
    if not keep: continue
    new_m = m.apply_selection_string(selection_string)
    print("%s (%s residues)  "%(selection_string,
      new_m.get_hierarchy().overall_counts().n_residues), file = log)
    # Now put all of new_m in a chain with chain.id = str(region_number)
    for model in new_m.get_hierarchy().models()[:1]: # only one model
      for chain in model.chains()[:1]: # only allowing one chain
        chain.id = region_name_dict[region_number]
    if full_new_model:
      full_new_model = add_model(full_new_model, new_m)
    else:
      full_new_model = new_m
  m = full_new_model

  # All done
  return group_args(
    group_args_type = 'model_info',
    model = m,
    chainid_list = chainid_list)

################################################################################
####################   split_model_into_compact_units   ########################
################################################################################

def split_model_into_compact_units(
     m,
     d_min = 15,
     grid_resolution = 6,
     close_distance = 15,
     minimum_domain_length = 10,
     maximum_fraction_close = 0.3,
     maximum_domains = None,
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
   d_min:  resolution used for low-res map.  Corresponds roughly to domain size.
   grid_resolution:  resolution of map used to define the gridding
   close_distance:  distance between two CA (or P) atoms considered close
                    NOTE: may be useful to double default for P compared to CA
   minimum_domain_length: typical size (CA or P) of the smallest segments to keep
   minimum_remainder_sequence_length: minimum length of a removed sequence
      segment to write out to a new sequence file
   bfactor_min: smallest bfactor for atoms to include in calculations
   bfactor_max: largest bfactor for atoms to include in calculations
   maximum_domains:  If more than this many domains, merge closest ones until
     reaching this number

   Output:
   group_args object with members:
    m:  new model with chainid values from 0 to N where there are N domains
      chainid 1 to N are the N domains, roughly in order along the chain.
    chainid_list:  list of all the chainid values

   On failure:  returns None
  """
  print("\nSelecting domains as compact chains",
     file = log)
  d_min = min(50, d_min) # limitation in fmodel

  # Make sure the model has P1 crystal_symmetry.  Put a box around it that is
  #  big (do not use original crystal symmetry because it might  be too big
  #  or too small)

  m = m.deep_copy()  # don't modify original

  box_cushion = 0.5 * d_min  # big box
  original_crystal_symmetry = m.crystal_symmetry()
  original_uc_crystal_symmetry = m.unit_cell_crystal_symmetry()
  m.add_crystal_symmetry_if_necessary(box_cushion = box_cushion, force = True)
  m.set_shift_cart(None)
  m.set_unit_cell_crystal_symmetry(m.crystal_symmetry())

  # Select CA and P atoms with B-values in range
  selection_string = '(name ca or name p)'
  m_ca= m.apply_selection_string(selection_string)

  # Put the model inside a box and get a map_model_manager
  put_model_inside_cell(m_ca, grid_resolution)
  mmm = m_ca.as_map_model_manager()

  # Generate map at medium_res for this model and use it to get domains

  print("Generating map at resolution of %.1f A to identify domains" %(
    d_min), file = log)
  mmm.set_resolution(d_min)
  mmm.generate_map(d_min, resolution_factor = 0.25 * grid_resolution/d_min )

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
  regions_list = assign_ca_to_region(co_info, m_ca, minimum_domain_length,
     close_distance,
     maximum_domains = maximum_domains,
     maximum_fraction_close = maximum_fraction_close,
     log = log)

  info = set_chain_id_by_region(m, m_ca, regions_list, log = log)
  if original_crystal_symmetry and info and info.model:
    info.model.set_crystal_symmetry(original_crystal_symmetry)
  return info

def get_region_name_dict(m, unique_regions, keep_list = None):
  region_name_dict = {}
  chainid_list = []
  chainid = m.first_chain_id().strip()
  if not keep_list:
    keep_list = len(unique_regions) * [True]
  assert len(keep_list) == len(unique_regions)

  i = 0
  for region_number, keep in zip(unique_regions, keep_list):
      if not keep: continue
      i += 1
      if len(chainid) == 1 and i < 10:
          region_name = "%s%s" %(chainid,i)
      else:
          region_name = str(i)
      region_name_dict[region_number] = region_name
      chainid_list.append(region_name)
  return region_name_dict, chainid_list

def set_chain_id_by_region(m, m_ca, regions_list, log = sys.stdout):
  # Set chainid based on regions_list

  atoms = m_ca.get_hierarchy().atoms()  # new
  unique_regions = get_unique_values(regions_list)

  region_name_dict, chainid_list = get_region_name_dict(m, unique_regions)

  region_dict = {}
  for at, region_number in zip(atoms, regions_list):
    resseq_int = at.parent().parent().resseq_as_int()
    region_dict[resseq_int] = region_number

  # And apply to full model
  full_region_list = flex.int()
  for at in m.get_hierarchy().atoms():
    resseq_int = at.parent().parent().resseq_as_int()
    region = region_dict.get(resseq_int,0)
    full_region_list.append(region)

  # Now create new model with chains based on region list
  full_new_model = None
  print("\nSelection list based on domains:", file =log)
  for region_number in unique_regions:
    sel = (full_region_list == region_number)
    new_m = m.select(sel)
    selection_string = selection_string_from_model(
       new_m.apply_selection_string("name ca or name P"))

    print("%s (%s residues)  "%(selection_string,
      new_m.get_hierarchy().overall_counts().n_residues), file = log)
    # Now put all of new_m in a chain with chain.id = str(region_number)
    for model in new_m.get_hierarchy().models()[:1]: # only one model
      for chain in model.chains()[:1]: # only allowing one chain
        chain.id = region_name_dict[region_number]
    if full_new_model:
      full_new_model = add_model(full_new_model, new_m)
    else:
      full_new_model = new_m
  full_new_model.reset_after_changing_hierarchy()
  m = full_new_model

  # All done
  return group_args(
    group_args_type = 'model_info',
    model = m,
    chainid_list = chainid_list)

def selection_string_from_model(model):
    resno_list = get_residue_numbers_in_model(model)
    from mmtbx.domains_from_pae import cluster_as_selection
    selection_string = cluster_as_selection(resno_list)
    return selection_string


def assign_ca_to_region(co_info,
    m,
    minimum_domain_length,
    close_distance,
    maximum_domains = None,
    maximum_fraction_close = None,
    n_cycles = 10,
    log = sys.stdout):
  region_id_map = co_info.region_id_map
  id_list = co_info.id_list
  regions_list = flex.int()
  sites_frac = m.crystal_symmetry().unit_cell().fractionalize(m.get_sites_cart())
  for sf in sites_frac:
    regions_list.append(int(region_id_map.value_at_closest_grid_point(sf)))
  # Now remove occasional ones out of place
  for cycle in range(n_cycles):
    regions_list = replace_lone_sites(regions_list)
    regions_list = replace_short_segments(regions_list, minimum_domain_length)
    for i in range(len(get_unique_values(regions_list))):
      new_regions_list = swap_close_regions(
        m.get_sites_cart(), regions_list, minimum_domain_length, close_distance)
      if new_regions_list:
         regions_list = new_regions_list
      else:
       break
  # Merge close regions if there are too many
  for k in range(len(get_unique_values(regions_list))):
    if maximum_domains and \
         len((get_unique_values(regions_list))) > maximum_domains:
      regions_list = merge_closest_regions(m.get_sites_cart(), regions_list,
        close_distance, log = log)
    else:
      break

  # Merge close regions if they are really close
  for k in range(len(get_unique_values(regions_list))):
    regions_list = merge_very_close_regions(m.get_sites_cart(), regions_list,
        close_distance,
        minimum_domain_length =  minimum_domain_length,
        maximum_fraction_close = maximum_fraction_close,
        log = log)

  # Finally check for any short fragments not attached to neighbors
  regions_list = remove_short_fragments_obscured_by_gap(regions_list,
    m, minimum_domain_length)
  return regions_list

def remove_short_fragments_obscured_by_gap(regions_list,
    m, minimum_domain_length):
  # Find any regions that are very short (<minimum_domain_length) and merge
  # with adjacent sequence if available
  # This catches cases where there was a gap in sequence.

  region_dict = {}
  for at, region_number in zip(m.get_hierarchy().atoms(), regions_list):
    resseq_int = at.parent().parent().resseq_as_int()
    region_dict[resseq_int] = region_number

  # Find all cases where regions go like:
  #  1 1 1 2 2 (gap)   -> 1 1 1 1 1 (gap)
  #  1 1 1 2 2 3 3 3 3  -> 1 1 1 1 1 3 3 3 3 or 1 1 1 3 3 3 3 3 3
  #  (gap) 2 2 3 3 3 3 -> (gap) 3 3 3 3 3 3
  # First split up resseq_int into ranges..
  residues_as_groups = get_indices_as_ranges(list(region_dict.keys()))
  for r in residues_as_groups:
    # Find all the places where region_number changes
    working_regions = []
    working_region = None
    for i in range(r.start, r.end+1):
      region_number = region_dict[i]
      if not working_region or region_number !=working_region.region_number:
        working_region = group_args(
          group_args_type = 'working region',
          region_number = region_number,
          start = i,
          end = i,
          )
        working_regions.append(working_region)
      else:
        working_region.end = i
    for previous_region,working_region,next_region in zip(
      [None]+working_regions[:-1], working_regions, working_regions[1:]+[None]):
      if (working_region.end - working_region.start + 1) < \
           minimum_domain_length:
        if previous_region:
          working_region.region_number = previous_region.region_number
        elif next_region:
          working_region.region_number = next_region.region_number
        else:  # skip as nothing to do
          pass
    # And update dictionary
    for working_region in working_regions:
      for i in range(working_region.start, working_region.end+1):
        region_dict[i] = working_region.region_number
  # And use dictionary to update regions_list

  new_regions_list = regions_list.deep_copy()
  i = -1
  for at in m.get_hierarchy().atoms():
    i += 1
    resseq_int = at.parent().parent().resseq_as_int()
    new_regions_list[i] = region_dict[resseq_int]

  return new_regions_list

def merge_very_close_regions(sites_cart, regions_list, close_distance,
     minimum_domain_length= None,
     maximum_fraction_close = None,
     log = sys.stdout):
  unique_values = get_unique_values(regions_list)
  n = len(unique_values)
  close_to_other_info  = get_close_to_other_list(sites_cart, regions_list,
     close_distance)
  if close_to_other_info.n_close_list:  # there are close pairs
    # Get best pair
    n_close_list = sorted(close_to_other_info.n_close_list,
      key = lambda c: c.n_close, reverse = True)
    best_pair = n_close_list[0]
    n_close = best_pair.n_close
    n_possible = best_pair.n_possible
    if n_possible >=  minimum_domain_length and \
        n_close >= n_possible * maximum_fraction_close and \
        best_pair.i is not None:
      best_i = best_pair.i
      best_j = best_pair.j
      print("Merging groups with %s sets of close residues" %(
        best_pair.n_close), file = log)

      sel = (regions_list == best_j)
      regions_list.set_selected(sel, best_i)
      update_regions_list(regions_list)

  return regions_list
def merge_closest_regions(sites_cart, regions_list, close_distance,
     log = sys.stdout):
  unique_values = get_unique_values(regions_list)
  n = len(unique_values)
  close_to_other_info  = get_close_to_other_list(sites_cart, regions_list,
     close_distance)
  if close_to_other_info.n_close_list:  # there are close pairs
    # Get best pair
    n_close_list = sorted(close_to_other_info.n_close_list,
      key = lambda c: c.n_close, reverse = True)
    best_pair = n_close_list[0]
    best_i = best_pair.i
    best_j = best_pair.j
    if best_i is not None:
      print("Merging groups with %s sets of close residues" %(
        best_pair.n_close), file = log)
  else:  # just take the pair that comes closest
    best_i = None
    best_j = None
    best_dist = None
    for k in range(len(unique_values)):
      i = unique_values[k]
      sc1 = sites_cart.select(regions_list == i)
      for l in range(k+1, len(unique_values)):
        j = unique_values[l]
        sc2 = sites_cart.select(regions_list == j)
        dist, i1, i2 = sc1.min_distance_between_any_pair_with_id(sc2)
        if dist and (best_dist is None or (dist < best_dist)):
          best_i = i
          best_j = j
          best_dist = dist
    if best_i is not None:
      print("Merging groups with distance of %.2f A" %(best_dist), file = log)

  if best_i is not None:
    sel = (regions_list == best_j)
    regions_list.set_selected(sel, best_i)
    update_regions_list(regions_list)

  return regions_list

def get_unique_values(regions_list):
  unique_values = []
  for x in regions_list:
    if not x in unique_values:
      unique_values.append(x)
  return unique_values

def swap_close_regions(sites_cart, regions_list, minimum_domain_length,
    close_distance = None):

  # Count number of residues in each pair that are close to the other
  # Split a group if some residues are close to other and not to self

  close_to_other_info  = get_close_to_other_list(sites_cart, regions_list,
     close_distance)
  close_to_other_list = close_to_other_info.close_to_other_list

  closer_to_other_swaps = get_closer_to_other(close_to_other_list,
      minimum_domain_length)
  found_something = False
  # Apply close swaps
  for s in closer_to_other_swaps:
    c = s.k_list[0]
    for k in range(c.start, c.end+1):
      regions_list[k] = s.j
      found_something = True

  update_regions_list(regions_list)

  if found_something:
    return regions_list

def get_close_to_other_list(sites_cart, regions_list, close_distance):

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
  for i in id_list:
    for j in id_list:
      if i==j: continue
      n_close = 0 # number in i close to j
      n_possible = 0

      for k in range(sites_dict[i].size()):
         index = index_dict[i][k]
         distances = (sites_dict[j] - col(sites_dict[i][k])).norms()
         local_n_close = (distances < close_distance).count(True)
         if local_n_close > 0:
           n_close += 1
         n_possible += 1
         distances_self = (sites_dict[i] - col(sites_dict[i][k])).norms()
         self_local_n_close = (distances_self < close_distance).count(True) - 1
         if local_n_close > self_local_n_close:
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
        n_possible = n_possible,
        i = i,
        j = j,))
  return group_args(
    group_args_type = 'close to other and close list',
      n_close_list = n_close_list,
      close_to_other_list = close_to_other_list)

def get_closer_to_other(close_to_other_list, minimum_domain_length):
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
      if close_dict[i][j] >= 0: # minimum_domain_length//2:
        pass
      else:
        del close_dict[i][j]
        if not close_dict[i]:
          del close_dict[i]
  all_k_list = []
  for i in close_dict.keys():
    for j in close_dict[i].keys():
      k_list = get_k_list(i,j,close_to_other_list)
      k_list = merge_k_list(k_list, minimum_domain_length)
      if not k_list:continue
      all_k_list.append(group_args(
        group_args_type = 'k list',
        i = i,
        j = j,
        k_list = k_list,
       ))
  return all_k_list


def merge_k_list(k_list, minimum_domain_length):
  n = len(k_list)
  for i in range(n):
    last_n = len(k_list)
    k_list = merge_k_list_once(k_list)
    if len(k_list) == last_n:
      break

  new_k_list = []
  for k1 in k_list:
    n1 = k1.end - k1.start + 1
    if n1 >= minimum_domain_length//3:
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

def replace_short_segments(regions_list, minimum_domain_length):
  id_list = get_unique_values(regions_list)
  new_regions_list = regions_list.deep_copy()
  for co_id in id_list:
    indices = (regions_list == co_id).iselection()
    indices_as_ranges = get_indices_as_ranges(indices)
    for r in indices_as_ranges:
      if r.end - r.start + 1 < minimum_domain_length:
        value = regions_list[r.start - 1] if r.start > 0 else \
           regions_list[min(regions_list.size() - 1, r.end + 1)]
        for i in range(r.start,r.end+1):
          new_regions_list[i] = value

  regions_list = new_regions_list
  update_regions_list(regions_list)
  return regions_list

def update_regions_list(regions_list):
  id_list = get_unique_values(regions_list)
  new_id_dict = {}
  i = 0
  for id_value in id_list:
    i += 1
    new_id_dict[id_value] = i
  new_id_list = list(new_id_dict.keys())

  for i in range(regions_list.size()):
    regions_list[i] = new_id_dict[regions_list[i]]

def get_indices_as_ranges(indices):
  indices = sorted(indices)
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
  avg_value = map_data.as_1d().min_max_mean().mean

  # Find max number of clusters in range of 0.5 to 1.0 * max
  n = 100
  max_clusters = None
  cutoff= None
  for t in range(int(min_cutoff*n),n+1):
    threshold = avg_value + t * (max_value-avg_value)/n
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

################################################################################
####################   Convenience functions          ##########################
################################################################################
def set_high_pae_for_missing(pae_matrix, pae_cutoff,
      residues_remaining):
   n,n = tuple(pae_matrix.shape)
   pae_1d = pae_matrix.flatten().tolist()
   skip_this_one = []
   for i in range(n):
     if i in residues_remaining:
        skip_this_one.append(False)
     else:
        skip_this_one.append(True)
   n_skipped = 0
   ii = -1
   for i in range(n):
     for j in range(n):
       ii += 1
       if skip_this_one[i] or skip_this_one[j]:
         pae_1d[ii] = pae_cutoff + 10
         n_skipped += 1
   import numpy

   matrix = numpy.empty((n, n))

   matrix.ravel()[:] = pae_1d
   return matrix


def get_residue_numbers_in_model(m_ca, remove_offset_of = None):
  residue_numbers = []
  for at in m_ca.get_hierarchy().atoms():
    resseq_int = at.parent().parent().resseq_as_int()
    if remove_offset_of is not None:
      resseq_int = resseq_int - remove_offset_of
    residue_numbers.append(resseq_int)
  return residue_numbers

def add_model(s1, s2):
  ''' add chains from s2 to existing s1'''
  s1 = s1.deep_copy()
  s1_ph = s1.get_hierarchy() # working hierarchy
  existing_chain_ids = s1_ph.chain_ids()
  for model_mm_2 in s2.get_hierarchy().models()[:1]:
    for chain in model_mm_2.chains():
      assert chain.id not in existing_chain_ids # duplicate chains in add_model
      new_chain = chain.detached_copy()
      for model_mm in s1_ph.models()[:1]:
        model_mm.append_chain(new_chain)

  s1.reset_after_changing_hierarchy()
  return s1
################################################################################
####################   END Convenience functions          ######################
################################################################################

if __name__ == "__main__":
  # run a simple version by default to demo usage
  args = sys.argv[1:]

  master_phil = iotbx.phil.parse(master_phil_str)
  params = master_phil.extract()
  master_phil.format(python_object=params).show(out=sys.stdout)
  p = params.process_predicted_model

  if len(args) < 2:
    print("libtbx.python process_predicted_model.py input.pdb output.pdb")
  else:
    input_file_name = args[0]
    output_file_name = args[1]
    if len(args) > 2:
       p.b_value_field_is = args[2]
    else:
       p.b_value_field_is = 'lddt'
    if len(args) > 3:
       p.domain_size = float(args[3])
    else:
       p.domain_size = 15
    from iotbx.data_manager import DataManager
    dm = DataManager()
    dm.set_overwrite(True)
    m = dm.get_model(input_file_name)

    p.remove_low_confidence_residues = True
    p.maximum_rmsd = 1.5
    p.split_model_by_compact_regions = True

    print("\nProcessing and splitting model into domains")
    model_info = process_predicted_model(m,  params)

    chainid_list = model_info.chainid_list
    print("Segments found: %s" %(" ".join(chainid_list)))

    mmm = model_info.model.as_map_model_manager()
    mmm.write_model(output_file_name)
    for chainid in chainid_list:
      selection_string = "chain %s" %(chainid)
      ph = model_info.model.get_hierarchy()
      asc1 = ph.atom_selection_cache()
      sel = asc1.selection(selection_string)
      m1 = model_info.model.select(sel)
      dm.write_model_file(m1, '%s_%s.pdb' %(output_file_name[:-4],chainid))
