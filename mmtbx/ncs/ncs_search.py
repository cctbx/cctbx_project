from __future__ import division
from scitbx.linalg import eigensystem
from scitbx.array_family import flex
from scitbx.math import superpose
from operator import itemgetter
import mmtbx.ncs.ncs_utils as nu
from libtbx.utils import Sorry
from scitbx import matrix
import math

class NCS_groups_container(object):

  def __init__(self):
    """
    A Container for ncs groups
    Note that the first copy is the master ncs

    Attributes:
      iselections (list of flex.size_t):selection array for the complete ASU
      residue_index_list (list): list of list of matching residues indices
      copies (list of lists):List of lists of the ncs copies chain IDs
      transforms (list of transform objects):
    """
    self.iselections = []
    self.residue_index_list = []
    self.copies = []
    self.transforms = []

class Transform(object):

  def __init__(self,
               rotation = None,
               translation = None,
               serial_num = None,
               coordinates_present = None,
               ncs_group_id = None,
               rmsd = 0):
    """
    Basic transformation properties

    Args:
      rotation : Rotation matrix object
      translation: Translation matrix object
      serial_num : (int) Transform serial number
      coordinates_present: equals 1 when coordinates are presents in PDB file
      ncs_group_id : (int) The group ID of the group containing this transform
      rmsd (float): RMS distance between ncs copies
    """
    self.r = rotation
    self.t = translation
    self.serial_num = serial_num
    self.coordinates_present = bool(coordinates_present)
    self.ncs_group_id = ncs_group_id
    self.rmsd = rmsd

class Score_record(object):

  def __init__(self,
               min_contig_length=10,score=0,
               gap_penalty=1,origin=(0,0)):
    """
    score object used when aligning sequences

    Attributes:
      score (int)
      consecutive_matches (int): num of consecutive matches in current segment
      match_count (int): number of matching residues
      min_contig_length (int): segments < min_contig_length rejected
      gap_penalty (int): penalty for gaps
      origin (tuple): (row,col) of the matrix cell we from the previous
        alignment step. Used to trace back optimal alignment
    """
    self.score = score
    self.consecutive_matches = 0
    self.match_count = 0
    self.min_contig_length = min_contig_length
    self.gap_penalty = gap_penalty
    self.origin = origin

  def reset_match_counting(self,seq_length):
    """
    Args:
      seq_length (int): length of shortest of the compared lists

    When there is a gap in the matching, reset the consecutive_matches and
    adjust the score if consecutive_matches < min_contig_length
    """
    min_matches = min(self.min_contig_length,seq_length)
    if self.consecutive_matches < min_matches:
      self.score -= self.consecutive_matches * self.gap_penalty
      self.match_count -= self.consecutive_matches
    self.consecutive_matches = 0

def find_ncs_in_hierarchy(ph,
                          min_contig_length=10,
                          min_fraction_domain=0.2,
                          use_minimal_master_ncs=True,
                          rmsd_eps=5.0):
  """
  Find NCS relation in hierarchy

  Args:
    ph (object): hierarchy
    min_contig_length (int): minimum length of matching chain segments
    min_fraction_domain (float): Threshold for similarity between chains
      similarity define as:
      (number of matching res) / (number of res in longer chain)
    use_minimal_master_ncs (bool): use maximal or minimal common chains
        in master ncs groups
    rmsd_eps (float): limit of rms difference chains when aligned together

  Return:
    groups_list (list of NCS_groups_container objects)
    group_dict (dict):
      keys: tuple of master chain IDs
      values: NCS_groups_container objects
  """
  # Get the list of matching chains
  chain_match_list = search_ncs_relations(
    ph=ph,
    min_contig_length=min_contig_length,
    min_fraction_domain=min_fraction_domain)
  #
  match_dict = clean_chain_matching(chain_match_list,ph,rmsd_eps)
  #
  if use_minimal_master_ncs:
    transform_to_group,match_dict = minimal_master_ncs_grouping(match_dict)
  else:
    transform_to_group,match_dict = minimal_ncs_operators_grouping(match_dict)
  #
  group_dict = build_group_dict(transform_to_group,match_dict)
  return group_dict

def minimal_ncs_operators_grouping(match_dict):
  """
  Look for NCS groups with the smallest number of chains in the master copy
  This is not the minimal number of NCS operations

  Args:
    match_dict(dict):
      key:(chains_id_1,chains_id_2)
      val:[selection_1,selection_2,res_list_1,res_list_2,rot,trans,rmsd]
      chain_ID (str), selection_1/2 (flex.size_t)
      res_list_1/2 (list): list of matching residues indices
      rot (matrix obj): rotation matrix
      tran (matrix obj): translation vector
      rmsd (float)

  Returns:
    group_dict (dict):
      keys: tuple of master chain IDs
      values: NCS_groups_container objects
    match_dict (dict): updated match_dict
  """
  # get a sorted list of match_dict key (by chain names)
  sorted_chain_groups_keys = sorted(match_dict)
  # collect all chain IDs
  chain_ids = {c_id for key in sorted_chain_groups_keys for c_id in key}
  chain_ids = sorted(chain_ids)
  # for each group in chain_groups find the transforms and group them
  transform_to_group = {}
  copy_to_transform = {}
  chains_in_groups = set()
  tr_sn = 0
  #
  for (master_id, copy_id) in sorted_chain_groups_keys:
    [sel_1,sel_2,res_1,res_2,r,t,rmsd] = match_dict[master_id,copy_id]
    # check if master is not a copy in another group
    if master_id in chains_in_groups: continue
    if copy_id in chains_in_groups: continue
    # check if the chains related by ncs operation
    tr_num, is_transpose = find_same_transform(r,t,transform_to_group)
    if is_transpose:
      # master and copy related by transpose of the rotation r
      if master_id in transform_to_group[tr_num][0]:
        # master_id is already a master for transform tr_num
        tr_num = None
        temp_master, temp_copy = master_id, copy_id
      else:
        temp_master, temp_copy = copy_id, master_id
        # flip the key in match_dict
        match_dict.pop((master_id, copy_id))
        # replace rotation and translation
        r,t = inverse_transform(r,t)
        match_dict[copy_id,master_id] = [sel_2,sel_1,res_2,res_1,r,t,rmsd]
    else:
      temp_master, temp_copy = master_id, copy_id
    if tr_num:
      # get chain IDs of all copies sharing this transform
      copies_ids = transform_to_group[tr_num][1]
      if not (temp_master in copies_ids):
        if copy_to_transform.has_key(temp_copy):
          # get the transforms that are leading to temp_copy
          tr_list = copy_to_transform.pop(temp_copy,None)
          for tr in tr_list:
            # since we increase the number of copies using the transform
            # we can remove the transform from its previous group
            transform_to_group.pop(tr,None)
          # update with new transform number
          copy_to_transform[temp_copy] = [tr_num]
        # remove transforms for which the new master was a copy,
        if copy_to_transform.has_key(temp_master):
          tr_list = copy_to_transform.pop(temp_master,None)
          for tr in tr_list:
            transform_to_group.pop(tr,None)
        # update dictionaries with additions to master and copies
        transform_to_group[tr_num][0].append(temp_master)
        transform_to_group[tr_num][1].append(temp_copy)
        chains_in_groups.update(set(transform_to_group[tr_num][1]))
    else:
      # create a new transform object and update groups
      tr_sn += 1
      transform_to_group[tr_sn] = [[temp_master],[temp_copy],(r,t)]
      # for each copy, collect all the transforms used to make it
      if copy_to_transform.has_key(temp_copy):
        copy_to_transform[temp_copy].append(tr_sn)
      else:
        copy_to_transform[temp_copy] = [tr_sn]

  # remove masters if appear in copies
  for k,v in transform_to_group.iteritems():
    masters = v[0]
    copies = v[1]
    common_ids = set(masters) & set(copies)
    while common_ids:
      ch_id = common_ids.pop()
      i = masters.index(ch_id)
      masters.pop(i)
      copies.pop(i)
      common_ids = set(masters) & set(copies)
    transform_to_group[k][0] = masters
    transform_to_group[k][1] = copies

  #clean up singles
  multiple_tr_num = [k for k,v in transform_to_group.iteritems() if
                    (len(v[0]) > 1)]
  singles_tr_num = [k for k,v in transform_to_group.iteritems() if
                    (len(v[0]) == 1)]
  for m_tr_num in multiple_tr_num:
    m_m = transform_to_group[m_tr_num][0]
    m_c = transform_to_group[m_tr_num][1]
    for s_tr_num in singles_tr_num:
      if transform_to_group.has_key(s_tr_num):
        s_m = transform_to_group[s_tr_num][0][0]
        s_c = transform_to_group[s_tr_num][1][0]
        # test if the single transform is part of a larger one
        t1 = (s_m in m_m) and (s_c in m_c)
        t2 = (s_c in m_m) and (s_m in m_c)
        t3 = (s_c in m_c) and (s_m in m_c)
        t4 = (s_c in m_m) and (s_m in m_m)
        t5 = (s_c in chains_in_groups) and (s_m in chains_in_groups)
        if t1 or t2 or t3 or t4 or t5:
          transform_to_group.pop(s_tr_num,None)

  # remove overlapping selection, keep largest groups
  master_size = [[len(v[0]),k] for k,v in transform_to_group.iteritems()]
  master_size.sort()
  chain_left_to_add = set(chain_ids)
  transform_to_use = set()
  chains_in_master = set()
  chains_in_copies = set()
  while bool(chain_left_to_add) and bool(master_size):
    [n,k] = master_size.pop()
    [masters,copies,(rot,trans)] = transform_to_group[k]
    # check that all chains in copies still need to be added
    test1 = len(set(copies) & chain_left_to_add) == len(copies)
    # check that copies are not in masters
    test2 = not bool(set(copies).intersection(chains_in_master))
    # check that masters are not in copies
    test3 = not bool(set(masters).intersection(chains_in_copies))
    if test1 and test2 and test3:
      transform_to_use.add(k)
      chains_in_master.update(set(masters))
      chains_in_copies.update(set(copies))
      chain_left_to_add -= set(copies)
      chain_left_to_add -= set(masters)

  # delete all but transform_to_use
  keys = transform_to_group.keys()
  for key in keys:
    if not key in transform_to_use:
      transform_to_group.pop(key)

  return transform_to_group,match_dict

def minimal_master_ncs_grouping(match_dict):
  """
  Look for NCS groups with the smallest number of chains in the master copy
  This is not the minimal number of NCS operations

  Args:
    match_dict(dict):
      key:(chains_id_1,chains_id_2)
      val:[selection_1,selection_2,res_list_1,res_list_2,rot,trans,rmsd]
      chain_ID (str), selection_1/2 (flex.size_t)
      res_list_1/2 (list): list of matching residues indices
      rot (matrix obj): rotation matrix
      tran (matrix obj): translation vector
      rmsd (float)

  Returns:
    transform_to_group (dict):
      keys: tuple of master chain IDs
      values: NCS_groups_container objects
    match_dict (dict): updated match_dict
  """
  # get a sorted list of match_dict key (by chain names)
  sorted_chain_groups_keys = sorted(match_dict)
  # collect all chain IDs
  chain_ids = {c_id for key in sorted_chain_groups_keys for c_id in key}
  chain_ids = sorted(chain_ids)
  # for each group in chain_groups find the transforms and group them
  transform_to_group = {}
  chains_in_groups = set()
  tr_sn = 0
  #
  for (master_id, copy_id) in sorted_chain_groups_keys:
    [sel_1,sel_2,res_1,res_2,r,t,rmsd] = match_dict[master_id,copy_id]
    # check if master is not a copy in another group
    if (master_id in chains_in_groups) or (copy_id in chains_in_groups):
      match_dict.pop((master_id, copy_id),None)
      continue
    # check if the chains related by ncs operation
    tr_num, is_transpose = find_same_transform(r,t,transform_to_group)
    if is_transpose:
      # master and copy related by transpose of the rotation r
      if master_id in transform_to_group[tr_num][0]:
        # master_id is already a master for transform tr_num
        tr_num = None
        temp_master, temp_copy = master_id, copy_id
      else:
        temp_master, temp_copy = copy_id, master_id
        # flip the key in match_dict
        match_dict.pop((master_id, copy_id),None)
        # replace rotation and translation
        r,t = inverse_transform(r,t)
        match_dict[copy_id,master_id] = [sel_2,sel_1,res_2,res_1,r,t,rmsd]
    else:
      temp_master, temp_copy = master_id, copy_id
    #
    chains_in_groups.add(temp_copy)
    if tr_num:
      # update dictionaries with additions to master and copies
      transform_to_group[tr_num][0].append(temp_master)
      transform_to_group[tr_num][1].append(temp_copy)
    else:
      # create a new transform object and update groups
      tr_sn += 1
      transform_to_group[tr_sn] = [[temp_master],[temp_copy],(r,t)]

  # If a chain is a master or a copy for more than one transform,
  # choose only the one with the smaller number of chains in master
  master_size = [[len(v[0]),k] for k,v in transform_to_group.iteritems()]
  # sort list by the number of chains in master for each transform
  master_size.sort()
  chain_left_to_add = set(chain_ids)
  transform_to_use = set()
  # Collect the transforms, starting with the smallest master
  while bool(chain_left_to_add) and bool(master_size):
    [n,k] = master_size.pop(0)
    tr = transform_to_group[k]
    if bool(set(tr[1]) & chain_left_to_add):
      transform_to_use.add(k)
      chain_left_to_add -= set(tr[1])
      chain_left_to_add -= set(tr[0])

  # delete all transforms that are not being used
  tr_sns = transform_to_group.keys()
  for tr_sn in tr_sns:
    if not tr_sn in transform_to_use:
      transform_to_group.pop(tr_sn)

  return transform_to_group,match_dict

def build_group_dict(transform_to_group,match_dict):
  """
  find all transforms with common masters and organize the chains
  in the same order as the master and build groups transform dictionary

  Args:
    transform_to_group (dict):
      key: temporary transform number
      values:[masters],[copies],(rotation,translation)]
    match_dict(dict):
      key:(chains_id_1,chains_id_2)
      val:[selection_1,selection_2,rot,trans,rmsd]
          chain_ID (str), selection_1/2 (flex.size_t)
          residue_index_list_1/2 (list): matching residues indices
          rot (matrix obj): rotation matrix
          tran (matrix obj): translation vector
          rmsd (float)

  Returns:
    group_dict (dict):
      keys: tuple of master chain IDs
      values: NCS_groups_container objects
  """
  group_dict = {}
  group_id = 0
  tr_sn = 0
  for _,v in transform_to_group.iteritems():
    [masters,copies,(r,t)] = v
    # sort masters chain IDs
    sorted_masters = sorted(masters)
    # reorder the copies order to match the new master order
    sort_map = [masters.index(x) for x in sorted_masters]
    copies = [copies[i] for i in sort_map]
    # use the chain IDs of the master as a key for group_dict
    key = tuple(sorted_masters)
    m_sel_list,c_sel_list = get_iselection(sorted_masters,copies,match_dict)
    if group_dict.has_key(key):
      zipped_param = zip(m_sel_list,c_sel_list,group_dict[key].iselections[0])
      for i,rec in enumerate(zipped_param):
        (m_sel,c_sel,current_master) = rec
        # for each chain, check if the master have the same selection
        current_master_set = set(current_master - min(current_master))
        m_sel_set = set(m_sel - min(m_sel))
        if current_master_set != m_sel_set:
          updated_iselection,c_sel = update_selections(
            group_dict[key].iselections,m_sel,c_sel,i)
          c_sel_list[i] = c_sel
          group_dict[key].iselections = updated_iselection
      #
      tr_sn += 1
      [_,_,res_l_m,res_l_c,r,t,rmsd] = match_dict[sorted_masters[0],copies[0]]
      tr = Transform(
        rotation=r,
        translation=t,
        serial_num=tr_sn,
        coordinates_present=True,
        ncs_group_id=group_id,
        rmsd=rmsd)
      group_dict[key].iselections.append(c_sel_list)
      group_dict[key].residue_index_list.append(res_l_c)
      group_dict[key].copies.append(copies)
      group_dict[key].transforms.append(tr)
    else:
      tr_sn += 1
      group_id += 1
      [_,_,res_l_m,res_l_c,r,t,rmsd] = match_dict[sorted_masters[0],copies[0]]
      # add master as first copy (identity transform)
      new_ncs_group = NCS_groups_container()
      #
      tr = Transform(
        rotation=matrix.sqr([1,0,0,0,1,0,0,0,1]),
        translation=matrix.col([0,0,0]),
        serial_num=tr_sn,
        coordinates_present=True,
        ncs_group_id=group_id,
        rmsd=0)
      new_ncs_group.iselections.append(m_sel_list)
      new_ncs_group.residue_index_list.append(res_l_m)
      new_ncs_group.copies.append(sorted_masters)
      new_ncs_group.transforms.append(tr)
      # add the copy
      tr_sn += 1
      tr = Transform(
        rotation=r,
        translation=t,
        serial_num=tr_sn,
        coordinates_present=True,
        ncs_group_id=group_id,
        rmsd=rmsd)
      new_ncs_group.iselections.append(c_sel_list)
      new_ncs_group.residue_index_list.append(res_l_c)
      new_ncs_group.copies.append(copies)
      new_ncs_group.transforms.append(tr)
      #
      group_dict[key] = new_ncs_group
  return group_dict

def update_selections(iselections,m_sel,c_sel,i):
  """
  Update the total master and copies selection with the smallest common
  selection, among all ncs copies

  Args:
    iselections (list): list of flex.size_t iselection for all NCS copies
    m_sel,c_sel (flex.size_t): iselection for master and a single copy
    i (int): chain number in iselections

  Returns:
    iselections : updated iselections
    c_sel : updated c_sel
  """
  current_master = iselections[0][i]
  m_sel_set = set(m_sel - min(m_sel))
  master_offset = min(current_master)
  current_master_set = set(current_master - master_offset)
  #
  remove_from_old_ncs_copies = current_master_set - m_sel_set
  remove_from_new_copy = m_sel_set - current_master_set
  #
  if remove_from_old_ncs_copies:
    new_iselections = []
    for isel_list in iselections:
      isel = isel_list[i]
      offset = min(isel)
      isel_set = set(isel - offset)
      isel = list(isel_set - remove_from_old_ncs_copies)
      isel = flex.size_t(isel) + offset
      isel_list[i] = isel
      new_iselections.append(isel_list)
  else:
    new_iselections = iselections
  if remove_from_new_copy:
    offset = min(c_sel)
    c_sel_set = set(c_sel - offset)
    c_sel = list(c_sel_set - remove_from_new_copy)
    c_sel = flex.size_t(c_sel) + offset
  return new_iselections,c_sel

def get_iselection(sorted_masters,copies,match_dict):
  """
  Combine iselection of all chains in NCS master and NCS copy to a single
  iselelction

  Args:
    sorted_masters (list): list of master chain IDs
    copies (list): list of copies chain IDs
    match_dict(dict):
      key:(chains_id_1,chains_id_2)
      val:[selection_1,selection_2,rot,trans,rmsd]
      chain_ID (str), selection_1/2 (flex.size_t)
      residue_index_list_1/2 (list): matching residues indices
      rot (matrix obj): rotation matrix
      tran (matrix obj): translation vector
      rmsd (float)

  Returns:
    m_sel,c_sel (list of flex.size_t): list of selected atoms per chain
  """
  m_sel = []
  c_sel = []
  for key in zip(sorted_masters,copies):
    [sel_1,sel_2,_,_,_,_,_] = match_dict[key]
    m_sel.append(sel_1)
    c_sel.append(sel_2)
  # Note: consider sorting the selections
  return m_sel,c_sel

def clean_chain_matching(chain_match_list,ph,rmsd_eps=5):
  """
  Remove all bad matches from chain_match_list

  Args:
    ph (object): hierarchy
    chain_match_list (list): list of
      [chain_ID_1, chain_ID_2, sel_1, sel_1,res_m/res_c similarity]
      chain_ID (str), sel_1/2 (flex.bool)
      res_m/res_c (lists): indices of the aligned components
      similarity (float): similarity between chains
    rmsd_eps (float): limit of rms difference chains

  Returns:
    match_dict(dict): key:(chains_id_1,chains_id_2)
                      val:[selection_1,selection_2,
                           res_list_1,res_list_2,rot,trans,rmsd]
  """
  # remove all non-matching pairs, where similarity == 0
  match_list = [x for x in chain_match_list if x[4] > 0]
  # keep only best (or 95% of best) matches
  best_matches = {}
  # Get rmsd
  match_dict = {}
  for match in match_list:
    [ch_a_id, ch_b_id, sel_a, sel_b,res_list_1,res_list_2,similarity] = match
    best_matches,match_dict = update_best_matches_dict(
      best_matches,match_dict,ch_a_id,ch_b_id,similarity)
    other_sites = ph.select(sel_a).atoms().extract_xyz()
    ref_sites = ph.select(sel_b).atoms().extract_xyz()
    lsq_fit_obj = superpose.least_squares_fit(
      reference_sites = ref_sites,
      other_sites     = other_sites)
    r = lsq_fit_obj.r
    t = lsq_fit_obj.t
    rmsd = ref_sites.rms_difference(lsq_fit_obj.other_sites_best_fit())
    rmsd = round(rmsd,4)
    match_dict[ch_a_id, ch_b_id] = [sel_a,sel_b,res_list_1,res_list_2,r,t,rmsd]
  # Clean using rmsd limit
  match_dict = {k:v for (k,v) in match_dict.iteritems() if v[-1] <= rmsd_eps}
  return match_dict

def update_best_matches_dict(best_matches,match_dict,ch_a_id,ch_b_id,similarity):
  """
  Updates the best_matches dictionaries best_matches,match_dict to keep only
  matches that are at least 95% of best match

  Args:
    best_matches (dict):
      key: chain_ID (str)
      Val: [[similarity,(chain_ID,chain_ID2)],[similarity,(chain_ID3,chain_ID)]]
    match_dict (dict):
      key: (chain_id_a,chain_id_b)
      val: [selection_1,selection_2,res_list_1,res_list_2,rot,trans,rmsd]
    ch_a_id,ch_b_id (str): chain IDs
    similarity (float): similarity between chains

  Returns:
    updated best_matches,match_dict
  """
  records_to_remove = set()
  for ch_id in [ch_a_id,ch_b_id]:
    if best_matches.has_key(ch_id):
      temp_rec = []
      # records with largest similarity is last
      max_sim = best_matches[ch_id][-1][0]
      if similarity > max_sim:
        for s,(a,b) in best_matches[ch_id]:
          if s/similarity >= 0.95:
            temp_rec.append([s,(a,b)])
          else:
            records_to_remove.add((a,b))
        temp_rec.append([similarity,(ch_a_id,ch_b_id)])
      elif similarity/max_sim >= 0.95:
        temp_rec.insert(0,[similarity,(ch_a_id,ch_b_id)])
      else:
        temp_rec = best_matches[ch_id]
      best_matches[ch_id] = temp_rec
    else:
      best_matches[ch_id] = [[similarity,(ch_a_id,ch_b_id)]]
  # clean match_dict
  if records_to_remove:
    for key in records_to_remove:
      match_dict.pop(key,None)
  return best_matches,match_dict


def find_same_transform(r,t,transforms,eps=0.1,angle_eps=5):
  """
  Check if the rotation r and translation t exist in the transform dictionary.

  Comparing rotations and the result of applying rotation and translation on
  a test vector

  Args:
    r (matrix.sqr): rotation
    t (matrix.col): translation
    transforms (dict): dictionary of all transforms
    angle_eps (float): allowed difference in similar rotations
    eps (float): allowed difference in the average distance between

  Returns:
    tr_num: (str) transform serial number
    is_transpose: (bool) True if the matching transform is transpose
  """
  is_transpose = False
  tr_num = None
  for k,v in transforms.iteritems():
    if hasattr(v,'r'):
      rr = v.r
      tt = v.t
    else:
      (rr,tt) = v[2]
    is_the_same, is_transpose = is_same_transform(
      r, t, rr, tt, eps=eps, angle_eps=angle_eps)
    if is_the_same:
      return k, is_transpose
  return tr_num, is_transpose

def search_ncs_relations(ph,
                         min_contig_length=10,
                         min_fraction_domain=0.2):
  """
  Search for NCS relations between chains or parts of chains, in a protein
  hierarchy

  Args:
    ph (object): hierarchy
    min_contig_length (int): segments < min_contig_length rejected
    min_fraction_domain (float): Threshold for similarity between chains.
      similarity define as:
      (number of matching res) / (number of res in longer chain)

  Returns:
    chain_match_list (list): list of
      [chain_ID_1,chain_ID_2,sel_1,sel_2,res_sel_m, res_sel_c,similarity]
      chain_ID (str), sel_1 (flex.bool), sel_1 (flex.bool),
      res_sel_m/c (lists): indices of the aligned components
      similarity (float): similarity between chains
    We use selection_1 and selection_2 because in some PDB records a residue
    might contain side chain in one chain but not in another
  """
  # collect all chain IDs
  cache = ph.atom_selection_cache()
  chain_match_list = []
  model  = ph.models()[0]
  chain_ids = list({x.id for x in model.chains()})
  chain_ids = sorted(chain_ids)
  n_chains = len(chain_ids)
  ph_chains = model.chains()
  # loop over all chains
  for i in xrange(n_chains-1):
    master_ch_id = ph_chains[i].id
    m_selection = cache.selection('chain ' + master_ch_id)
    m_ph = ph.select(m_selection)
    m_chain = ph_chains[i]
    for j in xrange(i+1,n_chains):
      copy_ch_id = ph_chains[j].id
      c_selection = cache.selection('chain ' + copy_ch_id)
      c_ph = ph.select(c_selection)
      c_chain = ph_chains[j]
      # res_a/b are lists of (residue name, residue number)
      chain_id_m, res_m = get_residue_sequence(m_ph)
      chain_id_c, res_c = get_residue_sequence(c_ph)
      # get residue lists
      seq_m = [x for (x,y) in res_m]
      seq_c = [x for (x,y) in res_c]
      res_sel_m, res_sel_c, similarity = res_alignment(
        seq_a=seq_m,seq_b=seq_c,
        min_contig_length=min_contig_length,
        min_fraction_domain=min_fraction_domain)
      sel_m, sel_c,res_sel_m,res_sel_c = get_matching_atoms(
        m_ph,c_ph,res_sel_m,res_sel_c,m_chain,c_chain)
      if res_sel_m:
        # add only non empty matches
        sel_m = sel_m.iselection()
        sel_c = sel_c.iselection()
        # adjust the selection in the chain to the complete hierarchy
        sel_m += min(m_selection.iselection())
        sel_c += min(c_selection.iselection())
        rec = [master_ch_id,copy_ch_id,sel_m,sel_c,res_sel_m,res_sel_c,
               similarity]
        chain_match_list.append(rec)
  return chain_match_list

def res_alignment(seq_a, seq_b,
                  min_contig_length=10,
                  min_fraction_domain=0.2):
  """
  Align two chains hierarchies.

  Penalize misalignment, gaps, when contiguous section is less
  than min_contig_length and when the
  Do not give any points for alignment. (score only change when "bad" things
  happens)
  Stop if the alignment is less than the "similarity_threshold"
  After alignment verify that aligned parts are >= similarity_threshold

  Args:
    seq_a, seq_b (lists of str): list of characters to compare
    min_fraction_domain (float): min percent of similarity between hierarchies
      similarity define as:
      (number of matching res) / (number of res in longer chain)
    min_contig_length (int): domain < min_contig_length rejected

  Returns:
    aligned_sel_a (list): the indices of the aligned components of seq_a
    aligned_sel_b (list): the indices of the aligned components of seq_b
    similarity (float): actual similarity between hierarchies
  """
  a = len(seq_a)
  b = len(seq_b)
  if (a == 0) or (b == 0): return [],[],0
  gap_penalty = 1
  misalign_penalty = 2 * (a + b)
  # limit the number of mis-alignments
  max_mis_align = int((1 - min_fraction_domain) * (a + b))
  # Starting score according to the required similarity
  score = max_mis_align
  # build score matrix
  R = build_score_matrix(
    row=a,col=b,max_score=score,
    gap_penalty=gap_penalty,min_contig_length=min_contig_length)
  # populate score matrix
  for j in xrange(1,b + 1):
    for i in xrange(1,a + 1):
      not_aligned = (seq_a[i-1].lower() != seq_b[j-1].lower())
      s1 = R[i-1,j-1].score - misalign_penalty * not_aligned
      s2 = R[i-1,j].score - gap_penalty
      s3 = R[i,j-1].score - gap_penalty
      s = max(s1,s2,s3)
      R[i,j] = Score_record(min_contig_length=min_contig_length,score=s)
      if s1 == s:
        R[i,j].consecutive_matches = R[i-1,j-1].consecutive_matches + 1
        R[i,j].origin = (i-1,j-1)
        R[i,j].match_count = R[i-1,j-1].match_count + 1
      else:
        if s2 == s:
          (i_temp,j_temp) = (i-1,j)
        else:
          (i_temp,j_temp) = (i,j-1)
        R[i,j].origin = (i_temp,j_temp)
        R[i,j].match_count = R[(i_temp,j_temp)].match_count
        R[i,j].consecutive_matches = R[(i_temp,j_temp)].consecutive_matches
        R[i,j].reset_match_counting(min(a,b))
  R[i,j].reset_match_counting(min(a,b))
  #
  aligned_sel_a, aligned_sel_b, similarity = get_matching_res_indices(
    R=R,row=a,col=b,min_fraction_domain=min_fraction_domain)
  return aligned_sel_a, aligned_sel_b, similarity

def get_matching_res_indices(R,row,col,min_fraction_domain):
  """
  Trace back best solution and collect the indices of matching pairs

  Args:
    R (list of lists of int): the score matrix
    row (int): number of rows
    col (int): number of columns
    min_fraction_domain (float): min percent of similarity between hierarchies
      similarity define as:
      (number of matching res) / (number of res in longer chain)

  Returns:
    sel_a (list): matching indices sequence a
    sel_b (list): matching indices sequence b
    similarity (float): actual similarity between hierarchies
  """
  # index of best score from last row
  j_max = col
  i_max = row
  sel_a = []
  sel_b = []
  # test alignment quality
  if not R.has_key((i_max,j_max)):
    # alignment process was halted, return empty arrays
    return flex.size_t(sel_a), flex.size_t(sel_b),0
  #
  similarity = R[i_max,j_max].match_count / max(i_max,j_max)
  if similarity < min_fraction_domain:
    # chains are to different, return empty arrays
    return flex.size_t(sel_a), flex.size_t(sel_b), 0
  #
  while (i_max > 0) and (j_max > 0):
    if R[i_max,j_max].origin == (i_max - 1, j_max - 1):
      sel_a.append(i_max - 1)
      sel_b.append(j_max - 1)
      i_max -= 1
      j_max -= 1
    elif R[i_max,j_max].origin == (i_max - 1, j_max):
      i_max -= 1
    elif R[i_max,j_max].origin ==(i_max, j_max - 1):
      j_max -= 1
    else:
      raise Sorry('Sequence alignment error')
  sel_a.reverse()
  sel_b.reverse()
  assert len(sel_a) == len(sel_b)
  return sel_a, sel_b, similarity


def get_matching_atoms(ph_a,ph_b,res_num_a,res_num_b,chain_a,chain_b):
  """
  Get selection of matching chains, match residues atoms
  We keep only residues with continuous matching atoms

  Args:
    ph_a (hierarchy): first chain
    ph_b (hierarchy): second chain
    res_num_a/b (list of int): indices of matching residues numbers, ph_a/b
    chain_a/b (chain hierarchy): chains hierarchies of chains a and b

  Returns:
    sel_a, sel_b (flex.bool): matching atoms selection
    res_num_a/b (list of int): updated res_num_a/b
  """
  atom_cache_a = ph_a.atom_selection_cache().selection
  atom_cache_b = ph_b.atom_selection_cache().selection
  l_a = ph_a.atoms().size()
  l_b = ph_b.atoms().size()
  sel_a = flex.bool([False,] * l_a)
  sel_b = flex.bool([False,] * l_b)
  # collect the residue IDs according to the res_num
  res_ids_a = [x.resid() for x in chain_a.residues()]
  res_ids_b = [x.resid() for x in chain_b.residues()]
  #
  res_num_a_updated = []
  res_num_b_updated = []
  for (resid_a,resid_b,na,nb) in zip(res_ids_a,res_ids_b,res_num_a,res_num_b):
    sa = atom_cache_a('resid ' + resid_a)
    sb = atom_cache_b('resid ' + resid_b)
    res_a = ph_a.select(sa)
    res_b = ph_b.select(sb)
    # select only atoms that exist in both residues
    atoms_names_a = list(res_a.atoms().extract_name())
    atoms_names_b = list(res_b.atoms().extract_name())
    atoms_a,atoms_b,similarity = res_alignment(
      seq_a=atoms_names_a, seq_b=atoms_names_b,
      min_contig_length=100)
    sa = flex.size_t(atoms_a)+min(sa.iselection())
    sb = flex.size_t(atoms_b)+min(sb.iselection())
    # keep only residues with continuous matching atoms
    if sa.size() != 0:
      res_num_a_updated.append(na)
      res_num_b_updated.append(nb)
    sa = flex.bool(l_a,sa)
    sb = flex.bool(l_b,sb)
    sel_a = sel_a | sa
    sel_b = sel_b | sb
  return sel_a,sel_b,res_num_a_updated,res_num_b_updated

def get_residue_sequence(chain_ph):
  """
  Get a list of pairs (residue name, residue number) and the chain ID
  Exclude water

  Args:
    chain_ph (iotbx_pdb_hierarchy_ext): hierarchy of a single chain

  Returns:
    chain_id (str)
    res_seq (list)
  """
  res_list = [x.resname for x in chain_ph.atom_groups()]
  resid_list = [x.resseq for x in chain_ph.residue_groups()]
  # Exclude water
  res_seq = [(x,y) for (x,y) in zip(res_list, resid_list) if x.lower() != 'hoh']
  #
  atoms = chain_ph.atoms()
  chain_id = atoms[0].chain().id
  return chain_id, res_seq

def build_score_matrix(row, col,
                       max_score=1000,
                       gap_penalty=1,
                       min_contig_length=10):
  """
  Build a matrix in a dictionary form
  The matrix values are initialized according the the gap penalties.

  We only initializing the zero row and column, the rest of the items are
  created while doing the alignment

  Args:
    row (int): number of rows
    col (int): number of columns
    max_score (int): the score assigned for perfect alignment
    gap_penalty (int): penalty for gaps
    min_contig_length (int): segments < min_contig_length rejected

  Returns:
    R (dict): the score matrix in a form of a dictionary of Score_record objects
  """
  R = {}
  # populate zero cases
  for i in xrange(row + 1):
    score = max_score - (i * gap_penalty)
    R[i,0] = Score_record(
      min_contig_length=min_contig_length,score=score,origin=(i-1,0))
  for i in xrange(col + 1):
    score = max_score - (i * gap_penalty)
    R[0,i] = Score_record(
      min_contig_length=min_contig_length,score=score,origin=(0,i-1))
  return R

def inverse_transform(r,t):
  """ inverse rotation and translation """
  r = r.transpose()
  t = - r*t
  return r,t

def angle_between_rotations(v1,v2):
  """ get angle between two vectors"""
  cos_angle = v1.dot(v2)
  result = math.acos(min(1,cos_angle))
  result *= 180/math.pi
  return result

def get_rotation_vec(r):
  """ get the eigen vector associated with the eigen value 1"""
  eigen = eigensystem.real_symmetric(r.as_sym_mat3())
  eigenvectors = eigen.vectors()
  eigenvalues = eigen.values()
  i = list(eigenvalues.round(4)).index(1)
  return eigenvectors[i:(i+3)]

def is_same_transform(r1,t1,r2,t2,eps=0.1,angle_eps=5):
  """
  Args:
    r1, r2: Rotation matrices
    t1, t2: Translation vectors
    eps, angle_eps: allowed numerical difference

  Returns:
    (bool,bool) (is_the_same, is_transpose)
  """
  if (not r1.is_zero()) and (not r2.is_zero()):
    assert r1.is_r3_rotation_matrix(rms_tolerance=0.001)
    assert r2.is_r3_rotation_matrix(rms_tolerance=0.001)
    # test vector
    xyz = flex.vec3_double([(11,103,523),(-500.0,2.0,10.0),(0.0,523.0,-103.0)])
    a_ref = (r1.elems * xyz + t1).as_double()
    rt, tt = inverse_transform(r1,t1)
    a_ref_transpose = (rt.elems * xyz + tt).as_double()
    v1 = get_rotation_vec(r1)
    v2 = get_rotation_vec(r2)
    a = (r2.elems * xyz + t2).as_double()
    d = (a_ref-a)
    d = (d.dot(d))**.5/a.size()
    dt = (a_ref_transpose-a)
    dt = (dt.dot(dt))**.5/a.size()
    ang = angle_between_rotations(v1,v2)
    d_ang = min(ang, (180 - ang))
    if (d_ang < angle_eps) and (d < eps):
      return True, False
    elif (d_ang < angle_eps) and (dt < eps):
      return True, True
    else:
      return False, False
  else:
    return False, False
