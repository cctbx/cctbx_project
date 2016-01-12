from __future__ import division
from scitbx.linalg import eigensystem
from scitbx.array_family import flex
from scitbx.math import superpose
from libtbx.utils import Sorry
from scitbx import matrix
import math
import sys

__author__ = 'Youval'

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

  def __init__(self,score=-10000,origin=(0,0)):
    """
    score object used when aligning sequences

    Attributes:
      score (int)
      consecutive_matches (int): num of consecutive matches in current segment
      match_count (int): number of matching residues
      origin (tuple): (row,col) of the matrix cell we from the previous
        alignment step. Used to trace back optimal alignment
      no_altloc (list of bool): False when residue has alternate location
    """
    self.score = score
    self.consecutive_matches = 0
    self.match_count = 0
    self.gap_penalty = 1
    self.origin = origin

class Chains_info(object):
  """ Container for hierarchy analysis """
  def __init__(self):
    self.res_names = []
    self.resid = []
    self.atom_names = []
    self.atom_selection = []
    self.chains_atom_number = 0
    self.no_altloc = []
    self.center_of_coordinates = None

  def __str__(self):
    from StringIO import StringIO
    assert 0
    res = StringIO()
    print >> res, "res_names:", self.res_names
    print >> res, "self.resid", self.resid
    print >> res, "self.atom_names", self.atom_names
    print >> res, "self.atom_selection", self.atom_selection
    print >> res, "self.chains_atom_number", self.chains_atom_number
    print >> res, "self.no_altloc", self.no_altloc
    print >> res, "self.center_of_coordinates", self.center_of_coordinates
    return res.getvalue()


def find_ncs_in_hierarchy(ph,
                          min_contig_length=10,
                          min_percent=0.95,
                          use_minimal_master_ncs=True,
                          max_rmsd=5.0,
                          write=False,
                          log=None,
                          check_atom_order=False,
                          allow_different_size_res=True,
                          exclude_misaligned_residues=False,
                          similarity_threshold=0.95,
                          match_radius=4.0,
                          ignore_chains=None):
  """
  Find NCS relation in hierarchy

  Args:
    ph (object): hierarchy
    min_contig_length (int): minimum length of matching chain segments
    min_percent (float): Threshold for similarity between chains
      similarity define as:
      (number of matching res) / (number of res in longer chain)
    use_minimal_master_ncs (bool): use maximal or minimal common chains
        in master ncs groups
    max_rmsd (float): limit of rms difference chains when aligned together
    write (bool): when true, write ncs search messages to log
    check_atom_order (bool): check atom order in matching residues.
        When False, matching residues with different number of atoms will be
        excluded from matching set
    exclude_misaligned_residues (bool): check and exclude individual residues
      alignment quality
    match_radius (float): max allow distance difference between pairs of matching
      atoms of two residues
    allow_different_size_res (bool): keep matching residue with different
      number of atoms
    similarity_threshold (float): min similarity between matching chains
    ignore_chains (set of str): set of chain IDs to exclude

  Return:
    groups_list (list of NCS_groups_container objects)
    group_dict (dict):
      keys: tuple of master chain IDs
      values: NCS_groups_container objects
  """
  if not log: log = sys.stdout
  chains_info = get_chains_info(ph,ignore_chains=ignore_chains)
  # Get the list of matching chains
  chain_match_list = search_ncs_relations(
    chains_info=chains_info,
    min_contig_length=min_contig_length,
    min_percent=min_percent,
    write=write,
    log=log,
    check_atom_order=check_atom_order,
    use_minimal_master_ncs=use_minimal_master_ncs,
    allow_different_size_res=allow_different_size_res,
    ignore_chains=ignore_chains)
  #
  match_dict = clean_chain_matching(
    chain_match_list=chain_match_list,
    ph=ph,
    max_rmsd=max_rmsd,
    exclude_misaligned_residues=exclude_misaligned_residues,
    match_radius=match_radius,
    similarity_threshold=similarity_threshold)
  # print "match_dict.keys()", match_dict.keys()
  # print "match_dict"
  # for k, v in match_dict.iteritems():
  #   print "  ", k, list(v[0]), list(v[1])
  # assert 0
  #

  # new, the basic way of processing, by Oleg.
  return ncs_grouping_and_group_dict(match_dict, ph)


  # this is not used anymore
  if use_minimal_master_ncs:
    transform_to_group,match_dict = minimal_master_ncs_grouping(match_dict, ph)
  else:
    transform_to_group,match_dict = minimal_ncs_operators_grouping(match_dict)
  #
  # print "match_dict"
  # for k, v in match_dict.iteritems():
  #   print "  ", k, list(v[0]), list(v[1])


  group_dict = build_group_dict(transform_to_group,match_dict,chains_info)
  # print "group_dict"
  # for k, v in group_dict.iteritems():
  #   print "  ", k, v
  # return group_dict

def minimal_ncs_operators_grouping(match_dict):
  """
  Look for NCS groups with the smallest number of chains in the master copy
  This is not the minimal number of NCS operations

  Args:
    match_dict(dict):
      key:(chains_id_1,chains_id_2)
      val:[select_1,select_2,res_list_1,res_list_2,rot,trans,rmsd]
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
          tr_list = copy_to_transform.pop(temp_copy,[])
          for tr in tr_list:
            # since we increase the number of copies using the transform
            # we can remove the transform from its previous group
            transform_to_group.pop(tr,None)
          # update with new transform number
          copy_to_transform[temp_copy] = [tr_num]
        # remove transforms for which the new master was a copy,
        if copy_to_transform.has_key(temp_master):
          tr_list = copy_to_transform.pop(temp_master,[])
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
  # Clean solution
  transform_to_group = remove_masters_if_appear_in_copies(transform_to_group)
  transform_to_group = clean_singles(transform_to_group,chains_in_groups)
  transform_to_group = remove_overlapping_selection(transform_to_group,chain_ids)
  return transform_to_group,match_dict

def get_rmsds2(master_xyz, copy_xyz, cur_ttg):
  """
  This function is for debugging purposes and not called.
  """
  xyz = cur_ttg[2][0].elems * master_xyz + cur_ttg[2][1]
  # rmsd1 = 0
  # if copy_xyz.size() == xyz.size():
  rmsd1 = copy_xyz.rms_difference(xyz)
  xyz = cur_ttg[2][0].elems * master_xyz + cur_ttg[2][1]
  # rmsd2 = 0
  # if copy_xyz.size() == xyz.size():
  rmsd2 = copy_xyz.rms_difference(xyz)
  # print "rmsds:", rmsd1, rmsd2
  return rmsd1, rmsd2

def get_rmsds(hierarchy, cache, cur_ttg, master, copy):
  """
  This function is for debugging purposes and not called.
  Similar check will be performed later in execution and in case of
  wrong grouping will raise Sorry: bad phil records.
  """
  str_sel_m = "chain "+" or chain ".join(cur_ttg[0]+[master])
  str_sel_c = "chain "+" or chain ".join(cur_ttg[1]+[copy])
  sel1 = cache.selection("chain "+" or chain ".join(cur_ttg[0]+[master]))
  sel2 = cache.selection("chain "+" or chain ".join(cur_ttg[1]+[copy]))
  # print "sel1, sel2", str_sel_m, "|", str_sel_c
  master_xyz = hierarchy.select(sel1).atoms().extract_xyz()
  copy_xyz = hierarchy.select(sel2).atoms().extract_xyz()
  xyz = cur_ttg[2][0].elems * master_xyz + cur_ttg[2][1]
  rmsd1 = 0
  if copy_xyz.size() == xyz.size():
    rmsd1 = copy_xyz.rms_difference(xyz)

  str_sel_m = "chain "+" or chain ".join(cur_ttg[0]+[copy])
  str_sel_c = "chain "+" or chain ".join(cur_ttg[1]+[master])
  # print "sel1, sel2", str_sel_m, "|", str_sel_c
  sel1 = cache.selection("chain "+" or chain ".join(cur_ttg[0]+[copy]))
  sel2 = cache.selection("chain "+" or chain ".join(cur_ttg[1]+[master]))
  # print "sel1, sel2", sel1, sel2
  master_xyz = hierarchy.select(sel1).atoms().extract_xyz()
  copy_xyz = hierarchy.select(sel2).atoms().extract_xyz()
  xyz = cur_ttg[2][0].elems * master_xyz + cur_ttg[2][1]
  rmsd2 = 0
  if copy_xyz.size() == xyz.size():
    rmsd2 = copy_xyz.rms_difference(xyz)
  return rmsd1, rmsd2


def get_info_from_match_dict(match_dict, key, chain):
  assert chain in key, "Mismatch between key and chain"
  [sel_1,sel_2,res_1,res_2,_,_,rmsd] = match_dict[key]
  # print "sel_1,sel_2,res_1,res_2,_,_,rmsd", sel_1,sel_2,res_1,res_2,rmsd
  if chain == key[0]:
    return sel_1, res_1, rmsd
  else:
    return sel_2, res_2, rmsd


def get_bool_selection_to_keep(big_selection, small_selection):
  """
  given 2 iselections (they are sorted), returns bool selection of size
  big selection showing what are the matches with small selection.
  Rather fast algorithm but may be beneficial to transfer to C++
  O(n+m), where n,m - sizes of selections
  """
  assert big_selection.size >= small_selection.size()
  result = flex.bool(big_selection.size(), False)
  i_in_big = 0
  i_in_small = 0
  size_small = small_selection.size()
  size_big = big_selection.size()
  n_matches = 0
  nw = 0
  while (i_in_big < size_big) and (i_in_small < size_small):
    if big_selection[i_in_big] == small_selection[i_in_small]:
      result[i_in_big] = True
      i_in_big += 1
      i_in_small += 1
      n_matches += 1
    elif big_selection[i_in_big] > small_selection[i_in_small]:
      i_in_small += 1
      nw += 1
    else:
      i_in_big += 1
  # this assert is optional, in general case it is not guaranteed that
  # all numbers from small selection are present in big selection.
  assert n_matches == size_small, "%d %d" % (n_matches, size_small)
  return result

def ncs_grouping_and_group_dict(match_dict, hierarchy):
  """
  The implementation of simplest way to do NCS grouping. Maximum one chain
  in selection.
  Do the job of minimal_master_ncs_grouping/minimal_ncs_operators_grouping
  and build_group_dict.
  """
  group_dict = {}
  # temp storages
  # [{chain id: key to match_dict}, {...} ... ]
  preliminary_ncs_groups = []
  for chain_pair, info in match_dict.iteritems():
    [sel_1,sel_2,res_1,res_2,r,t,rmsd] = info
    # if sel_1.size() < len_of_smallest_selection:
    #   len_of_smallest_selection = sel_1.size()
    #   key_with_smallest_selection = chain_pair
    # print "selection sizes:",chain_pair, sel_1.size(), sel_2.size()
    i_existing = None
    i_found_gr = None
    i = 0
    while i < 2 and i_existing is None:
      for n_gr, prel_ncs_group in enumerate(preliminary_ncs_groups):
        # print "checking ", chain_pair[i], prel_ncs_group.keys()
        if chain_pair[i] in prel_ncs_group.keys():
          i_existing = i
          i_found_gr = n_gr
          break
      i += 1
    if i_existing is None:
      assert i_found_gr is None
      # add new preliminary ncs group
      preliminary_ncs_groups.append({
          chain_pair[0]:chain_pair,
          chain_pair[1]:chain_pair})
    else:
      # add other chain to found ncs group
      preliminary_ncs_groups[i_found_gr][chain_pair[1-i_existing]] = chain_pair
  # print "preliminary_ncs_groups", preliminary_ncs_groups


  # now we need to just transform preliminary_ncs_groups using match_dict
  # into group_dict. This means that for every dict in preliminary_ncs_groups
  # we need to determine master, and find out rot and transl functions for all
  # the rest chains (selections). Master is going to be the first in
  # alphabetical order.

  group_id = 0
  tr_sn = 1

  for prel_gr_dict in preliminary_ncs_groups:
    # print "==============="
    sorted_gr_chains = sorted(prel_gr_dict.keys())

    # master should be the chain with minimal number of selected atoms
    # just to make it easier filter out the rest of chains
    # print "sorted_gr_chains", sorted_gr_chains
    # print "prel_gr_dict", prel_gr_dict
    min_n_atoms = 1e100
    master = None
    for ch in sorted_gr_chains:
      sel, _,_ = get_info_from_match_dict(match_dict, prel_gr_dict[ch], ch)
      if sel.size() < min_n_atoms:
        min_n_atoms = sel.size()
        master = ch
    assert master is not None
    # print "selected master first:", master

    # second option to master selection:
    # let's try to select common chain to be a master. I'm not sure that this
    # will be always possible though
    # also, we should try to determine the smallest selection for the master
    # chain straight away
    all_pairs = prel_gr_dict.values()
    left = set(all_pairs[0])
    # print "left", left
    # print "all_pairs", all_pairs
    for i in all_pairs[1:]:
      left = left & set(i)
    # should be 1 (a lot of chains) or 2 (if there only 2 chains)
    assert len(left) > 0
    # print "left", left
    if len(left) > 1:
      master = sorted(left)[0]
    else:
      master = left.pop()

    # print "selected master second:", master

    # selecting smallest master key - for no reason actually
    key_with_smallest_selection = None
    len_of_smallest_selection = 1e100
    for ch, key in prel_gr_dict.iteritems():
      master_sel, master_res, master_rmsd = get_info_from_match_dict(
              match_dict, key, master)
      if master_sel.size() < len_of_smallest_selection:
        len_of_smallest_selection = master_sel.size()
        key_with_smallest_selection = key
    # print "key_with_smallest_selection, len_of_smallest_selection",key_with_smallest_selection, len_of_smallest_selection

    assert master is not None
    assert master in key_with_smallest_selection, "%s, %s" % (master, key_with_smallest_selection)

    #
    # Let's do intersection of all master selection to determine
    # the minimum selection suitable to all copies.
    min_master_selection = None
    for ch, key in prel_gr_dict.iteritems():
      master_sel, master_res, master_rmsd = get_info_from_match_dict(
              match_dict, key, master)
      if min_master_selection is None:
        min_master_selection = master_sel
      else:
        min_master_selection = min_master_selection.intersection(master_sel)
    # print "size of min_master_selection", min_master_selection.size()

    #
    #
    # create a new group
    new_ncs_group = NCS_groups_container()
    tr = Transform(
        rotation=matrix.sqr([1,0,0,0,1,0,0,0,1]),
        translation=matrix.col([0,0,0]),
        serial_num=tr_sn,
        coordinates_present=True,
        ncs_group_id=group_id,
        rmsd=0)
    tr_sn += 1

    # master_sel, master_res, master_rmsd = get_info_from_match_dict(
    #     match_dict,key_with_smallest_selection, master)
    new_ncs_group.iselections.append([min_master_selection])
    new_ncs_group.residue_index_list.append([master_res])
    new_ncs_group.copies.append([master])
    new_ncs_group.transforms.append(tr)

    for ch_copy in sorted_gr_chains:
      # master_size = new_ncs_group.iselections[0][0].size()
      master_size = min_master_selection.size()
      if ch_copy == master:
        continue
      # check whether the master is the same XXX Do this later
      # otherwise just calculate new r,t
      # key = prel_gr_dict[ch_copy]
      copy_sel, copy_res, copy_rmsd = get_info_from_match_dict(
          match_dict,prel_gr_dict[ch_copy], ch_copy)
      from iotbx.pdb.atom_selection import selection_string_from_selection
      # print "master_sel:", selection_string_from_selection(hierarchy, master_sel)
      # print "copy_sel:", selection_string_from_selection(hierarchy, copy_sel)
      # print "master_sel:", list(master_sel)
      # print "copy_sel:", list(copy_sel)
      new_copy_sel = copy_sel
      new_master_sel = min_master_selection
      m_sel, m_res, m_rmsd = get_info_from_match_dict(
          match_dict, prel_gr_dict[ch_copy], master)
      # print "m_sel:", list(m_sel)
      if copy_sel.size() > min_master_selection.size():
        # clean copy sel
        # print "copy is bigger", copy_sel.size(), min_master_selection.size()
        # print "sizes:", master_sel.size(), m_sel.size()
        filter_sel = get_bool_selection_to_keep(
            big_selection=m_sel,
            small_selection=min_master_selection)
        new_copy_sel = copy_sel.select(filter_sel)
      elif copy_sel.size() < min_master_selection.size():
        # clean master sel and all other copies...
        # should never be the case anymore
        # print "master is bigger", copy_sel.size(), master_sel.size()
        # print "sizes:", master_sel.size(), m_sel.size()
        # print "master:", list(master_sel)
        filter_sel = get_bool_selection_to_keep(
            big_selection=master_sel,
            small_selection=m_sel)
        # print list(filter_sel)
        new_master_sel = master_sel.select(filter_sel)
        # print "len new_master_sel", len(new_master_sel)
        for i in range(len(new_ncs_group.iselections)):
          # print "new_ncs_group.iselections", new_ncs_group.iselections
          new_ncs_group.iselections[i] = [new_ncs_group.iselections[i][0].select(filter_sel)]
        master_sel = new_master_sel
        master_size = master_sel.size()
        assert 0

      # STOP()
      r,t,copy_rmsd = my_get_rot_trans(
          ph=hierarchy,
          master_selection=new_master_sel,
          copy_selection=new_copy_sel)
      tr = Transform(
          rotation=r,
          translation=t,
          serial_num=tr_sn,
          coordinates_present=True,
          ncs_group_id=group_id,
          rmsd=copy_rmsd)
      assert master_size == new_copy_sel.size(), "%d %d" % (master_size, new_copy_sel.size())
      new_ncs_group.iselections.append([new_copy_sel])
      new_ncs_group.residue_index_list.append([copy_res])
      new_ncs_group.copies.append([ch_copy])
      new_ncs_group.transforms.append(tr)
      tr_sn += 1
    group_dict[tuple(master)] = new_ncs_group
    master_size = new_ncs_group.iselections[0][0].size()
    for isel_arr in new_ncs_group.iselections[1:]:
      assert master_size ==isel_arr[0].size(), "%d %d" % (master_size, isel_arr[0].size().size())

    # print "new_ncs_group.ise", new_ncs_group.iselections
    # for isele_arr in new_ncs_group.iselections:
    #   print "final selections are:", list(isele_arr[0])
    # print "new_ncs_group.copies", new_ncs_group.copies
    # print "new_ncs_group.residue_index_list", new_ncs_group.residue_index_list
    group_id += 1

  # print "group_dict", group_dict
  # STOP()
  return group_dict


def minimal_master_ncs_grouping(match_dict, hierarchy):
  """
  XXX
  XXX still used in ncs_preprocess.py:validate_ncs_phil_groups()
  XXX
  Look for NCS groups with the smallest number of chains in the master copy
  This is not the minimal number of NCS operations

  Args:
    match_dict(dict):
      key:(chains_id_1,chains_id_2)
      val:[select_1,select_2,res_list_1,res_list_2,rot,trans,rmsd]
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
  # sorted_chain_groups_keys = match_dict.keys()
  # print "sorted_chain_groups_keys", sorted_chain_groups_keys
  # collect all chain IDs
  chain_ids = {c_id for key in sorted_chain_groups_keys for c_id in key}
  chain_ids = sorted(chain_ids)
  # print "chain_ids", chain_ids
  # STOP()
  # for each group in chain_groups find the transforms and group them
  transform_to_group = {}
  chains_in_groups = set()
  tr_sn = 0
  cache = hierarchy.atom_selection_cache()
  #
  for (master_id, copy_id) in sorted_chain_groups_keys:
    [sel_1,sel_2,res_1,res_2,r,t,rmsd] = match_dict[master_id,copy_id]
    # check that master and copy are not a copies in another group
    # print "tr_sn:", tr_sn
    # print "  master, copy, chains_in_groups:", master_id, copy_id, chains_in_groups
    # print "  transform_to_group",
    # for k, v in transform_to_group.iteritems():
    #   print "    ", k, v[:2]
    if (master_id in chains_in_groups) or (copy_id in chains_in_groups):
      match_dict.pop((master_id, copy_id),None)
      # print "!!! Popping", master_id, copy_id
      continue
    # check if the chains related by ncs operation
    tr_num, is_transpose = find_same_transform(r,t,transform_to_group)
    # print "tr_num", tr_num, is_transpose
    if is_transpose and ((master_id in chains_in_groups) or (copy_id in chains_in_groups)):
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
    if tr_num and temp_master not in transform_to_group[tr_num][0]:
      # Here we want to make sure that this update won't violate
      # order of atoms, if not, we are changing the order...
      rmsd1, rmsd2 = get_rmsds2(
          master_xyz=hierarchy.select(sel_1).atoms().extract_xyz(),
          copy_xyz=hierarchy.select(sel_2).atoms().extract_xyz(),
          cur_ttg=transform_to_group[tr_num])
      # rmsd1, rmsd2 = get_rmsds(
      #     hierarchy,
      #     cache,
      #     transform_to_group[tr_num],
      #     temp_master,
      #     temp_copy)
      # print "RMSDs:", rmsd1, rmsd2
      if rmsd1 > 10 or rmsd2 > 10:
        tr_sn += 1
        transform_to_group[tr_sn] = [[temp_master],[temp_copy],(r,t)]
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
    [_,k] = master_size.pop(0)
    tr = transform_to_group[k]
    # make sure ALL copies from transform need to be added
    if not bool(set(tr[1]) - chain_left_to_add):
      transform_to_use.add(k)
      chain_left_to_add -= set(tr[1])
      chain_left_to_add -= set(tr[0])
  # delete all transforms that are not being used
  tr_sns = transform_to_group.keys()
  for tr_sn in tr_sns:
    if not tr_sn in transform_to_use:
      transform_to_group.pop(tr_sn)
  return transform_to_group,match_dict

def remove_masters_if_appear_in_copies(transform_to_group):
  """
  remove masters if appear in copies of the same transform
  and of another transform

  Args:
    transform_to_group (dict):
      keys: tuple of master chain IDs
      values: NCS_groups_container objects
    chain_ids (list): sorted list of chains IDs

  Returns:
    updated transform_to_group
  """
  # clean transforms
  for k,v in transform_to_group.iteritems():
    masters = v[0]
    copies = v[1]
    keep_items = []
    for i,ch_id in enumerate(masters):
      if not (ch_id in copies):
        keep_items.append(i)
    transform_to_group[k][0] = [masters[x] for x in keep_items]
    transform_to_group[k][1] = [copies[x] for x in keep_items]
  # clean other transforms once the transforms are clean
  all_copies = []
  remove_items = []
  keys = sorted(transform_to_group,key=lambda k:transform_to_group[k][0])
  for k in keys:
    copy_id = transform_to_group[k][1]
    master_id = transform_to_group[k][0]
    if (copy_id in all_copies) or (master_id in all_copies):
      remove_items.append(k)
    else:
      all_copies.append(copy_id)
  for k in remove_items:
    transform_to_group.pop(k)
  return transform_to_group

def remove_overlapping_selection(transform_to_group,chain_ids):
  """
  remove overlapping selection and keep largest groups

  Args:
    transform_to_group (dict):
      keys: tuple of master chain IDs
      values: NCS_groups_container objects
    chain_ids (list): sorted list of chains IDs

  Returns:
    updated transform_to_group
  """
  # remove overlapping selection, keep largest groups
  master_size = [[len(v[0]),k] for k,v in transform_to_group.iteritems()]
  master_size.sort()
  chain_left_to_add = set(chain_ids)
  transform_to_use = set()
  chains_in_master = set()
  chains_in_copies = set()
  while bool(chain_left_to_add) and bool(master_size):
    [_,k] = master_size.pop()
    [masters,copies,_] = transform_to_group[k]
    # check that all chains in copies still need to be added
    test1 = len(set(copies).intersection(chain_left_to_add)) == len(copies)
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
  return transform_to_group

def clean_singles(transform_to_group,chains_in_groups):
  """
  When a NCS copy of a transform contains a single chains, but those chains are
  also a part of a transform with several chains in the NCS group,
  remove the single transform

  Args:
    transform_to_group (dict):
      keys: tuple of master chain IDs
      values: NCS_groups_container objects
    chains_in_groups (set): Set of chain IDs that are in some NCS group

  Returns:
    updated transform_to_group
  """
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
  return transform_to_group

def build_group_dict(transform_to_group,match_dict,chains_info):
  """
  find all transforms with common masters and organize the chains
  in the same order as the master and build groups transform dictionary

  Args:
    transform_to_group (dict):
      key: temporary transform number
      values:[masters],[copies],(rotation,translation)]
    match_dict(dict):
      key:(chains_id_1,chains_id_2)
        chain_ID (str)
      val:[selection_a,selection_b,res_list_a,res_list_b,rot,trans,rmsd]
          select_1/select_2 (flex.size_t)
          residue_index_list_1/2 (list): matching residues indices
          rot (matrix obj): rotation matrix
          tran (matrix obj): translation vector
          rmsd (float)
    chains_info : object containing
      chains (str): chain IDs OR selections string
      res_name (list of str): list of residues names
      resid (list of str): list of residues sequence number, resid
      atom_names (list of list of str): list of atoms in residues
      atom_selection (list of list of list of int): the location of atoms in ph
      chains_atom_number (list of int): list of number of atoms in each chain

  Returns:
    group_dict (dict):
      keys: tuple of master chain IDs
      values: NCS_groups_container objects
  """
  group_dict = {}
  group_id = 0
  tr_sn = 0
  adjust_key_lists = set()
  for k,v in transform_to_group.iteritems():
    [masters,copies,_] = v
    key = tuple(masters)
    m_isel_list,c_isel_list,res_l_m,res_l_c,rmsd,r,t = \
      collect_info(masters,copies,match_dict)
    if group_dict.has_key(key):
      # update existing group master with minimal NCS selection
      c_isel_list = update_group_dict(
        group_dict,key,adjust_key_lists,m_isel_list,c_isel_list)
      tr_sn += 1
      tr = Transform(
        rotation=r,
        translation=t,
        serial_num=tr_sn,
        coordinates_present=True,
        ncs_group_id=group_id,
        rmsd=rmsd)
      group_dict[key].iselections.append(c_isel_list)
      group_dict[key].residue_index_list.append(res_l_c)
      group_dict[key].copies.append(copies)
      group_dict[key].transforms.append(tr)
    else:
      # create a new group
      tr_sn += 1
      group_id += 1
      _,_,res_l_m,res_l_c,rmsd,r,t = collect_info(masters,copies,match_dict)
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
      new_ncs_group.iselections.append(m_isel_list)
      new_ncs_group.residue_index_list.append(res_l_m)
      new_ncs_group.copies.append(masters)
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
      new_ncs_group.iselections.append(c_isel_list)
      new_ncs_group.residue_index_list.append(res_l_c)
      new_ncs_group.copies.append(copies)
      new_ncs_group.transforms.append(tr)
      #
      group_dict[key] = new_ncs_group
  # adjust residue_index_list according to the new iselections
  if adjust_key_lists:
    update_res_list(group_dict,chains_info,adjust_key_lists)
  return group_dict

def update_group_dict(group_dict,key,adjust_key_lists,m_isel_list,c_isel_list):
  """
  Updates group_dict with the minimal common selection for all NCS copies
  Updates adjust_key_lists, a set that collect all records that where
  adjusted

  Args:
    group_dict (dict):
      keys: tuple of master chain IDs
      values: NCS_groups_container objects
    key : a key of group_dict
    adjust_key_lists (set): set of key
    m_isel_list,c_isel_list (list of flex.size_t): list of master and copy
     iselections

  Returns:
    new_copy (flex.size_t): selection on the new NCS copy
  """
  new_master = []
  new_copy = []
  for i in xrange(len(m_isel_list)):
    m_sel = m_isel_list[i]
    c_sel = c_isel_list[i]
    current_master = group_dict[key].iselections[0][i]
    # for each chain, check if the master have the same selection
    if current_master.size() != m_sel.size():
      adjust_master = True
    else:
      # make sure the same atoms are selected
      temp = (m_sel == current_master)
      adjust_master = (temp.count(False) != 0)
    if adjust_master:
      # find atoms that are only in the old or new master and remove them
      remove_from_new = set(m_sel) - set(current_master)
      remove_from_old = set(current_master) - set(m_sel)
      #
      sel_to_keep = selection_to_keep(m_sel,remove_from_new)
      m_sel = m_sel.select(sel_to_keep)
      c_sel = c_sel.select(sel_to_keep)
      #
      new_master.append(flex.size_t(m_sel))
      new_copy.append(flex.size_t(c_sel))
      adjust_key_lists.add(key)
      # update all existing copies
      n = len(group_dict[key].iselections)
      sel_to_keep = selection_to_keep(current_master,remove_from_old)
      for j in range(n):
        isel = group_dict[key].iselections[j][i]
        isel = isel.select(sel_to_keep)
        group_dict[key].iselections[j][i] = isel
    else:
      new_master.append(m_sel)
      new_copy.append(c_sel)
  return new_copy

def selection_to_keep(m_sel,remove_selection):
  """
  Get the selection of the selection to keep from flex.size_t array

  Args:
    m_sel (flex.size_t): the array from which we want to remove items
    remove_selection (flex.size_t): the values of items to remove
  Returns:
    sel_to_keep (flex.size_t)
  """
  msel = list(m_sel)
  # find values location
  remove_sel = {msel.index(x) for x in remove_selection}
  sel_to_keep = set(range(m_sel.size())) - remove_sel
  sel_to_keep = list(sel_to_keep)
  sel_to_keep.sort()
  return flex.size_t(sel_to_keep)

def update_res_list(group_dict,chains_info,group_key_lists):
  """
  Remove residues from residues list according the the atom selection
  (Atom selection might have changed to reflect the smallest common set of
  NCS related atoms in all NCS copies)

  Args:
    group_dict (dict):
      keys: tuple of master chain IDs
      values: NCS_groups_container objects
    chains_info : object containing
      chains (str): chain IDs OR selections string
      res_name (list of str): list of residues names
      resid (list of str): list of residues sequence number, resid
      atom_names (list of list of str): list of atoms in residues
      atom_selection (list of list of list of int): the location of atoms in ph
      chains_atom_number (list of int): list of number of atoms in each chain
    group_key_lists (list): set of group_dict keys to update
  """
  res_list = []
  # process all groups
  for key in group_key_lists:
    gr = group_dict[key]
    # iterate over the NCS copies in the group
    copies = []
    iselections = []
    transforms = []
    for i,ch_keys in enumerate(gr.copies):
      c_res_list = []
      atoms_in_copy = set()
      {atoms_in_copy.update(x) for x in gr.iselections[i]}
      # keep only none-zero copies
      if len(atoms_in_copy) > 0:
        copies.append(gr.copies[i])
        iselections.append(gr.iselections[i])
        transforms.append(gr.transforms[i])
      else: continue
      copy_res_lists = gr.residue_index_list[i]
      # iterate over the chains in each NCS group
      n_ch = len(ch_keys)
      for i_ch in range(n_ch):
        ch_key = ch_keys[i_ch]
        ch_res_list = copy_res_lists[i_ch]
        ch_info = chains_info[ch_key]
        c_res = []
        for res_num in ch_res_list:
          # iterate over residues and add them if they are in atoms_in_copy
          atoms_in_rs = set(ch_info.atom_selection[res_num])
          if bool(atoms_in_rs.intersection(atoms_in_copy)):
            # if some atoms in the residue present, include residue
            c_res.append(res_num)
        c_res_list.append(c_res)
      res_list.append(c_res_list)
    if len(res_list) > 0:
      group_dict[key].residue_index_list = res_list
      group_dict[key].copies = copies
      group_dict[key].iselections = iselections
      group_dict[key].transforms = transforms
    else:
      group_dict.pop(key,None)

def collect_info(sorted_masters,copies,match_dict):
  """
  Combine iselection of all chains in NCS master and NCS copy to a single
  iselelction

  Args:
    sorted_masters (list): list of master chain IDs
    copies (list): list of copies chain IDs
    match_dict(dict):
      key:(chains_id_1,chains_id_2)
        chain_ID (str)
      val:[selection_a,selection_b,res_list_a,res_list_b,rot,trans,rmsd]
          select_1/select_2(flex.size_t)
          residue_index_list_1/2 (list): matching residues indices
          rot (matrix obj): rotation matrix
          tran (matrix obj): translation vector
          rmsd (float)

  Returns:
    m_sel,c_sel (list of flex.size_t): list of selected atoms per chain
    m_res,c_res (list of lists): list of residue number per chain
    worst_rmsd (float): worst rmsd of the matching chains pairs
    r (3x3 matrix): rotation matrix
    t (3x1 matrix): translation vector
  """
  m_sel = []
  c_sel = []
  c_res = []
  m_res = []
  worst_rmsd = 0.0
  for key in zip(sorted_masters,copies):
    [sel_1,sel_2,res_1,res_2,_,_,rmsd] = match_dict[key]
    m_sel.append(sel_1)
    c_sel.append(sel_2)
    m_res.append(res_1)
    c_res.append(res_2)
    worst_rmsd = max(worst_rmsd,rmsd)
  [_,_,_,_,r,t,_] = match_dict[sorted_masters[0],copies[0]]
  return m_sel,c_sel,m_res,c_res,worst_rmsd,r,t

def clean_chain_matching(chain_match_list,ph,
                         max_rmsd=10.0,
                         exclude_misaligned_residues=False,
                         match_radius=4.0,similarity_threshold=0.95):
  """
  Remove all bad matches from chain_match_list

  Args:
    ph (object): hierarchy
    chain_match_list (list): list of
      [chain_ID_1, chain_ID_2, sel_1, sel_2,res_m/res_c similarity]
      chain_ID (str), sel_1/2 (list of lists)
      res_m/res_c (lists): indices of the aligned components
      similarity (float): similarity between chains
    max_rmsd (float): limit of rms difference chains
    exclude_misaligned_residues (bool): check and exclude individual residues
      alignment quality
    match_radius (float): max allow distance difference between pairs of matching
      atoms of two residues
    similarity_threshold (float): min similarity between matching chains

  Returns:
    match_dict(dict): key:(chains_id_a,chains_id_b)
                      val:[selection_a,selection_b,
                           res_list_a,res_list_b,rot,trans,rmsd]
  """
  # remove all non-matching pairs, where similarity == 0
  match_list = [x for x in chain_match_list if x[4] > 0]
  # keep only best (or 95% of best) matches
  best_matches = {}
  # Get rmsd
  match_dict = {}
  for match in match_list:
    [ch_a_id,ch_b_id,list_a,list_b,res_list_a,res_list_b,similarity] = match
    update_match_dicts(
      best_matches,match_dict,ch_a_id,ch_b_id,similarity,similarity_threshold)
    sel_a = make_selection_from_lists(list_a)
    sel_b = make_selection_from_lists(list_b)
    other_sites = ph.select(sel_a).atoms().extract_xyz()
    ref_sites = ph.select(sel_b).atoms().extract_xyz()
    lsq_fit_obj = superpose.least_squares_fit(
      reference_sites = ref_sites,
      other_sites     = other_sites)
    r = lsq_fit_obj.r
    t = lsq_fit_obj.t
    # todo: find r_2*A = r*A + t (where the translation is zero)
    # use B = r*A + t, r_2*A = B , r_2 = B*A.inverse()
    other_sites_best = lsq_fit_obj.other_sites_best_fit()
    rmsd = round(ref_sites.rms_difference(other_sites_best),4)
    if rmsd <= max_rmsd:
      if exclude_misaligned_residues:
        # get the chains atoms and convert selection to flex bool
        sel_a,sel_b,res_list_a,res_list_b,ref_sites,other_sites_best = \
          remove_far_atoms(
            list_a, list_b,
            res_list_a,res_list_b,
            ref_sites,lsq_fit_obj.other_sites_best_fit(),
            match_radius=match_radius)
      if sel_a.size() > 0:
        match_dict[ch_a_id,ch_b_id]=[sel_a,sel_b,res_list_a,res_list_b,r,t,rmsd]
  return match_dict

def remove_far_atoms(list_a, list_b,
                     res_list_a,res_list_b,
                     ref_sites,other_sites,
                     match_radius=4.0):
  """
  When comparing lists of matching atoms, remove residues where some atoms are
  are locally misaligned, for example when matching residues are
  perpendicular to each other rather than being close to parallel.

  The criteria used:
  For each matching residues, the difference between distance of farthest
  matching atoms pair and the distance of closest pair mast be < match_radius

  Args:
    list_a, list_a (list of list): list of residues atoms
    res_list_a,res_list_b (list): list of residues in chains
    ref_sites,other_sites (flex.vec3): atoms coordinates
    match_radius (float): max allow distance difference

  Returns:
    Updated arguments:
      sel_a,sel_b,
      res_list_a_new,res_list_b_new,
      ref_sites_new,other_sites_new
  """
  # check every residue for consecutive distance
  res_list_a_new = []
  res_list_b_new = []
  ref_sites_new = flex.vec3_double([])
  other_sites_new = flex.vec3_double([])
  sel_a = flex.size_t([])
  sel_b = flex.size_t([])
  current_pos = 0
  for i in xrange(len(res_list_a)):
    # find the matching atoms form each residue (work on small sections)
    res_len = list_a[i].size()
    res_ref_sites = ref_sites[current_pos:current_pos+res_len]
    res_other_sites = other_sites[current_pos:current_pos+res_len]
    current_pos += res_len
    xyz_diff = abs(res_ref_sites.as_double() - res_other_sites.as_double())
    (min_d,max_d,_) = xyz_diff.min_max_mean().as_tuple()
    if (max_d - min_d) <= match_radius:
      ref_sites_new.extend(res_ref_sites)
      other_sites_new.extend(res_other_sites)
      sel_a.extend(list_a[i])
      sel_b.extend(list_b[i])
      res_list_a_new.append(res_list_a[i])
      res_list_b_new.append(res_list_b[i])
  return sel_a,sel_b,res_list_a_new,res_list_b_new,ref_sites_new,other_sites_new

def update_match_dicts(best_matches,match_dict,
                       ch_a_id,ch_b_id,similarity,
                       similarity_threshold=0.95):
  """
  Updates the best_matches dictionaries best_matches,match_dict to keep only
  matches that are at least 95% of best match and

  Args:
    best_matches (dict):
      key: chain_ID (str)
      Val: [[similarity,(chain_ID,chain_ID2)],[similarity,(chain_ID3,chain_ID)]]
    match_dict (dict):
      key: (chain_id_a,chain_id_b)
      val: [selection_1,selection_2,res_list_1,res_list_2,rot,trans,rmsd]
    ch_a_id,ch_b_id (str): chain IDs
    similarity (float): similarity between chains
    similarity_threshold (float): min similarity between matching chains
      Note that smaller value cause more chains to be grouped together and
      can lower the number of common residues
  """
  records_to_remove = set()
  for ch_id in [ch_a_id,ch_b_id]:
    if best_matches.has_key(ch_id):
      temp_rec = []
      # records with largest similarity are last
      max_sim = best_matches[ch_id][-1][0]
      if similarity > max_sim:
        for s,(a,b) in best_matches[ch_id]:
          if s/similarity >= similarity_threshold:
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

def find_same_transform(r,t,transforms):
  """
  Check if the rotation r and translation t exist in the transform dictionary.
  Note that there can be both inverse and regular match. Return the
  non-transpose if exist.

  Args:
    r (matrix.sqr): rotation
    t (matrix.col): translation
    transforms (dict): dictionary of all transforms

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
    is_the_same, is_transpose_flag = is_same_transform(r, t, rr, tt)
    if is_the_same:
      if is_transpose_flag:
        # when transpose is found, keep it but continue the search
        tr_num = k
        is_transpose  = True
      else:
        # found non-transform match
        return k, False
  return tr_num, is_transpose

def get_sequence_from_array(arr):
  from iotbx.pdb import amino_acid_codes
  aa_3_as_1 = amino_acid_codes.one_letter_given_three_letter
  res = ""
  for r in arr:
    res += aa_3_as_1.get(r)
  return res

def search_ncs_relations(ph=None,
                         chains_info = None,
                         min_contig_length=10,
                         min_percent=0.95,
                         write=False,
                         log=None,
                         check_atom_order=False,
                         allow_different_size_res=True,
                         use_minimal_master_ncs=True,
                         ignore_chains=None):
  """
  Search for NCS relations between chains or parts of chains, in a protein
  hierarchy

  Args:
    ph (object): hierarchy
    chains_info (dict): values are object containing
      chains (str): chain IDs OR selections string
      res_name (list of str): list of residues names
      resid (list of str): list of residues sequence number, resid
      atom_names (list of list of str): list of atoms in residues
      atom_selection (list of list of list of int): the location of atoms in ph
      chains_atom_number (list of int): list of number of atoms in each chain
    min_contig_length (int): segments < min_contig_length rejected
    min_percent (float): Threshold for similarity between chains.
      similarity define as:
      (number of matching res) / (number of res in longer chain)
    write (bool): when true, write ncs search messages to log
    check_atom_order (bool): check atom order in matching residues.
        When False, matching residues with different number of atoms will be
        excluded from matching set
    use_minimal_master_ncs (bool): use maximal or minimal common chains
        in master ncs groups (when True, search does not need to be done all
        all chains -> can be significantly faster on large structures)
    allow_different_size_res (bool): keep matching residue with different
      number of atoms
    ignore_chains (set of str): set of chain IDs to exclude

  Returns:
    msg (str): message regarding matching residues with different atom number
    chain_match_list (list): list of
      [chain_ID_1,chain_ID_2,sel_1,sel_2,res_sel_m, res_sel_c,similarity]
      chain_ID (str), sel_1/2 (flex.size_t),
      res_sel_m/c (lists of lists): indices of the aligned components
      similarity (float): similarity between chains
    We use sel_2 to avoid problems when residues have different number of atoms
  """
  # print "searching ncs relations..."
  if not log: log = sys.stdout
  if not chains_info:
    assert bool(ph)
    chains_info = get_chains_info(ph,ignore_chains=ignore_chains)
  # collect all chain IDs
  chain_match_list = []
  msg = ''
  if use_minimal_master_ncs:
    sorted_ch = sorted(chains_info)
  else:
    sorted_ch = sort_by_dist(chains_info)
  n_chains = len(sorted_ch)
  chains_in_copies = set()
  # print "sorted_ch", sorted_ch
  # loop over all chains
  for i in xrange(n_chains-1):
    if use_minimal_master_ncs:
      if (sorted_ch[i] in chains_in_copies): continue
      # update sorted_ch list according to distance between master copies
      sorted_ch = update_chain_ids_search_order(
        chains_info,sorted_ch,chains_in_copies,i)
    m_ch_id = sorted_ch[i]
    master_n_res = len(chains_info[m_ch_id].res_names)
    seq_m = chains_info[m_ch_id].res_names
    if master_n_res == 0: continue
    # get residue lists for master
    for j in xrange(i+1,n_chains):
      c_ch_id = sorted_ch[j]
      copy_n_res = len(chains_info[c_ch_id].res_names)
      frac_d = min(copy_n_res,master_n_res)/max(copy_n_res,master_n_res)
      # print "  copy_n_res,master_n_res", copy_n_res,master_n_res
      if frac_d < min_percent:
        if (min_percent == 1):
          msg = 'NCS copies are not identical'
          break
        else:
          # print "Strange exit"
          continue
      seq_c = chains_info[c_ch_id].res_names
      # get residue lists for copy
      res_sel_m, res_sel_c, similarity = res_alignment(
        seq_a=seq_m,seq_b=seq_c,
        min_contig_length=min_contig_length,
        min_percent=min_percent)
      sel_m, sel_c,res_sel_m,res_sel_c,new_msg = get_matching_atoms(
        chains_info,m_ch_id,c_ch_id,res_sel_m,res_sel_c,
        check_atom_order=check_atom_order,
        allow_different_size_res=allow_different_size_res)
      msg += new_msg
      if res_sel_m:
        # add only non empty matches
        rec = [m_ch_id,c_ch_id,sel_m,sel_c,res_sel_m,res_sel_c,similarity]
        chain_match_list.append(rec)
      # Collect only very good matches, to allow better similarity search
      # print "similarity", similarity
      if similarity > 0.9: chains_in_copies.add(c_ch_id)
  if write and msg:
    print >> log,msg
  if (min_percent == 1) and msg:
    # must be identical
    raise Sorry('NCS copies are not identical')
  return chain_match_list

def res_alignment(seq_a, seq_b,
                  min_contig_length=10,
                  min_percent=0.95):
  """
  Align two lists of strings (All lower case characters).

  Penalize misalignment, gaps, when contiguous section is less
  than min_contig_length and when the
  Do not give any points for alignment. (score only change when "bad" things
  happens)

  Args:
    seq_a, seq_b (lists of str): list of characters to compare
    min_percent (float): min percent of similarity between hierarchies
      similarity define as:
      (number of matching res) / (number of res in longer chain)
    min_contig_length (int): domain < min_contig_length rejected

  Returns:
    aligned_sel_a (list): the indices of the aligned components of seq_a
    aligned_sel_b (list): the indices of the aligned components of seq_b
    similarity (float): actual similarity between hierarchies
  """
  # print "inputs in res_alignment, seq_a:", seq_a
  # print "inputs in res_alignment, seq_b:", seq_b
  # print "inputs in res_alignment:", min_contig_length, min_percent
  # print "seq_a", get_sequence_from_array(seq_a)
  # print "====="
  # print "seq_b", get_sequence_from_array(seq_b)

  a = len(seq_a)
  b = len(seq_b)
  # print "a,b", a,b
  # Check for the basic cases
  if (a == 0) or (b == 0): return [],[],0
  if seq_a == seq_b: return range(a),range(a),1.0
  min_matches = min(min_contig_length,min(a,b))
  misalign_penalty = 100000
  # limit the number of mis-alignments
  max_mis_align = int((1 - min_percent) * max(a,b))
  # Starting score according to the required similarity
  score = max_mis_align + 1
  # build score matrix
  R = initialize_score_matrix(row=a,col=b,max_score=score)
  # populate score matrix
  for j in xrange(1,b + 1):
    for i in xrange(1,a + 1):
      not_aligned = (seq_a[i-1].lower() != seq_b[j-1].lower())
      s1 = R[i-1,j-1].score - misalign_penalty * not_aligned
      s2,mc2 = gap_score(R[i-1,j],min_matches)
      s3,mc3 = gap_score(R[i,j-1],min_matches)
      if (s1 >= s2) and (s1 >= s3):
        s = s1
        i_temp,j_temp = i-1,j-1
        match_count = R[(i_temp,j_temp)].match_count + 1
        consecutive_matches = R[(i_temp,j_temp)].consecutive_matches + 1
      elif (s2 >= s1) and (s2 >= s3):
        i_temp,j_temp = i-1,j
        s = s2
        match_count = mc2
        consecutive_matches = 0
      else:
        i_temp,j_temp = i,j-1
        s = s3
        match_count = mc3
        consecutive_matches = 0
      R[i,j].score = s
      R[i,j].origin = (i_temp,j_temp)
      R[i,j].match_count = match_count
      R[i,j].consecutive_matches = consecutive_matches
  #
  aligned_sel_a, aligned_sel_b, similarity = get_matching_res_indices(
    R=R,row=a,col=b,min_percent=min_percent,min_contig_length=min_matches)
  # print "aligned_sel_a, aligned_sel_b, similarity", list(aligned_sel_a), list(aligned_sel_b), similarity
  return aligned_sel_a, aligned_sel_b, similarity

def gap_score(r,min_matches):
  """
  Calculate the score in case of a gap, taking into account consecutive matches
  and the minimum matches required.
  If we have a short matching section, we need to adjust the match_count,
  so not to count it

  Args:
    r (Score_record object)
    min_matches (int): minimum matches required between to sequences

  Returns:
    s (int): score
    min_matches (int): adjusted match count
  """
  s = r.score - 1
  match_count = r.match_count
  if r.consecutive_matches < min_matches:
    s -= r.consecutive_matches
    match_count -= r.consecutive_matches
  return s,match_count

def get_matching_res_indices(R,row,col,min_percent,min_contig_length):
  """
  Trace back best solution and collect the indices of matching pairs

  Args:
    R (list of lists of int): the score matrix
    row (int): number of rows
    col (int): number of columns
    min_percent (float): min percent of similarity between hierarchies
      similarity define as:
      (number of matching res) / (number of res in longer chain)
    min_contig_length (int): domain < min_contig_length rejected

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
    return flex.size_t([]), flex.size_t([]),0
  #
  similarity = R[i_max,j_max].match_count / max(i_max,j_max)
  if similarity < min_percent:
    # chains are to different, return empty arrays
    return flex.size_t([]), flex.size_t([]), 0
  # retrieve the matching sequences from the score matrix
  stop_test = 1
  domain_length = 0
  temp_sel_a = []
  temp_sel_b = []
  while stop_test > 0:
    if R[i_max,j_max].origin == (i_max - 1, j_max - 1):
      temp_sel_a.append(i_max - 1)
      temp_sel_b.append(j_max - 1)
      domain_length +=1
    else:
      # domain ended, if it is long enough keep it
      if domain_length >= min_contig_length:
        sel_a.extend(temp_sel_a)
        sel_b.extend(temp_sel_b)
      # restart domain length counting
      domain_length = 0
      temp_sel_a = []
      temp_sel_b = []
    i_max, j_max = R[i_max,j_max].origin
    # i_max or j_max reach to zero -> stop
    stop_test = i_max * j_max
  # check if the last domain is long enough to keep
  if domain_length >= min_contig_length:
    sel_a.extend(temp_sel_a)
    sel_b.extend(temp_sel_b)
  sel_a.reverse()
  sel_b.reverse()
  assert len(sel_a) == len(sel_b)
  return sel_a, sel_b, similarity

def get_matching_atoms(chains_info,a_id,b_id,res_num_a,res_num_b,
                       check_atom_order=False,
                       allow_different_size_res=True):
  """
  Get selection of matching chains, match residues atoms
  We keep only residues with continuous matching atoms

  Residues with alternative locations and of different size are excluded

  Args:
    chains_info (object): object containing
      chains (str): chain IDs or selections string
      res_name (list of str): list of residues names
      resid (list of str): list of residues sequence number, resid
      atom_names (list of list of str): list of atoms in residues
      atom_selection (list of list of list of int): the location of atoms in ph
      chains_atom_number (list of int): list of number of atoms in each chain
    a_id,b_id (str): Chain IDs
    res_num_a/b (list of int): indices of matching residues position
    check_atom_order (bool): check atom order in matching residues.
        When False, matching residues with different number of atoms will be
        excluded from matching set
    allow_different_size_res (bool): keep matching residue with different
      number of atoms

  Returns:
    sel_a/b (list of lists): matching atoms selection
    res_num_a/b (list of int): updated res_num_a/b
    msg (str): message regarding matching residues with different atom number
  """
  sel_a = []
  sel_b = []
  # check if any of the residues has alternate locations
  a_altloc = bool(chains_info[a_id].no_altloc)
  if a_altloc:
    a_altloc = chains_info[a_id].no_altloc.count(False) > 0
  b_altloc = bool(chains_info[b_id].no_altloc)
  if b_altloc:
    b_altloc = chains_info[b_id].no_altloc.count(False) > 0
  test_altloc = a_altloc or b_altloc
  #
  res_num_a_updated = []
  res_num_b_updated = []
  residues_with_different_n_atoms = []
  for (i,j) in zip(res_num_a,res_num_b):
    # iterate over atoms in residues
    sa = flex.size_t(chains_info[a_id].atom_selection[i])
    sb = flex.size_t(chains_info[b_id].atom_selection[j])
    dif_res_size = sa.size() != sb.size()
    atoms_names_a = chains_info[a_id].atom_names[i]
    atoms_names_b = chains_info[b_id].atom_names[j]
    resid_a = chains_info[a_id].resid[i]
    force_check_atom_order = dif_res_size and allow_different_size_res
    altloc = False
    if test_altloc:
      if a_altloc:
        altloc |= (not chains_info[a_id].no_altloc[i])
      if b_altloc:
        altloc |= (not chains_info[b_id].no_altloc[j])
    if check_atom_order or force_check_atom_order:
      # select only atoms that exist in both residues
      atoms_a,atoms_b,similarity = res_alignment(
        seq_a=atoms_names_a, seq_b=atoms_names_b,
        min_contig_length=1,min_percent=0.2)
      # get the number of the atom in the chain
      sa = flex.size_t(atoms_a) + sa[0]
      sb = flex.size_t(atoms_b) + sb[0]
    if dif_res_size or altloc:
      residues_with_different_n_atoms.append(resid_a)
      if (not allow_different_size_res) or altloc:
        sa = flex.size_t([])
        sb = flex.size_t([])
    # keep only residues with continuous matching atoms
    if sa.size() != 0:
      res_num_a_updated.append(i)
      res_num_b_updated.append(j)
      sel_a.append(sa)
      sel_b.append(sb)
  if residues_with_different_n_atoms:
    problem_res_nums = {x.strip() for x in residues_with_different_n_atoms}
    msg = "NCS related residues with different number of atoms, selection"
    msg += a_id + ':' + b_id + '\n['
    msg += ','.join(problem_res_nums) + ']\n'
  else:
    msg = ''
  return sel_a,sel_b,res_num_a_updated,res_num_b_updated,msg

def make_selection_from_lists(sel_list):
  """ Convert a list of lists to flex.size_t selection array  """
  sel_list_extended = [x for y in sel_list for x in y]
  sel_set = set(sel_list_extended)
  assert len(sel_list_extended) == len(sel_set)
  sel_list_extended.sort()
  return flex.size_t(sel_list_extended)

def get_chains_info(ph,selection_list=None,exclude_water=True,
                    ignore_chains=None):
  """
  Collect information about chains or segments of the hierarchy according to
  selection strings
  Exclude water atoms
  When there are alternate conformations, we use the first one

  Args:
    ph : protein hierarchy
    selection_list (list of str): specific selection of hierarchy segments
    ignore_chains (set of str): set of chain IDs to exclude

  Returns:
    chains_info (dict): values are object containing
      chains (str): chain IDs OR selections string
      res_name (list of str): list of residues names
      resid (list of str): list of residues sequence number, resid
      atom_names (list of list of str): list of atoms in residues
      atom_selection (list of list of list of int): the location of atoms in ph
      chains_atom_number (list of int): list of number of atoms in each chain
    exclude_water (bool): exclude water
  """
  use_chains = not bool(selection_list)
  if not ignore_chains: ignore_chains = set()
  #
  chains_info =  {}
  cache = ph.atom_selection_cache().selection
  if use_chains:
    model  = ph.models()[0]
    # build chains_info from hierarchy
    for ch in model.chains():
      if ch.id in ignore_chains: continue
      if not chains_info.has_key(ch.id):
        chains_info[ch.id] = Chains_info()
        ph_sel = ph.select(cache('chain ' + ch.id))
        coc = flex.vec3_double([ph_sel.atoms().extract_xyz().mean()])
        chains_info[ch.id].center_of_coordinates = coc
      chains_info[ch.id].chains_atom_number += ch.atoms().size()
      resids = chains_info[ch.id].resid
      res_names = chains_info[ch.id].res_names
      atom_names = chains_info[ch.id].atom_names
      atom_selection = chains_info[ch.id].atom_selection
      no_altloc = chains_info[ch.id].no_altloc
      # check for alternative conformers (also when chains are split)
      if len(ch.conformers()) > 1:
        # process cases with alternative locations
        conf = ch.conformers()[0]
        for res in conf.residues():
          x = res.resname
          if exclude_water and (x.lower() == 'hoh'): continue
          resids.append(res.resid())
          res_names.append(x)
          atoms = res.atoms()
          atom_names.append(list(atoms.extract_name()))
          atom_selection.append(list(atoms.extract_i_seq()))
          no_altloc.append(res.is_pure_main_conf)
      else:
        for res in ch.residue_groups():
          for atoms in res.atom_groups():
            x = atoms.resname
            if exclude_water and (x.lower() == 'hoh'): continue
            resids.append(res.resid())
            res_names.append(x)
            atom_names.append(list(atoms.atoms().extract_name()))
            atom_selection.append(list(atoms.atoms().extract_i_seq()))
            no_altloc.append(True)
      #
      chains_info[ch.id].resid = resids
      chains_info[ch.id].res_names = res_names
      chains_info[ch.id].atom_names = atom_names
      chains_info[ch.id].atom_selection = atom_selection
      chains_info[ch.id].no_altloc = no_altloc
  else:
    # build chains_info from selection
    l = ph.atoms().size()
    selection_test = flex.bool([False,] * l)
    chain_ids = sorted(selection_list)

    for sel_str in chain_ids:
      ph_sel = ph.select(cache(sel_str))
      model =  ph_sel.models()
      if model:
        # check that none-empty selection
        model = model[0]
        chains_info[sel_str] = Chains_info()
        chains_info[sel_str].chains_atom_number = ph_sel.atoms().size()
        coc = flex.vec3_double([ph_sel.atoms().extract_xyz().mean()])
        chains_info[sel_str].center_of_coordinates = coc
        for ch in model.chains():
          for res in ch.residue_groups():
            for atoms in res.atom_groups():
              x = atoms.resname
              if exclude_water and (x.lower() == 'hoh'): continue
              chains_info[sel_str].resid.append(res.resid())
              chains_info[sel_str].res_names.append(x)
              atom_list = list(atoms.atoms().extract_name())
              chains_info[sel_str].atom_names.append(atom_list)
              i_seq_list = list(atoms.atoms().extract_i_seq())
              chains_info[sel_str].atom_selection.append(i_seq_list)
              # test that there are not overlapping selection
              sel = flex.bool(l,flex.size_t(i_seq_list))
              if (selection_test & sel).count(True) > 0:
                raise Sorry('Overlapping NCS group selections!')
              else:
                selection_test |= sel
      else:
        raise Sorry('Empty NCS group selections!')
  return chains_info

def initialize_score_matrix(row, col,max_score=1000):
  """
  initialize a matrix in a dictionary form
  The matrix values are initialized according the the gap penalties.

  We initializing the zero row and column with descending score and the reset
  of the matrix with -10000 score

  Args:
    row (int): number of rows
    col (int): number of columns
    max_score (int): the score assigned for perfect alignment

  Returns:
    R (dict): the score matrix in a form of a dictionary of Score_record objects
  """
  R = {}
  for i in xrange(row + 1):
    score_row = max_score - i
    R[i,0] = Score_record(score=score_row,origin=(i-1,0))
    for j in xrange(1,col + 1):
      if i == 0:
        score_col = max_score - j
        R[i,j] = Score_record(score=score_col,origin=(0,j-1))
      else:
        R[i,j] = Score_record()
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

def is_same_transform(r1,t1,r2,t2):
  """
  Check if transform is the same by comparing rotations and the result of
  applying rotation and translation on
  a test vector

  Args:
    r1, r2: Rotation matrices
    t1, t2: Translation vectors

  Returns:
    (bool,bool) (is_the_same, is_transpose)
  """
  # Allowed deviation for values and angle
  eps=0.1
  angle_eps=5.0
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

def sort_by_dist(chains_info):
  """
  The process of grouping chains is not exhaustive, meaning that is does not
  consider every possible chains grouping, due to potential computation time
  such a search might take. To improve the the probability of getting
  preferred grouping, instead of searching using the order of the chain
  names, we arrange chains by distance.
  Starting from the first chain look for the closest chain to it, then
  the closest to the second chain, and so one.

  Note that this only improve but does not guarantee optimal grouping

  This function is used only when minimizing number of NCS operators

  Args:
    chains_info (obj): Chain_info object

  Return:
    ch_id_list (lst): ordered list of chain IDs
  """
  sorted_list = sorted(chains_info)
  if sorted_list:
    ch_id_list = [sorted_list.pop(0)]
    while sorted_list:
      coc1 = chains_info[ch_id_list[-1]].center_of_coordinates
      min_d = 10000000
      min_ch_id = ''
      min_indx = None
      for i,ch_id in enumerate(sorted_list):
        coc2 = chains_info[ch_id].center_of_coordinates
        d = coc1.max_distance(coc2)
        if d < min_d:
          min_d = d
          min_ch_id = ch_id
          min_indx = i
      if min_ch_id:
        ch_id_list.append(sorted_list.pop(min_indx))
    return ch_id_list
  else:
    return []


def  update_chain_ids_search_order(chains_info,sorted_ch,chains_in_copies,i):
  """
  This function is used only when minimizing number of chains in master NCS.
  It design to improve grouping by choosing masters that are near by.

  Replace chain ID i in the list sorted_ch with the closest chain i-1 in
  sorted_ch, ignoring chains that are already in NCS groups

  Args:
    chains_info (obj): Chain_info object
    sorted_ch (list): chain IDs list
    chains_in_copies (set): chain already in NCS groups
    i (int): position in the list


  Return:
    ch_id_list (list): ordered list of chain IDs

  23/09/2015 Oleg. Disabled because arbitrary tossing chains may cause
  selection problems in future execution. The test exercising it is in:
  cctbx_project/iotbx/regression/ncs/tst_ncs_input.py, exersice_23().
  In nutshell:
  given chain order in pdb file: ADEFCB and groups are ABC and DEF, where
  matches are A-D, B-E, C-F
  selections are going to extract atoms in the following order:
  ACB, DEF and they will be not atom-to-atom matched (the resulting matching
  is A-D, C-E, B-F)
  """
  return sorted_ch
  if i == 0:
    return sorted_ch
  else:
    ref = chains_info[sorted_ch[i-1]].center_of_coordinates
    min_d = 10000000
    min_ch_id = None
    # use only potential masters
    test_list = list(set(sorted_ch[i:]) - chains_in_copies)
    while test_list:
      ch_id = test_list.pop(0)
      c_of_c = chains_info[ch_id].center_of_coordinates
      d = ref.max_distance(c_of_c)
      if d < min_d:
        min_d = d
        min_ch_id = ch_id
    # flip the i position with the min_indx position
    min_indx = sorted_ch.index(min_ch_id)
    # alternatively, just place chain in front of current instead of flipping
    # Doesn't really help.
    # sorted_ch.remove(min_ch_id)
    # sorted_ch.insert(i, min_ch_id)
    sorted_ch[i], sorted_ch[min_indx] = sorted_ch[min_indx], sorted_ch[i]
    return sorted_ch


def my_get_rot_trans(
    ph,
    master_selection,
    copy_selection):
  m_atoms = ph.select(master_selection).atoms()
  c_atoms = ph.select(copy_selection).atoms()
  # print "selection master:", list(master_selection)
  # print "selection copy:", list(copy_selection)
  # Check that master and copy are identical
  assert len(m_atoms) == len(c_atoms), "%d, %d" % (len(m_atoms), len(c_atoms))
  # master
  other_sites = m_atoms.extract_xyz()
  # copy
  ref_sites = c_atoms.extract_xyz()
  if ref_sites.size() > 0:
    lsq_fit_obj = superpose.least_squares_fit(
        reference_sites = ref_sites,
        other_sites     = other_sites)
    r = lsq_fit_obj.r
    t = lsq_fit_obj.t
    rmsd = ref_sites.rms_difference(lsq_fit_obj.other_sites_best_fit())
    return r,t,rmsd
  else:
    assert 0, "strange thing..."

def get_rot_trans(ph,
                  master_selection,
                  copy_selection,
                  max_rmsd=0.02):
  """
  Get rotation and translation using superpose.

  This function is used only when phil parameters are provided. In this case
  we require the selection of NCS master and copies to be correct.
  Correct means:
    1) residue sequence in master and copies is exactly the same
    2) the number of atoms in master and copies is exactly the same

  One can get exact selection strings by ncs_object.show(verbose=True)

  Args:
    ph : pdb hierarchy input object
    master/copy_selection (str): master and copy selection strings
    max_rmsd (float): limit of rms difference between chains to be considered
      as copies

  Returns:
    r: rotation matrix
    t: translation vector
    rmsd (float): RMSD between master and copy
    msg (str): error messages
  """
  msg = ''
  r_zero = matrix.sqr([0]*9)
  t_zero = matrix.col([0,0,0])
  #
  if hasattr(ph,'hierarchy'):
    ph = ph.hierarchy
  if ph:
    cache = ph.atom_selection_cache().selection
    master_ncs_ph = ph.select(cache(master_selection))
    ncs_copy_ph = ph.select(cache(copy_selection))
    seq_m,res_ids_m  = get_residue_sequence(master_ncs_ph)
    seq_c,res_ids_c = get_residue_sequence(ncs_copy_ph)
    res_sel_m, res_sel_c, similarity = res_alignment(
      seq_a=seq_m,seq_b=seq_c,
      min_contig_length=0,min_percent=0)
    #
    m_atoms = master_ncs_ph.atoms()
    c_atoms = ncs_copy_ph.atoms()
    # Check that master and copy are identical
    if (similarity != 1) or (m_atoms.size() != c_atoms.size()) :
      return r_zero,t_zero,0,'Master and Copy selection do not exactly match'
    # master
    other_sites = m_atoms.extract_xyz()
    # copy
    ref_sites = c_atoms.extract_xyz()
    if ref_sites.size() > 0:
      lsq_fit_obj = superpose.least_squares_fit(
          reference_sites = ref_sites,
          other_sites     = other_sites)
      r = lsq_fit_obj.r
      t = lsq_fit_obj.t
      rmsd = ref_sites.rms_difference(lsq_fit_obj.other_sites_best_fit())
      if rmsd > max_rmsd:
        return r_zero,t_zero,0,msg
    else:
      return r_zero,t_zero,0,'No sites to compare.\n'
    return r,t,round(rmsd,4),msg
  else:
    return r_zero,t_zero,0,msg

def get_residue_sequence(ph):
  """
  Get a list of residues numbers and names from hierarchy "ph", excluding
  water molecules

  Args:
    ph (hierarchy): hierarchy of a single chain

  Returns:
    res_list_new (list of str): list of residues names
    resid_list_new (list of str): list of residues number
  """
  res_list_new = []
  resid_list_new = []
  for res_info in ph.atom_groups():
    x = res_info.resname
    if x.lower() != 'hoh':
      # Exclude water
      res_list_new.append(x)
      resid_list_new.append(res_info.parent().resseq)
  return res_list_new,resid_list_new
