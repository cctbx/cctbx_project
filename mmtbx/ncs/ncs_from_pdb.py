"""
  This code is a portion of simple_ncs_from_pdb.py (tct 2006-12-12)
  moved into cctbx (Aug-2014) to make the core functionality available
  without phenix dependencies

  The original file is at phenix\phenix\command_line\simple_ncs_from_pdb.py

  The original functionality and interface is maintained
  via simple_ncs_from_pdb.py
*****************************************************************************

Purpose:
  Figure out the ncs relations from the chains in a pdb object

Major assumption: residue numbers are consistent among chains

Approach:
  Use residue numbers to align the residue names, identify
  pairs of chains that can match. Choose groupings of chains that maximize the
  smallest number of matching residues between each member of a group and the
  first (reference) member of the group. Within a pair of chains, allow some
  segments to match and others not. Each pair of segments must have a
  length >= min_length and percent identity >=min_percent.  A pair of segments
  may not end in a mismatch. An overall pair of chains must have an rmsd
  of CA atoms of <= rmsd_max.

  If find_invariant_domain is specified then once all chains that can be matched
  with the above algorithm are identified, all remaining chains are matched,
  allowing the break-up of chains into invariant domains. The invariant
  domains each get a separate NCS group.

  If residue numbers are not the same for corresponding chains, but
  they are simply offset by a constant for each chain, this will be
  recognized and the chains will be aligned.
"""
from __future__ import division
from mmtbx.invariant_domain import find_domain
from cctbx.array_family import flex
from libtbx import adopt_init_args
from scitbx.math import superpose
from scitbx.math import matrix
from libtbx.utils import Sorry
from mmtbx.ncs.ncs import ncs
import string
import math
import sys


class ncs_from_pdb(object):
  """ find ncs relation in a biomolecules structure """

  def __init__(self,
               verbose=False,
               log=sys.stdout,
               njump=1,
               min_length=1,
               min_percent=0.95,
               suggested_ncs_groups=[],
               require_equal_start_match=True,
               maximize_size_of_groups=True,
               required_chains=[],
               min_fraction_domain=0.2,
               initial_rms=0.5,
               match_radius=2.0,
               similarity_threshold=0.75,
               min_contig_length=3,
               smooth_length=0,
               max_rmsd_domain=2,
               min_fraction_represented=0.10):
    """  Setup common variables

    When two chains match but starts at a different residue number, use
    offset_dict to indicate that difference. For example Chains A and B
    match, but chain A first residue number is 1 and B first residue number
    is 10, so offset_dict[AB] = 9

    if njump is not 1, then we are only going to take residues with
    numbers such that mod(resno,njump)=0
    so if njump=10 then we take residue 0,10,20,30,40 and put them
    in position 1,2,3,4,5

    Args:
      verbose (bool): print out during ncs processing
      log: indicate where to log output
      njump (int): residue sampling factor, reduces the number of residues used
      min_length (int): minimum number of matching residues in a segment for
        recursive call
      min_percent (float): min percent identity of matching residues
      suggested_ncs_groups (list):
      require_equal_start_match (bool): require that all matching
        segments start at the same relative residue number for all members
        of an NCS group, trimming the matching region as necessary. This is
        same, but not otherwise
      maximize_size_of_groups (bool): request that the scoring be set up to
        maximize the number of members in NCS groups
        (maximize_size_of_groups=True) or that scoring is set up to maximize
        the length of the matching segments in the NCS group
        (maximize_size_of_groups=False)
      required_chains (list):
      min_fraction_domain (float):domain must be this fraction of a chain
      initial_rms (float): Guess of RMS among chains
      match_radius (float): Keep atoms that are within match_radius of
        NCS-related atoms
      similarity_threshold (float): Threshold for similarity between segments
      min_contig_length (int): segments < min_contig_length rejected
      smooth_length (int): two segments separated by smooth_length or less get
        connected
      max_rmsd_domain (int): max rmsd of domains
      min_fraction_represented (float): Minimum fraction of residues
        represented by NCS to keep. If less...skip ncs entirely
    """
    self.offset_dict = {}
    self.total_residues = 0
    adopt_init_args(self, locals())

  def get_chain_list(self,hierarchy, ignore_chains):
    """
    get list of chains. Each chain list is a list of (blanks or residue names)
    a chain list has a slot for each residue in a chain, it can be occupied
    or not.  Two chains have the same parent chain if all residues that are
    present match exactly.

    each chain starts at a specific residue number, stored as
    starting_residue_numbers for that chain

    Args:
      hierarchy
      ignore_chains

    Returns:
      chains (list of lists): list of the residue sequence list (one letter
        characters) of each chain
      chain_ids (list of str): list of chains IDs
      starting_residue_numbers (list): list of chains first residue numbers
      offset_dict (dict): offset_dict[id1+id2] to start2 to match seq of chain1
      total_residues (int): Total number of residues in the hierarchy
    """
    chains = []
    chain_ids = []
    starting_residue_numbers = []
    model_number = 0
    total_residues = 0
    for model in hierarchy.models():
      model_number+=1
      if model_number>1: continue  # take first model only
      for chain in model.chains():
        if chain.id in ignore_chains: continue
        chain_id=chain.id
        if chain_id in chain_ids:
          msg = "NOTE: duplicate chain ID (ignoring): "
          if self.verbose:print >>self.log,msg,chain_id
          continue
        if self.verbose: print >>self.log,"Examining chain ",chain_id
        conformer_number = 0
        for conformer in chain.conformers():
          conformer_number+=1
          if conformer_number>1:continue  # take first conformer only
          highest_residue_number=-999999
          lowest_residue_number=  9999999
          number_of_good_residues = 0
          for residue in conformer.residues():  # get highest residue number
            resseq_int = residue.resseq_as_int()
            total_residues+=1
            if self.skip(resseq_int): continue  # based on njump
            if resseq_int>highest_residue_number:
               highest_residue_number=resseq_int
            if resseq_int<lowest_residue_number:
               lowest_residue_number=resseq_int
          if self.verbose:
            msg = "Residue range found: "
            print >>self.log,msg,lowest_residue_number,highest_residue_number
          seq = []
          res_dif = highest_residue_number-lowest_residue_number
          res_range = res_dif//self.njump+1
          for i in xrange(res_range):
             seq.append(None)
          residues_seen_already = []
          for residue in conformer.residues():
            resseq_int = residue.resseq_as_int()
            if self.skip(resseq_int): continue  # based on njump
            if resseq_int in residues_seen_already:continue
            residues_seen_already.append(resseq_int)
            one_char = self.one_char(residue.resname)
            if one_char: # self.skip anything that is not protein and not DNA/RNA
              found_ca_or_p=False  # require ' CA ' or ' P  ': see get_xyz_list
              for atom in residue.atoms():
                 if (atom.name == " CA " or atom.name == " P  "):
                   found_ca_or_p=True
                   break
              if found_ca_or_p:
                number_of_good_residues+=1
                #    position of this residue
                i = (resseq_int - lowest_residue_number)//self.njump
                seq[i]=one_char
                if self.verbose:
                  print >>self.log,"Residue",residue.resseq,one_char
        if number_of_good_residues >= self.min_length//self.njump:
          chains.append(seq)
          chain_ids.append(chain_id)
          starting_residue_numbers.append(lowest_residue_number)

    # Now if some chains are same as other chains but residue numbers are
    # offset, figure that out...
    offset_dict={}  # add offset_dict[id1+id2] to start2 to match seq of chain1
    zipped_data = zip(chains,chain_ids,starting_residue_numbers)
    l = len(chains)
    for i in xrange(l):
      chain1,id1,start1 = zipped_data[i]
      for j in xrange(i+1):
        chain2,id2,start2 = zipped_data[j]
        offset, offset_key=self.find_offset(chain1,id1,start1,chain2,id2,start2)
        if offset is None:
          offset = 0
          if self.verbose:print >>self.log,"FAILED OFFSET ",offset,id1,id2,id1+id2
        else:
          if self.verbose:print >>self.log,"OFFSET ",offset,id1,id2,id1+id2
        offset_dict[offset_key]=offset
    self.offset_dict = offset_dict
    self.total_residues = total_residues
    return chains,chain_ids,starting_residue_numbers,offset_dict,total_residues

  def find_offset(self,chain1,id1,start1,chain2,id2,start2):
    """
    return offset that makes chain1 same as chain2 (no mismatches allowed)
    list chain1:
    first residue of chain2 can match any one of residues of chain1.
    offset can be:
    start2+offset=start1-len(chain2)+1  to start2+offset=start1+len(chain1)

    Args:
      chain1 (list): Residue single character list
      id1 (str): chain1 ID
      start1 (int):  number of first residue in chain 1
      chain2 (list): Residue single character list
      id1 (str): chain1 ID
      start2 (int): number of first residue in chain 2

    Returns:
      best_offset (dict): relative chains offset
      offset_key  (str): sorted by chain names, dictionary key
    """
    temp = sorted([[id1,1],[id2,-1]],key=lambda x: x[0])
    reverse_offset = temp[0][1]
    offset_key = temp[0][0] + temp[1][0]
    best_overlap = 0
    best_offset=None
    for offset_a in xrange(-len(chain2)+1,len(chain1)+1):
     offset = start1 - start2 + offset_a * self.njump
     overlap = self.get_overlap(
       chain1=chain1,start1=start1,chain2=chain2,start2=start2,offset=offset)
     if overlap and (overlap > best_overlap):
       best_overlap = overlap
       best_offset = offset * reverse_offset
    return best_offset, offset_key

  def get_overlap(self,chain1,start1,chain2,start2,offset):
    """
    Args:
      chain1 (list): Residue single character list
      start1 (int):  number of first residue in chain 1
      chain2 (list): Residue single character list
      start2 (int): number of first residue in chain 2
      offset (int):

    Returns:
      n_match
    """
    offset1 = 0
    offset2 = offset
    n_match = 0
    mismatches = 0
    for i in xrange(len(chain1)):
      resno = i * self.njump + start1 + offset1
      res1 = chain1[i]
      #resno=j*self.njump+start2+offset2
      j = (resno - start2 - offset2)//self.njump
      if (j >= 0) and (j < len(chain2)):
        res2 = chain2[j]
      else:
        res2 = None
      if (res1 != None) and (res2 != None):
        if res1 == res2:
          n_match += 1
        else:
          mismatches += 1
    if 100.*float(n_match)/float(max(1,n_match+mismatches)) < self.min_percent:
      return None
    if n_match >= self.min_length//self.njump:
      return n_match
    else:
      return None

  def skip(self,resno):
    """  skip unless resno (residue number) falls on n * self.njump    """
    if self.njump <=1:
      return False
    if int(0.5+math.fmod(100*self.njump+resno,self.njump)) != 0:
      return True

  def one_char(self,residue_name):
    """
    return one char for DNA/RNA/Protein
    070508 catch +G or DG and call both of them G:
    2010-12-01 also Ar == A (ribo) call it A

    Args:
      residue_name (str): 3 letters residue code

    Returns:
      (str): One letter residue code
    """
    edited_residue_name=residue_name.lower().replace(" ","")
    if len(edited_residue_name)==2 and edited_residue_name[0] in ['d','+']:
      edited_residue_name=edited_residue_name[1]
    elif len(edited_residue_name)==2 and edited_residue_name[1] in ['r']:
      edited_residue_name=edited_residue_name[0] #2010-12-01
    elif edited_residue_name in ['ade','cyt','gua','thy','uri']:
       for one,three in zip(
           ['a','c','g','t','u'],
           ['ade','cyt','gua','thy','uri']):
         if edited_residue_name==three:
            edited_residue_name=one
            break
    aa_names= [
      'GLY','ALA','SER','VAL',
      'ILE','LEU','MET','CYS',
      'PHE','TYR','LYS','ARG',
      'TRP','HIS','GLU','ASP',
      'GLN','ASN','PRO','THR',
      'MSE', # treat MSE as MET for now
      'G','C','A','T','U']
    aa_letters=[
       'G','A','S','V',
       'I','L','M','C',
       'F','Y','K','R',
       'W','H','E','D',
       'Q','N','P','T',
       'M',                   # treat MSE as MET for now
       'G','C','A','T','U']

    for name,letter in zip(aa_names,aa_letters):
      if string.replace(string.lower(edited_residue_name)," ","") == \
         string.replace(string.lower(name)," ",""):
        return string.lower(letter)
    return None

  def find_groups(self,
                  hierarchy,
                  chains,
                  chain_ids,
                  starting_residue_numbers,
                  min_length=5,
                  min_percent=100.,
                  max_rmsd=5.,
                  max_rmsd_user=15,
                  called_by_self=False,
                  exact_match_only=False,
                  quick=True):
    """
    Args:
      hierarchy (object)
      chains (list of lists): list of the residue sequence list (one letter
        characters) of each chain
      chain_ids (list of str): list of chains IDs
      starting_residue_numbers (list): list of chains first residue numbers
      min_length (int): minimum number of matching residues in a segment for
        recursive call
      min_percent (float): min percent identity of matching residues
      max_rmsd  (float): max rmsd of domains
      max_rmsd_user (int): max rmsd of chains suggested by user
        (i.e., if called from phenix.refine with suggested ncs groups)
      called_by_self (bool):
      exact_match_only (bool):
      quick (bool):

    Returns:
      groups
      list_of_residue_range_list
    """
    match_list = []
    zipped_data = zip(chains,chain_ids,starting_residue_numbers)
    data_length = len(chains)
    for i in xrange(data_length):
      chain1,chain_id1,start1 = zipped_data[i]
      # if quick is set and all chains match to the first one we have enough
      if quick and self.have_enough(chain_ids,match_list): continue
      for j in xrange(i,data_length):
        chain2,chain_id2,start2 = zipped_data[j]
        if self.verbose:
           print >>self.log,"Comparing ",chain_id1,chain_id2
        matches,keep_range_list = self.matches_approx(
          chain1,start1,chain2,start2,chain_id1,chain_id2,
          min_length=min_length,min_percent=min_percent,
          exact_match_only=exact_match_only)
        if matches:
          # test for max_rmsd now and simply reject if too big
          group = [chain_id1,chain_id2]
          in_suggested_ncs_groups = self.are_in_suggested_ncs_groups(group)
          # reject if in >1 suggested NCS groups:
          if self.are_in_different_suggested_ncs_groups(group): continue
          residue_range=keep_range_list # list of ranges: [[1,3],[7-9]]
          residue_range_list=[residue_range,residue_range]
          # Look for Residue range used... it is 10-180
          rmsd_list,r_list,trans_list,center_list,residues_in_common_list=\
          self.get_rmsd(group,hierarchy,chains,chain_ids,
             starting_residue_numbers,residue_range_list)
          if (not rmsd_list) or (len(rmsd_list)<2):
            pass
          elif (not in_suggested_ncs_groups) and (rmsd_list[1]>max_rmsd):
            if self.verbose:
              msg = "Alignment ",group," rejected with rmsd =",rmsd_list[1]
              print >>self.log,msg
          elif in_suggested_ncs_groups and (rmsd_list[1]>max_rmsd_user) and \
                  (rmsd_list[1]>max_rmsd):
            print >>self.log,"NOTE: User-suggested alignment ",group, \
                " rejected with rmsd =",rmsd_list[1]
            msg = "Set max_rmsd_user="+str(int(rmsd_list[1]+0.9999))
            print >>self.log, msg + " to keep it"
          else:
            match_list.append([matches,chain_id1,chain_id2,residue_range])
    if self.verbose: print >>self.log,"MATCH LIST: ",match_list

    # now figure out what are in clusters together, and which is the best
    #  reference for each group ...
    remaining_ids=chain_ids
    still_finding_groups=True
    groups = []
    list_of_residue_range_list = []

    while still_finding_groups:
       best_id=None
       best_score=-99999
       best_group=None
       best_residue_range_list=None
       for trial_id in remaining_ids: # use this as seed to pull out all others
         score,group,residue_range_list=self.get_score_and_group(trial_id,
             remaining_ids,match_list)
         if score>best_score:
            best_score=score
            best_id=trial_id
            best_residue_range_list=residue_range_list
            best_group=group
       if best_score>0:
         if self.verbose:
            print >>self.log,"Best new group:",best_group,best_score,\
              best_residue_range_list
         groups.append(best_group)
         list_of_residue_range_list.append(best_residue_range_list)
         new_remaining_ids = []
         for id in remaining_ids:
           if not id in best_group:
             new_remaining_ids.append(id)
         remaining_ids=new_remaining_ids
         if self.verbose:
           print >>self.log,"Current remaining ids",remaining_ids
       else:
         if self.verbose:
            print >>self.log,"No new groups"
         still_finding_groups=False
    if not called_by_self and self.find_invariant_domains and remaining_ids:
       if self.verbose:
          print >>self.log,"Current groups: ",groups,list_of_residue_range_list
       # find groups from remaining_ids...this time taking ANY rmsd:
       remaining_chains = []
       remaining_chain_ids = []
       remaining_starting_residue_numbers = []
       for chain,chain_id,starting_residue_number in zip (
          chains,chain_ids,starting_residue_numbers):
         if chain_id in remaining_ids:
          remaining_chains.append(chain)
          remaining_chain_ids.append(chain_id)
          remaining_starting_residue_numbers.append(starting_residue_number)
       additional_groups,additional_list_of_residue_range_list=\
         self.find_groups(
           hierarchy,remaining_chains,remaining_chain_ids,
           remaining_starting_residue_numbers,
           min_length=self.min_length,min_percent=self.min_percent,
           max_rmsd=999999.,max_rmsd_user=15,called_by_self=True)

       if additional_groups:
         invariant_groups,invariant_list_of_residue_range_list=\
            self.find_invariant_groups(hierarchy,
             additional_groups,additional_list_of_residue_range_list,
             remaining_chains,remaining_chain_ids,
             remaining_starting_residue_numbers)
         if invariant_groups:
           print >>self.log,"Invariant groups found:",\
              invariant_groups,invariant_list_of_residue_range_list
           return groups+invariant_groups, \
             list_of_residue_range_list+invariant_list_of_residue_range_list

    return groups,list_of_residue_range_list

  def find_invariant_groups(self,hierarchy,
             additional_groups,additional_list_of_residue_range_list,
             remaining_chains,remaining_chain_ids,
             remaining_starting_residue_numbers):
    """
    for each group in additional groups, see if there is some way to break
    it down into invariant sub-groups. Return list of these and new
    residue ranges...

    Returns:
      invariant_groups,invariant_list_of_residue_range_list
    """
    invariant_groups=[]
    invariant_list_of_residue_range_list=[]
    if not additional_groups or not additional_list_of_residue_range_list:
      return invariant_groups,invariant_list_of_residue_range_list

    for group,residue_range_list in zip(
      additional_groups,additional_list_of_residue_range_list):
      if self.verbose:
        print >>self.log,"Looking for invariant domains for ...:",group,\
          self.add_offsets(residue_range_list,group)
      # we are going to require that our domains be the same for all members
      # of this group....Find the domains in pairs with first chain in
      #  list...and choose the set of domains that best fits
      #  ALL comparisons to this chain
      id1 = group[0]
      chain1,start1=self.get_chain(
        id1,remaining_chains,remaining_chain_ids,
        remaining_starting_residue_numbers)
      residue_range1=residue_range_list[0]
      best_rmsd=None
      best_score=None
      best_matches=[]
      best_domain_list=[]
      best_domains=[]
      best_overall_rmsd_list=[]
      best_number_of_residues_used=[]
      for id2,residue_range2 in zip(group,residue_range_list):
          if id1==id2: continue
          if self.verbose: print >>self.log,"Using ",id1,id2," to find domains..."
          chain2,start2=self.get_chain(id2,remaining_chains,
             remaining_chain_ids,remaining_starting_residue_numbers)
          if not chain1 or not chain2:
             raise Sorry("Chain id "+str(id1)+" "+str(id2)+
               " are not both present in remaining chain_ids: "+
               str(remaining_chain_ids))

          residue_list=self.get_residues_in_common(
            chain1,start1,chain2,start2,residue_range1,residue_range2,
            id1,id2)
          if not residue_list: continue  # 021107
          if self.verbose:
            print >>self.log,"residues in common",residue_list
          # get sites from first and second...
          offset1 = 0
          offset2=self.offset_dict[''.join(sorted([id1,id2]))]
          center1,xyz_list1=\
            self.get_xyz_list(id1,hierarchy,residue_list,offset1)
          center2,xyz_list2=\
            self.get_xyz_list(id2,hierarchy,residue_list,offset2)
          if self.verbose: print >>self.log,"Centers: ",center1,center2
          min_size=int(
             self.min_fraction_domain*float(len(chain1))/float(self.njump)
             )  # 2011-03-30
          if xyz_list1 is None: continue  # 021107
          if min_size>len(xyz_list1):
             min_size=len(xyz_list1)

          if len(xyz_list1)<4: continue

          domains = find_domain(
            xyz_list1,
            xyz_list2,
            initial_rms=self.initial_rms,
            match_radius=self.match_radius,
            overlap_thres=self.similarity_threshold,
            minimum_size=min_size)
          # make sure each residue is in only 1 domain...
          domain_list=[]
          used_list=[]
          for match in domains.matches:
            domain_residues=[]
            keep_list=[]
            iselect=match[0].as_1d()
            for i in iselect:
              if not i in used_list:
                keep_list.append(i)
                # save list of actual residue ranges
                domain_residues.append(residue_list[i])
            used_list+=keep_list
            domain_list.append(self.remove_small(self.connect_range(
             self.list_to_range(domain_residues)) ))
          # sum up and remove any domains that are too small or too large rmsd:
          if self.verbose:print >>self.log,"\nDomains in ",id1,id2
          final_domain_list=[]
          final_matches=[]
          for domain_residues,match in zip(domain_list,domains.matches):
            if not domain_residues or not match:
               if self.verbose:print >>self.log,"Domain completely removed..."
               continue
            if self.verbose:
              print >>self.log,"Domain: ",domain_residues
              print >>self.log,"RMSD: ",match[3]
              print >>self.log,"N: ",match[4]
            if match[3]>self.max_rmsd_domain:
               if self.verbose: print >>self.log,"RMSD too large...rejecting this domain"
               continue
            final_domain_list.append(domain_residues)
            final_matches.append(match)
          if not final_matches or not final_domain_list:
            if self.verbose:print >>self.log,"No domains from this pair"
            continue

          if self.verbose:
           print >>self.log,"Evaluating rmsd for all chain pairs using these final matches"
          overall_rmsd_list=[]
          # get rmsd of the whole group using residues in final_domain_list
          number_of_residues_used_list=[]
          for domain_residues,match in zip(final_domain_list,final_matches):
            if self.verbose: print >>self.log,"domain...",domain_residues
            list_of_domain_residues=[] # testing this definition of the domain..
            for id in group: list_of_domain_residues.append(domain_residues)
            rmsd_list,r_list,trans_list,center_list,residues_in_common_list=\
                self.get_rmsd(group,hierarchy,
                remaining_chains,remaining_chain_ids,
                remaining_starting_residue_numbers,list_of_domain_residues)
            if self.verbose:
              print >>self.log,"rmsd_list for ",domain_residues,':',rmsd_list

            if (not rmsd_list) or (len(rmsd_list)<2):
               overall_rmsd_list=[]
               break
            overall_rmsd = 0.
            mean_residues_used = 0
            for rmsd,res in zip(rmsd_list[1:],residues_in_common_list[1:]):
              overall_rmsd+=rmsd
              mean_residues_used+=res
            overall_rmsd=overall_rmsd/float(len(rmsd_list)-1)  # skip self-rmsd
            mean_residues_used=mean_residues_used/float(len(rmsd_list)-1)
            overall_rmsd_list.append(overall_rmsd)
            number_of_residues_used_list.append(mean_residues_used)
          test_rmsd = 0.
          score=None
          if overall_rmsd_list:
            test_mean_length = 0.
            for rmsd,res in zip(overall_rmsd_list,number_of_residues_used_list):
              test_rmsd+=rmsd
              test_mean_length+=res
            test_rmsd=test_rmsd/float(len(overall_rmsd_list))
            test_mean_length=test_mean_length/float(len(overall_rmsd_list))
            score = max(0.,(self.max_rmsd_domain-test_rmsd))*test_mean_length
          if self.verbose:
           print >>self.log,"rmsd for this match: ",test_rmsd,score
          if score is not None and (best_score is None or score>best_score):
            best_score=score
            best_rmsd=test_rmsd
            best_matches=final_matches
            best_domain_list=final_domain_list
            best_domains=domains
            best_number_of_residues_used_list=number_of_residues_used_list
            best_overall_rmsd_list=overall_rmsd_list
      if best_domains is not None and best_score is not None:
        if self.verbose:
          print >>self.log,"Best rmsd for group ",group,": ", best_rmsd,score
          print >>self.log,"Rmsd and residues used by domain: "
          for rmsd,res in zip(best_overall_rmsd_list,
              best_number_of_residues_used_list): print >>self.log,rmsd,res
        if self.verbose: best_domains.show(out=self.log)
        for domain_list in best_domain_list:
          invariant_groups.append(group)
          list_of_residue_ranges_for_this_group=[]
          for i in group:
            list_of_residue_ranges_for_this_group.append(domain_list)
          invariant_list_of_residue_range_list.append(
            list_of_residue_ranges_for_this_group)

    return invariant_groups,invariant_list_of_residue_range_list

  def remove_used_chains(self,chains,chain_ids,
                         starting_residue_numbers,used_ids):
    [new_chains,new_chain_ids,new_starting_residue_numbers]=[[],[],[]]
    for chain,id,start in zip(chains,chain_ids,starting_residue_numbers):
      if not id in used_ids:
        new_chains.append(chain)
        new_chain_ids.append(id)
        new_starting_residue_numbers.append(start)
    return [new_chains,new_chain_ids,new_starting_residue_numbers]

  def add_ids(self,groups,used_ids):
    for group in groups:
      for id in group:
        if not id in used_ids: used_ids.append(id)
    return used_ids

  def too_few_residues_represented(self,
                                   ncs_object=None,
                                   total_residues=None):
    """
    count up how many residues are represented in the NCS and return True
    if not a fraction self.min_fraction_represented
    """
    total_represented = 0
    for ncs_group in ncs_object.ncs_groups():
      for common in ncs_group.residues_in_common_list():
        total_represented+=common

    if total_residues is not None and total_residues>0:
      fraction=float(total_represented)/float(total_residues)
      if fraction < self.min_fraction_represented:
         return True # too few
      else:
         return False # not too few
    return False # do not know

  def list_to_range(self,list):
    """ [1,2,3,6,7] -> [[1,3],[6,7]]
    print >>self.log,"LIST",list  """
    range_list=[]
    work_list=[]
    i_prev=None
    for i in list:
      if i_prev==None or  i==i_prev+1:
        i_prev=i
        work_list.append(i)
      else:
        if work_list:
          range_list.append(work_list[0:1]+work_list[-1:])
          work_list=[i]
          i_prev=i
    if work_list:
          range_list.append(work_list[0:1]+work_list[-1:])
          work_list=[]
    #print >>self.log,"range_list: ",range_list
    return range_list

  def get_chain(self,id1,remaining_chains,remaining_chain_ids,
                remaining_starting_res_numbers):
    for chain,id,start in zip(remaining_chains,remaining_chain_ids,
        remaining_starting_res_numbers):
      if id1==id:
         return chain,start
    return [],0

  def remove_small (self,range_list):
    """
    remove all ranges that are smaller than self.min_contig_length//self.njump
    """
    contig_length = self.min_contig_length//self.njump
    new_range_list=[]
    for range in range_list:
       if range[1]-range[0]+1 < contig_length:
          if self.verbose: print >>self.log,"Rejecting short range ",range
       else:
          new_range_list.append(range)
    return new_range_list

  def connect_range(self,range_list):
    """ If two ranges are within self.smooth_length, connect them """
    smooth = self.smooth_length * self.njump
    new_range_list=[]
    work_range=[]
    for range in range_list:
      if not work_range:
        work_range=range
      elif range[0] <= work_range[1]+smooth+1:
        work_range[1]=range[1]
      else:
        new_range_list.append(work_range)
        work_range=range
    if work_range: new_range_list.append(work_range)
    return new_range_list

  def add_offsets(self,residue_range_list,group):
    """
    add offsets to all residue numbers in residue_range_list using
    information in self.offset_dict and group:
    """
    new_residue_range_list=[]
    id1=group[0]
    if self.verbose:print >>self.log,"OFFSETS:"
    for id2,residue_range in zip(group,residue_range_list):
      offset=self.offset_dict[''.join(sorted([id1,id2]))]
      if self.verbose: print >>self.log,residue_range,id1,id2,offset
      new_residue_range=[]
      for pair in residue_range:
        new_residue_range.append([pair[0]-offset,pair[1]-offset])
      new_residue_range_list.append(new_residue_range)
    if self.verbose:print >>self.log,"OLD, NEW: ",residue_range_list,new_residue_range_list
    return new_residue_range_list

  def subtract_offsets(self,residue_range_list,group):
    """
    remove offsets from all residue numbers in residue_range_list using
    information in self.offset_dict and group:
    """
    new_residue_range_list=[]
    id1=group[0]
    if self.verbose:print >>self.log,"OFFSETS:"
    for id2,residue_range in zip(group,residue_range_list):
      offset=self.offset_dict[''.join(sorted([id1,id2]))]
      if self.verbose: print >>self.log,residue_range,id1,id2,offset
      new_residue_range=[]
      for pair in residue_range:
        new_residue_range.append([pair[0]+offset,pair[1]+offset])
      new_residue_range_list.append(new_residue_range)
    if self.verbose:print >>self.log,"OLD, NEW: ",residue_range_list,new_residue_range_list
    return new_residue_range_list

  def get_suggested_group_list(self,suggested_ncs_groups,chain_ids):
    """
    121706: make sure that all suggested ids are actually in chain_ids!
    interpret suggested ncs groups and make a list of them in our format
    """
    suggested_group_list=[]
    for i_seq, group in enumerate(suggested_ncs_groups):
      if(group.reference is not None):
        selections=[group.reference]+group.selection
        group=[]
        range_list=[]
        first=True
        for sel in selections:
          selection_object=self.atom_selection(sel)
          if selection_object is None:
            print >>self.log,"No matching residues in suggested NCS groups"
            return None
          id,ranges=self.get_id_and_ranges(selection_object,chain_ids)
          if id is None or ranges is None:
            print >>self.log,"No matching residues in suggested NCS groups"
            return None
          if not id in chain_ids:
            raise Sorry("\nThe suggested chain "+str(id)+" was not found "+
                        " in the input file...it could \nbe the chain is too "+
                        "short (try decreasing min_length) \nor that it is not "+
                        "recognized as a chain. Note that ligands\n"+
                        " and solvent are not recognized as chains for NCS\n")
          group.append(id)
          range_list.append(ranges)
        if group and range_list:
          range_list_without_offsets=self.subtract_offsets(range_list,group)
          suggested_group_list+=[[group,range_list_without_offsets]]
          # allow list of lists
    return suggested_group_list

  def get_id_and_ranges(self,selection_object,chain_ids):
    """
    we want to return the chain ID and list of residue ranges referred
    to by this selection object. If there are >1 ID: raise an exception as
      we cannot be sure that it is handled correctly.
    assume that all residue numbers go up sequentially!
    In this routine there are no offsets.
    """
    id_list=[]
    residue_range_list=[]
    hierarchy = self.hierarchy
    selection = selection_object
    assert selection.size() == self.pdb_inp.atoms().size()
    for atom,is_selected in zip(self.pdb_inp.atoms(), selection):
      atom.tmp = int(is_selected)
    for model in hierarchy.models():
      first_in_model = True
      for chain in model.chains():
        first_in_chain = True
        conformer_number = 0
        for conformer in chain.conformers():
          conformer_number+=1
          if conformer_number>1: continue # take only first conformer
          first_in_conformer = True
          last_residue=None
          for residue in conformer.residues():
            if (self.is_selected(residue)):
              if (first_in_model):
                first_in_model = False
              if (first_in_chain):
                first_in_chain = False
              if (first_in_conformer):
                first_in_conformer = False

              if not chain.id in id_list: id_list.append(chain.id)
              resno=residue.resseq_as_int()
              if last_residue is None or resno != last_residue+1: # new segment
                if last_residue is not None:  # finished with a segment
                  residue_range_list.append([first_residue,last_residue])
                first_residue=resno
                last_residue=resno
              else:  # just increment last_residue
                last_residue=resno
          if last_residue is not None:
            residue_range_list.append([first_residue,last_residue])

    if len(id_list)==0 or not residue_range_list:
      return None,None

    if len(id_list)>1:
      raise Sorry("This version of simple_ncs_from_pdb requires any input NCS "+
                  "\nspecifications to have only one chain in each reference or selection"+
                  "\n('chain A or chain B' is not allowed)")

    return id_list[0],residue_range_list

  def is_selected(self,residue):
    for atom in residue.atoms():
      if (atom.tmp): return True
    return False

  def atom_selection(self, sel_str):
    try:
      selection=self.all_chain_proxies.selection(string = sel_str)
    except KeyboardInterrupt: raise
    except Exception:
      selection=None
    return selection

  def get_suggested_groups(self,suggested_ncs_groups,chain_ids):
    """ expect a phil object or text... """
    groups = []
    group = []
    suggested_group_list = []
    if not suggested_ncs_groups: return groups,suggested_group_list
    if not type(suggested_ncs_groups)==type("ABCD"):
      suggested_group_list=self.get_suggested_group_list(suggested_ncs_groups,
           chain_ids)
      return groups,suggested_group_list
    else:
     input_text=suggested_ncs_groups
     # expect text: "[ABC][DE]"
     # return a list of lists: [["A","B","C"],["D","E"]]
     for char in input_text:
      if char=="[":  # start a group
        if group:
          msg = 'Format for suggested_ncs_groups is "[ABCD][EFG]" or "ABCD"'
          raise Sorry(msg)
      elif char=="]":  # end a group
        if not group:
          raise Sorry(
            'Format for suggested_ncs_groups is "[ABCD][EFG]" or "ABCD"')
        group,groups=self.finish_suggested_group(group,groups)
      elif char in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
        if not char in chain_ids:
          raise Sorry ("The chain ID '"+str(char)+
         "' is not among the chain ids in this PDB file. \nAvailable chains: "+
         str(chain_ids))
        group.append(char)
     if group:
        group,groups=self.finish_suggested_group(group,groups)
     # check that each chain is used no more than once
     list_of_used_chains = []
     for group in groups:
       for id in group:
         if id in list_of_used_chains:
            raise Sorry(
            "Each chain can be used only once in suggested_ncs_groups"+\
            '\nChain '+str(id)+' is listed more than once in '+str(input_text))
         list_of_used_chains.append(id)
     return groups,suggested_group_list

  def finish_suggested_group(self,group,groups):
    if len(group)<2:
      raise Sorry('The group '+str(group)+
                  ' seems to have fewer than 2 members'+
                  '\nFormat for suggested_ncs_groups is "[ABCD][EFG]" or "ABCD"')
    groups.append(group)
    group=[]
    return group,groups

  def get_score_and_group(self,trial_id,remaining_ids,match_list):
    """
    score as the MINIMUM number of matches of anything to the trial_id
    if all the members of a user-specified group are present, add 10000
      to score
    return residue ranges that are not necessarily the same for all members
    of a group...
    090708: allow user to increase score for large groups with
    self.maximize_size_of_groups

    021107:
    Solution: need to use only those with id1==trial_id to keep numbering
    consistent. Then we also need matches to include both directions.

    021207 this didn't do it completely...we want to return offsets
    based not on trial_id but on the firstid in the list...

    the real problem is that here we use a seed (C) to pull out the
      set of chains (A C E G) and the offset matches C here. Then later
      we assume that the offset matched A. One solution might be
    to always use the quick option...using each set of chains
    identified with sequence. HERE JUST this is the best plan?

    021207 better: make sure that trial_id is first in list!
    090808 allow calling program to require a particular set of
    chains to always be present
    """
    score=None
    group = []
    residue_range_list = []

    for [matches,id1,id2,residue_range] in match_list:
      if id1==trial_id and id2==trial_id:
             group.append(trial_id)
             residue_range_list.append(residue_range)
    for [matches,id1,id2,residue_range] in match_list:
          if id1==trial_id and id2 in remaining_ids and not id2 in group:
             group.append(id2)
             residue_range_list.append(residue_range)
             if score is None:
               score=matches
             elif score > matches:
               score=matches  # decrease to lowest number of matches
    if self.require_equal_start_match:
      residue_range_list= \
         self.require_matches_to_start_at_same_place(group,residue_range_list)

    if self.maximize_size_of_groups and score is not None:
       score=score*math.sqrt(float(len(group)))

    if self.required_chains:
      found=False
      for id in group:
        if id in self.required_chains:
          found=True
          break
      if not found:
        score=None
        group = []
        residue_range_list = []

    if len(group)<2:
       score=None
       group = []
       residue_range_list = []
    if score and self.a_suggested_ncs_group_is_in_test_group(group):
      score+=10000

    return score,group,residue_range_list

  def a_suggested_ncs_group_is_in_test_group(self,test_group):
    """
    Returns:
      (bool): True -> all of a suggested_ncs_group are contained in test_group
    """
    for group in self.suggested_ncs_groups:
      are_in=True
      for id in group:
        if not id in test_group:
          are_in=False
          break
      if are_in:
        return True
    return False

  def require_matches_to_start_at_same_place(self,group,residue_range_list):
    """
    Returns:
      residue_range_list
    """
    new_residue_range_list=[]
    # find last start for each range
    if len(residue_range_list)<1: return 0,[]
    number_of_ranges=len(residue_range_list[0])
    last_start_list=number_of_ranges*[None]
    for rr,g in zip(residue_range_list,group):
      if number_of_ranges is None:
         number_of_ranges = len(rr)
      #print "group, range: ",g,rr
      new_start_list=[]
      rr_mod=rr[:len(last_start_list)]+(len(last_start_list)-len(rr))*[None]
      for r,s in zip(rr_mod,last_start_list):
        if r is not None and (s is None or r[0] > s):
           s = r[0]
        new_start_list.append(s)
      last_start_list=new_start_list
    #print "start_list; ",last_start_list
    # now adjust all the starting to be last_start_list
    for rr,g in zip(residue_range_list,group):
      new_rr=[]
      rr_mod=rr[:len(last_start_list)]+(len(last_start_list)-len(rr))*[None]
      for r,s in zip(rr_mod,last_start_list):
        if r is not None and s is not None and r[0] != s:
           r[0] = s  # modifies r inside rr inside residue_range_list
    return residue_range_list

  def get_rmsd(self,group,hierarchy,chains,chain_ids,starting_residue_numbers,
               residue_range_list):
    """
    get rmsd of all from the reference one, using only residues in
    residue_range_list for each chain
    residue_range_list is [residue_range1,residue_range2...]
    residue_range_1 is list of residue ranges: [[1,3],[7-9]...]

    NOTE: residue ranges are absolute
    """
    rmsd_list = []
    r_list = []
    t_list = []
    center_list = []
    residues_in_common_list = []
    if self.verbose:
      print >>self.log,"GET_RMSD: residue_range_list: ",residue_range_list
    if len(group)<2:
      return rmsd_list,r_list,t_list,center_list,residues_in_common_list
    chain_id1=group[0]
    residue_range1=residue_range_list[0]
    chain1,start1=self.find_chain_and_start(
       chain_id1,chains,chain_ids,starting_residue_numbers)
    if not chain1:
      return rmsd_list,r_list,t_list,center_list,residues_in_common_list
    for chain_id2,residue_range2 in zip(group,residue_range_list):
      chain2,start2 = self.find_chain_and_start(
        chain_id2,chains,chain_ids,starting_residue_numbers)
      if not chain2:
       return rmsd_list,r_list,t_list,center_list,residues_in_common_list
      msg = "getting CA values for ",chain_id1," and ",chain_id2
      if self.verbose: print >>self.log,msg
      # pull out CA coords from chain1 and chain2 for all residues where both are
      #   present (as reflected in chains[])
      # residue range1 and 2 are absolute
      residue_list = self.get_residues_in_common(
        chain1,start1,chain2,start2,residue_range1,
        residue_range2,chain_id1,chain_id2)
      if self.verbose: print >>self.log,"Residues in common: ",residue_list
      offset1 = 0
      offset2 = self.offset_dict[''.join(sorted([chain_id1,chain_id2]))]
      center1,xyz_list1=\
         self.get_xyz_list(chain_id1,hierarchy,residue_list,offset1)
      center2,xyz_list2=\
         self.get_xyz_list(chain_id2,hierarchy,residue_list,offset2)
      if self.verbose:print >>self.log,"Center: ",center1,center2
      if not xyz_list1 or len(xyz_list1)==0 or \
         not xyz_list2 or len(xyz_list2)==0 or \
         not len(xyz_list1)==len(xyz_list2):
        msg = "Note: No overlap of segments in getting NCS rmsd"
        if self.verbose: print >>self.log,msg
        return [],[],[],[],[],

      else:  # all ok...
        lsq_fit = superpose.least_squares_fit(
         reference_sites=xyz_list1,
         other_sites=xyz_list2)
        rmsd = xyz_list1.rms_difference(lsq_fit.other_sites_best_fit())
        rmsd_list.append(rmsd)
        residues_in_common_list.append(len(xyz_list1))
        r_list.append(lsq_fit.r)
        t_list.append(lsq_fit.t)
        center_list.append(center2)

        if self.verbose:
          xyz=flex.vec3_double()
          xyz.append(lsq_fit.r * matrix.col(xyz_list1[0]) + lsq_fit.t)
          print >>self.log,"Getting atom 1 of chain ",chain_id2,\
               " from atom 1 of chain ",chain_id1
          print >>self.log,"XYZ of atom 1 chain ",chain_id1,\
                  " after application of RT: ",xyz[0]
          print >>self.log,"xyz of atom 1 chain ",chain_id2," : ",xyz_list2[0]

    return rmsd_list,r_list,t_list,center_list,residues_in_common_list

  def get_residues_in_common(self,
                             chain1,start1,
                             chain2,start2,
                             residue_range1,
                             residue_range2,
                             id1,id2):
    """
    NOTE: residue ranges are absolute

    Returns:
      residue_list (list): list of residue numbers shared by chain1 and chain2
      in range defined by intersection of residue_range1 and residue_range2
      residue_range_1 is list of residue ranges: [[1,3],[7-9]...]
    """
    id = ''.join(sorted([id1,id2]))
    if id in self.offset_dict.keys():
      offset = self.offset_dict[id]
    else:
      print >>self.log,"NO ID : ",id
      offset = 0
    start2_use=start2+offset

    if self.verbose:
       print >>self.log,"Getting residues in common"
       print >>self.log,"Chain1, start1: ",chain1,start1
       print >>self.log,"Chain2, start2: ",chain2,start2,start2_use
       print >>self.log,"Residue_range1",residue_range1
       print >>self.log,"Residue_range2",residue_range2
    residue_list = []
    if not residue_range1 or not residue_range2:
      return residue_list
    all_residues1 = self.list_all_integers_in_range(residue_range1)
    all_residues2 = self.list_all_integers_in_range(residue_range2)
    all_to_use = []
    for res in all_residues1:
       if res in all_residues2: all_to_use.append(res)
    for res in all_to_use:
      i1 = (res-start1)//self.njump
      i2 = (res-start2_use)//self.njump
      # Matching
      if (i1 >= 0) and (i1 < len(chain1)) and (i2 >= 0) and (i2 < len(chain2)) \
              and chain1[i1] and chain2[i2] and chain1[i1]==chain2[i2]:
        residue_list.append(res)
    return residue_list

  def list_all_integers_in_range(self,residue_range_list):
    # NOTE: residue ranges are absolute
    all_int = []
    if not residue_range_list or len(residue_range_list)<1:return
    for range in residue_range_list:
       start=range[0]
       end=range[1]
       for i in xrange(start,end+1):
         if self.skip(i):continue
         all_int.append(i)
    return all_int

  def get_xyz_list(self,chain_id,hierarchy,residue_list,offset):
    """
    """
    if self.verbose: print >>self.log,"Getting xyz for ",chain_id,residue_list
    # NOTE: we are going to offset all the residue numbers in this chain
    # by offset...

    residue_found_list = []
    xyz_list=flex.vec3_double()
    model_number = 0
    center=flex.vec3_double()
    center.append( (0.,0.,0.))
    chain_ids_seen_already = []
    for model in hierarchy.models():
      model_number+=1
      # take first model only
      if model_number > 1: continue
      for chain in model.chains():
        if chain.id in chain_ids_seen_already: continue
        chain_ids_seen_already.append(chain.id)
        if chain.id != chain_id: continue # we only want chain_id...
        conformer_number = 0
        for conformer in chain.conformers():
          conformer_number+=1
          if conformer_number>1:continue  # take first conformer only
          residues_seen_already = []
          for residue in conformer.residues():
            resseq_int = residue.resseq_as_int()
            if self.skip(resseq_int):
              continue  # based on njump
            if resseq_int in residues_seen_already:continue
            residues_seen_already.append(resseq_int)
            if (resseq_int+offset) in residue_list:  # take CA from this
              found_ca_or_p=False
              for atom in residue.atoms():
                if (atom.name == " CA " or atom.name == " P  "):
                  xyz = atom.xyz
                  found_ca_or_p=True
                  break
              if found_ca_or_p:  # require CA: see get_chain_list
                xyz_list.append(xyz)
                residue_found_list.append(resseq_int+offset)
                center+=xyz
              else:
                if self.verbose:
                  print >>self.log,\
                  "NOTE: Missing CA/P for residue "+residue.resseq+ \
                  " of chain "+str(chain_id)
    if self.njump==1:   #  check only valid if njump==1
      for res in residue_list:
        if not res in residue_found_list:
          print >>self.log,"Residues expected: ",residue_list
          print >>self.log,"Residues found   : ",residue_found_list
          raise Sorry("Did not find expected residues in get_xyz_list"+
           "\nIt is possible that your PDB file has segments with one "+
           "chainID that are \ninterspersed with segments with another chainID")
    if len(xyz_list):
      center=center* (1./float(len(xyz_list)))
      return center[0],xyz_list
    else:
      return None,None

  def are_in_different_suggested_ncs_groups(self,test_group):
    """
    find a group that contains a member of test_group. If any other member
    is in another group, return True

    Returns:
      (bool): True  when members of test_group is contained
        in > 1 suggested ncs group
    """
    for group1 in self.suggested_ncs_groups:
     for id in test_group:
       # group1 contains an id of test_group
       if id in group1:
         for group2 in self.suggested_ncs_groups:
           if group1 == group2: continue
           for id in test_group:
             if id in group2:
               return True
    return False

  def find_chain_and_start(self,chain_id1,chains,chain_ids,
                           starting_residue_numbers):
    """
    get chain and starting residue number for chain with id chain_id
    """
    for chain,chain_id,start in zip(chains,chain_ids,starting_residue_numbers):
      if chain_id == chain_id1:
         return chain,start
    return None,None

  def are_in_suggested_ncs_groups(self,test_group):
    """
    Returns:
      (bool): True  when  all of test_group is contained in a suggested group
    """
    for group in self.suggested_ncs_groups:
      are_in=True
      for id in test_group:
        if not id in group:
          are_in=False
          break
      if are_in:
        return True  # all of test_group is contained in a suggested group
    return False

  def matches_approx(self,
                     chain1,start1,
                     chain2,start2,
                     id1,id2,
                     min_length=5,
                     min_percent=95,
                     exact_match_only=False):
    """
    return number of residues matching total, and a
    list of residue ranges over which chain1 and chain2 match
    within the tolerance min_percent and each with minimum length
    of min_length

    NOTE: residue range is in terms of position in chain (relative to start,
    and if njump>1 then the residue numbers are in  units of njump

    Returns:
      matching (int): number of matching residues
      keep_range_list (list): list of lists of the matching residues range,
        a list of ranges like [[1,3],[7-9]...]
    """
    id = ''.join(sorted([id1,id2]))
    if id in self.offset_dict.keys():
      offset = self.offset_dict[id]
    else:
      print >>self.log,"NO ID in matches_approx: ",id
      offset = 0

    start2_use=start2+offset  # equivalent residue number for start2
    keep_range_list = [] # a list of ranges like [[1,3],[7,9]...]

    if (not chain1) or (not chain2): return 0,keep_range_list
    end1 = start1 + (len(chain1)-1) * self.njump
    end2 = start2_use + (len(chain2)-1) * self.njump
    start = start1
    end = end1
    if start2_use > start1: start = start2_use
    if end2 < end1: end = end2
    # we start,end at residue number start,end, numbered according to chain 1

    start_njump = start//self.njump
    end_njump = end//self.njump
    if self.verbose:
       print >>self.log,"Finding matching segments in range from ",\
       start," to ",end," for chains ",chain1,chain2
    # iteratively select the longest segment that has ends that match with
    # at least 2 consecutive residues matching at each end...
    # and >= min_percent and >= min_length until there is none...
    # then sort in order of first residue number

    # if exact_match_only then skip if first residues or last do not match...

    matched_residues = []
    for i in xrange(start,end + self.njump, self.njump):
      res1 = chain1[(i-start1)//self.njump]
      res2 = chain2[(i-start2_use)//self.njump]
      if res1 and res2:
        if res1==res2:    #Matching
          matched_residues.append(True)
        else:
          matched_residues.append(False)
          if exact_match_only:
            return 0,keep_range_list
      else:
          matched_residues.append(False)

    segments = []
    still_looking=True
    possible=len(matched_residues)
    matching = 0
    while still_looking:
      best_segment,best_length = self.find_best_segment(
        matched_residues,segments,min_percent,min_length)
      if best_segment is None:
         still_looking=False
      else:
         segments.append(best_segment)
         matching += best_length
    segments.sort()
    if self.verbose:
       print >>self.log,"List of segments: ",segments

    # sort and translate segments into residue numbers...
    for [pos_start,pos_end] in segments:
      keep_range_list.append(
         [pos_start * self.njump + start,pos_end * self.njump + start])

    if self.verbose:
     print >>self.log,"Possible: ",possible," Matching: ",matching
     print >>self.log,"Residue range used: ",segments
    return matching,keep_range_list

  def have_enough(self,chain_ids,match_list):
    """ are all ids matched to a set of id #1? """
    if len(chain_ids)<2: return True
    for id in chain_ids[1:]:
      found=False
      for [matches,id1,id2,residue_range] in match_list:
        if id1==chain_ids[0] and id2==id: found=True
      if not found:
        if self.verbose: print >> self.log,"missing ",id,match_list
        return False
    if self.verbose: print >> self.log,"have ",chain_ids,match_list
    # all the id's are matched to id #1
    return True

  def find_best_segment(self,
                        matched_residues,
                        segments,
                        min_percent,
                        min_length):
    """ mask out all the residues that are used

    Returns:
      [best_start,best_end]: the start and end of a segment
      best_length (int): segment length
      """
    available_residues=matched_residues
    for [start,end] in segments:
      for i in xrange(start,end+1):
        available_residues[i]=False
    if self.verbose:
      print >>self.log,"Available_residues: ",available_residues

    # find longest segment in what remains, requiring matches on each
    # end and percent >= min_percent and length >= min_length
    # NOTE: min_length applies to absolute length; use min_length//njump
    best_start = 0
    best_end = 0
    best_length = 0
    for start in xrange(0,len(available_residues)):
      if not available_residues[start]:continue
      for end in xrange(start,len(available_residues)):
        if not available_residues[end]:continue
        n_match=available_residues[start:end+1].count(True)
        n_not_match=available_residues[start:end+1].count(False)
        n_total=n_match+n_not_match
        if n_total < min_length//self.njump: continue
        if not n_total:
          percent = 0.
        else:
          percent = 100. * float(n_match) / float(n_total)
        if percent < min_percent: continue
        if n_match>best_length:
           best_length = n_match
           best_start = start
           best_end = end
    if self.verbose:
      print >>self.log,"Best segment: ",best_start,best_end,best_length
    if best_length:
      return [best_start,best_end],best_length
    else:
      return None,0

def get_ncs_object_from_pdb(pdb_inp=None,hierarchy=None):
  """ Build and return ncs object from pdb hierarchy """
  if (not pdb_inp) and (not hierarchy):
    raise Sorry('No input is provided. Please provide pdb_inp or hierarchy.\n')
  elif pdb_inp and (not hierarchy):
    hierarchy=pdb_inp.hierarchy

  # build the ncs processing and output objects
  ncs_object = ncs(exclude_d=True,exclude_h=True)
  ncs_process = ncs_from_pdb()

  chains,chain_ids,starting_residue_numbers,offset_dict, total_residues \
    = ncs_process.get_chain_list(hierarchy=hierarchy,ignore_chains=[])

  # Get NCS groups using only sequence information and the suggested_ncs_groups
  sequence_groups,sequence_list_of_residue_range_list =\
    ncs_process.find_groups(
      hierarchy=hierarchy,
      chains=chains,
      chain_ids=chain_ids,
      starting_residue_numbers=starting_residue_numbers,
      min_length=1,
      min_percent=0.95,
      max_rmsd=2.0,
      max_rmsd_user=3.0,
      called_by_self=True,
      quick=False)

  # Get new groups with domains on all sequence groups
  invariant_groups,invariant_list_of_residue_range_list = \
    ncs_process.find_invariant_groups(
      hierarchy,sequence_groups,sequence_list_of_residue_range_list,
      chains,chain_ids,starting_residue_numbers)

  if invariant_groups:
    groups = invariant_groups
    list_of_residue_range_list = invariant_list_of_residue_range_list
    # build output object
    for group,residue_range_list in zip(groups,list_of_residue_range_list):
      rmsd_list,r_list,trans_list,center_list,residues_in_common_list = \
        ncs_process.get_rmsd(
          group,hierarchy,chains,chain_ids,
          starting_residue_numbers,residue_range_list)

      residue_range_list_with_offsets = \
        ncs_process.add_offsets(residue_range_list,group)

      chain_residue_id=[group,residue_range_list_with_offsets]
      ncs_object.import_ncs_group(ncs_rota_matr=r_list,
       rmsd_list=rmsd_list,
       residues_in_common_list=residues_in_common_list,
       center_orth=center_list,
       trans_orth=trans_list,
       chain_residue_id=chain_residue_id)
    return ncs_object

