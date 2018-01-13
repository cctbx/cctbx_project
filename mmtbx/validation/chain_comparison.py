from __future__ import division

# chain_comparison.py
# a tool to compare main-chain from two structures with or without crystal
# symmetry
#


import iotbx.phil
import sys,os
from libtbx.utils import Sorry,null_out
from scitbx.array_family import flex
from copy import deepcopy

master_phil = iotbx.phil.parse("""

  input_files {
    pdb_in = None
      .type = path
      .multiple = True
      .help = Input PDB file (enter once for target and once for query unless \
              query_dir is set)
      .short_caption = Input PDB file

    unique_query_only = False
      .type = bool
      .help = Use only unique chains in query. Normally use \
         unique_query_only=False and unique_part_of_target_only=True.
      .short_caption = Unique query only

    unique_target_pdb_in = None
      .type = path
      .help = Target model identifying which element is selected with \
           unique_query_only. NOTE: must be specified by keyword.
      .short_caption = Target model

    unique_part_of_target_only = True
      .type = bool
      .help = Use only unique chains in target (see also unique_query_only). \
      .short_caption = Unique target only

    test_unique_part_of_target_only = None
      .type = bool
      .help = Try both unique_part_of_target_only as True and False and \
             report result for whichever gives higher value of \
              fraction matching. Default is True if ncs_file is provided.
      .short_caption = Test unique target only

    allow_extensions = False
      .type = bool
      .help = If True, ignore parts of chains that do not overlap.  Normally \
              use False: identity is identity of overlapping part times \
              fraction of chain that overlaps.
      .short_caption = Allow extensions

    ncs_file = None
      .type = path
      .help = Select unique part of query. \
               If unique_query_only is False (typically) \
               apply NCS to it to generate full query.  Normally used with \
               test_unique_part_of_target_only=True.

    query_dir = None
      .type = path
      .help = directory containing query PDB files (any number)
      .short_caption = Query directory
  }
  output_files {

   match_pdb_file = None
     .type = path
     .help = Output file containing segments with specified match percentage
     .short_caption = Match PDB
  }
  crystal_info {
    chain_type = *PROTEIN RNA DNA
      .type = choice
      .short_caption = Chain type
      .help = Chain type.  All residues of other chain types ignored.

    use_crystal_symmetry = None
      .type = bool
      .short_caption = Use crystal symmetry in comparison
      .help = Default is True if space group is not P1.  \
              If set, use crystal symmetry to map atoms to closest positions
  }
  comparison {
    max_dist = 3.
      .type = float
      .short_caption = Maximum close distance
      .help= Maximum distance between atoms to be considered close

    distance_per_site = None
      .type = float
      .short_caption = Maximum distance spanned by a pair of residues
      .help =Maximum distance spanned by a pair of residues.  Set by \
            default as 3.8 A for protein and 8 A for RNA

    min_similarity = 0.9
      .type = float
      .short_caption = Minimum similarity in chains for uniqueness
      .help = When choosing unique chains, use min_similarity cutoff

    target_length_from_matching_chains = False
      .type = bool
      .short_caption = Use matching chains to get length
      .help= Use length of chains in target that are matched to \
                define full target length (as opposed to all unique \
                chains in target).

    minimum_percent_match_to_select = None
      .type = float
      .help = You can specify minimum_percent_match_to_select and \
              maximum_percent_match_to_select and match_pdb_file \
              in which case all segments in the query model that have \
              a percentage match (within max_dist of atom in target) \
              in this range will be written out to match_pdb_file.
    maximum_percent_match_to_select = None
      .type = float
      .help = You can specify minimum_percent_match_to_select and \
              maximum_percent_match_to_select and match_pdb_file \
              in which case all segments in the query model that have \
              a percentage match (within max_dist of atom in target) \
              in this range will be written out to match_pdb_file.
  }
  control {
      verbose = False
        .type = bool
        .help = Verbose output
        .short_caption = Verbose output

      quiet = False
        .type = bool
        .help = No printed output
        .short_caption = No printed output

  }
""", process_includes=True)
master_params = master_phil

class rmsd_values:
  def __init__(self):
    self.id_list=[]
    self.rmsd_list=[]
    self.n_list=[]
    self.match_percent_list=[]
    self.target_length_list=[]
    self.ratio_unique_to_total_target=None
    self.total_target=None
    self.total_query=None
    self.used_target=None
    self.used_query=None

  def add_match_percent(self,id=None,match_percent=None):
    ipoint=self.id_list.index(id)
    self.match_percent_list[ipoint]=match_percent

  def add_target_length(self,id=None,target_length=None):
    ipoint=self.id_list.index(id)
    self.target_length_list[ipoint]=target_length

  def add_rmsd(self,id=None,rmsd=None,n=None):
    self.id_list.append(id)
    self.rmsd_list.append(rmsd)
    self.n_list.append(n)
    self.match_percent_list.append(0)
    self.target_length_list.append(0)

  def get_match_percent(self,id=None):
    for local_id,local_match_percent in zip(
       self.id_list,self.match_percent_list):
      if id==local_id:
        return local_match_percent
    return 0

  def get_target_length(self,id=None):
    for local_id,local_target_length in zip(
       self.id_list,self.target_length_list):
      if id==local_id:
        return local_target_length
    return 0

  def get_close_to_target_percent(self,id=None):
    target_length=self.get_target_length(id=id)
    rmsd,n=self.get_values(id=id)
    if target_length is not None and n is not None:
      value=100.*n/max(1.,target_length)
      if self.ratio_unique_to_total_target:
         value=value/self.ratio_unique_to_total_target
      return value


    else:
      return 0.

  def get_values(self,id=None):
    for local_id,local_rmsd,local_n in zip(
       self.id_list,self.rmsd_list,self.n_list):
      if id==local_id:
        return local_rmsd,local_n
    return 0,0


def get_params(args,out=sys.stdout):
    command_line = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil=master_phil,
      pdb_file_def="input_files.pdb_in")

    params = command_line.work.extract()
    print >>out,"\nFind similarity between two main-chains"
    master_phil.format(python_object=params).show(out=out)
    return params

def best_match(sites1,sites2,crystal_symmetry=None,
     reject_if_too_far=None,distance_per_site=None):
  assert distance_per_site is not None
  # if reject_if_too_far and the centers of the two are further than can
  #  be reached by the remainders, skip

  unit_cell=crystal_symmetry.unit_cell()
  sps=crystal_symmetry.special_position_settings(min_distance_sym_equiv=0.5)
  from scitbx.array_family import flex

  # Match coordinates
  from cctbx import sgtbx

  # check central atoms if n>5 for each
  if sites1.size()>5 and sites2.size()>5:
    # what is distance?
    index1=sites1.size()//2
    index2=sites2.size()//2
    x1_ses=sps.sym_equiv_sites(site=sites1[index1])
    info=sgtbx.min_sym_equiv_distance_info(reference_sites=x1_ses,
           other=sites2[index2])
    dd=info.dist()

    # what is distance spannable by ends of each?
    max_dist=(index1+index2)*distance_per_site
    if dd > max_dist:
      info.i=index1
      info.j=index2
      return info  # hopeless

  best_info=None
  best_dist=None
  i=0
  for site in sites1:
    x1_ses=sps.sym_equiv_sites(site=site)
    j=0
    for site2 in sites2:
      info=sgtbx.min_sym_equiv_distance_info(reference_sites=x1_ses,
           other=site2)
      dd=info.dist()
      if best_dist is None or dd<best_dist:
         best_dist=dd
         best_info=info
         best_info.i=i  # just tack them on
         best_info.j=j
      j+=1
    i+=1
  return best_info

def select_atom_lines(hierarchy):
  lines=[]
  for line in hierarchy.as_pdb_string().splitlines():
    if line.startswith("ATOM "): # XXX NOTE USING PDB FORMAT HERE
      line=line.strip()
      lines.append(line)
  return lines

def get_best_match(xyz1,xyz2,crystal_symmetry=None,
    distance_per_site=None,used_j_list=None,removed_j=False):
  if crystal_symmetry:
    assert distance_per_site is not None
    info=best_match(
      xyz1,xyz2,
      crystal_symmetry=crystal_symmetry,
      distance_per_site=distance_per_site)
  else: # do it without symmetry
    (distance,i,j)=xyz1.min_distance_between_any_pair_with_id(xyz2)
    from libtbx import group_args
    info=group_args(i=i,j=j,distance=distance)

  if (not removed_j) and used_j_list and info.j in used_j_list:
    # move atom j away and try again
    if not distance_per_site:
      distance_per_site=4. # just need to move it away
    xyz2_new=xyz2.deep_copy()
    new_value=[]
    for x in xyz2_new[info.j]:
      new_value.append(x+distance_per_site*2.)
    xyz2_new[info.j]=tuple(new_value)

    return get_best_match(xyz1,xyz2_new,crystal_symmetry=crystal_symmetry,
         distance_per_site=distance_per_site,used_j_list=used_j_list+[info.j],
         removed_j=True)

  return info

def get_pdb_inp(text=None,file_name=None,source_info="string"):
  import iotbx.pdb
  if file_name:
    text=open(file_name).read()
    source_info="file %s" %(file_name)
  elif not text:
    text=""
  from cctbx.array_family import flex
  return iotbx.pdb.input(source_info=source_info,
       lines=flex.split_lines(text))

def get_chains_from_lines(lines):
  chains=[]
  for line in lines:
    h=get_pdb_inp(line).construct_hierarchy()
    id=h.models()[0].chains()[0].id
    if not id in chains:
      chains.append(id)
  return chains

def get_seq_from_lines(lines):
  seq=[]
  for line in lines:
    h=get_pdb_inp(line).construct_hierarchy()
    for res in h.models()[0].chains()[0].residues():
      seq.append(res.resname)
  return seq

def get_match_percent(seq1,seq2):
  assert len(seq1)==len(seq2)
  assert len(seq1)>0
  match_n=0
  for a,b in zip(seq1,seq2):
    if a.replace(" ","")==b.replace(" ",""): match_n+=1
  match_percent=100.*match_n/len(seq1)
  return match_n,match_percent

def apply_atom_selection(atom_selection,hierarchy=None):
  asc=hierarchy.atom_selection_cache()
  sel = asc.selection(string = atom_selection)
  return hierarchy.deep_copy().select(sel)  # deep copy is required
  #return hierarchy.select(sel)

def extract_representative_chains_from_hierarchy(ph,
    min_similarity=0.90,
    allow_extensions=False,
    allow_mismatch_in_number_of_copies=False,out=sys.stdout):

  unique_ph=extract_unique_part_of_hierarchy(ph,
    min_similarity=min_similarity,
    allow_mismatch_in_number_of_copies=allow_mismatch_in_number_of_copies,
    out=out)
  if ph==unique_ph:  # nothing is unique...
    print >>out,"No representative unique chains available"
    return None

  biggest_chain=None
  longest_chain=None
  for model in unique_ph.models()[:1]:
    for chain in model.chains():
      try:
        target_seq=chain.as_padded_sequence()  # has XXX for missing residues
        target_seq.replace("X","")
        chain_length=len(target_seq)
      except Exception, e:
        chain_length=0
      if chain_length and (longest_chain is None or chain_length>longest_chain):
        longest_chain=chain_length
        biggest_chain=chain
  if biggest_chain is None:
    print >>out,"Unable to extract unique part of hierarchy"
    return None

  biggest_chain_hierarchy=iotbx.pdb.input(
    source_info="Model",lines=flex.split_lines("")).construct_hierarchy()
  mm=iotbx.pdb.hierarchy.model()
  biggest_chain_hierarchy.append_model(mm)
  mm.append_chain(biggest_chain.detached_copy())

  # ready with biggest chain

  copies_of_biggest_chain_ph=extract_copies_identical_to_target_from_hierarchy(
     ph,target_ph=biggest_chain_hierarchy,
     allow_extensions=False,out=sys.stdout)
  return copies_of_biggest_chain_ph

def extract_copies_identical_to_target_from_hierarchy(ph,
     allow_extensions=False,
     min_similarity=None,target_ph=None,out=sys.stdout):
  new_hierarchy=iotbx.pdb.input(
    source_info="Model",lines=flex.split_lines("")).construct_hierarchy()
  mm=iotbx.pdb.hierarchy.model()
  new_hierarchy.append_model(mm)

  assert target_ph is not None
  target_seq=None
  for model in target_ph.models()[:1]:
    for chain in model.chains():
      try:
        target_seq=chain.as_padded_sequence()  # has XXX for missing residues
        target_seq.replace("X","")
        break
      except Exception, e:
        pass
  if not target_seq:
     raise Sorry("No sequence found in target sequence for "+
       "extract_copies_identical_to_target_from_hierarchy")

  matching_chain_list=[]
  for model in ph.models()[:1]:
    for chain in model.chains():
      try:
        seq=chain.as_padded_sequence()  # has XXX for missing residues
        seq=seq.replace("X","")
      except Exception, e:
        seq=""
      similar_seq=seq_it_is_similar_to(
         seq=seq,unique_sequences=[target_seq],
         allow_extensions=allow_extensions,
         min_similarity=min_similarity)  # check for similar...
      if similar_seq:
        matching_chain_list.append(chain)

  total_chains=0
  for chain in matching_chain_list:
    mm.append_chain(chain.detached_copy())
    total_chains+=1
  print >>out,"Total chains extracted: %s" %(total_chains)
  return new_hierarchy

def seq_it_is_similar_to(seq=None,unique_sequences=None,min_similarity=1.0,
   allow_extensions=False):
  from phenix.loop_lib.sequence_similarity import sequence_similarity
  for s in unique_sequences:
    sim=sequence_similarity().run(seq,s,use_fasta=True,verbose=False)
    if not allow_extensions and len(seq)!=len(s):
      fract_same=min(len(seq),len(s))/max(len(seq),len(s))
      sim=sim*fract_same
    if sim >= min_similarity:
      return s # return the one it is similar to
  return None

def extract_unique_part_of_sequences(sequence_list=None,
    allow_mismatch_in_number_of_copies=True,
    allow_extensions=False,
    min_similarity=1.0,out=sys.stdout):

  unique_sequences=[]
  best_chain_copies_dict={}

  unique_sequence_dict={}  # unique_sequence_dict[seq] is unique sequence sim
                           # to seq

  for seq in sequence_list:
      if not seq: continue
      similar_seq=seq_it_is_similar_to(
         seq=seq,unique_sequences=unique_sequences,
         allow_extensions=allow_extensions,
         min_similarity=min_similarity)  # check for similar...
      if similar_seq:
        unique_sequence_dict[seq]=similar_seq
        seq=similar_seq
      else:
        unique_sequences.append(seq)
        unique_sequence_dict[seq]=seq
      if not seq in best_chain_copies_dict:
        best_chain_copies_dict[seq]=0
      best_chain_copies_dict[seq]+=1

  copies_list=[]
  for seq in unique_sequences:
    copies_found=best_chain_copies_dict[seq]
    if not copies_found in copies_list: copies_list.append(copies_found)
  copies_list.sort()
  print >>out,"Numbers of copies of sequences: %s" %(str(copies_list))
  copies_base=copies_list[0]
  all_ok=True
  copies_in_unique={}
  if len(copies_list)==1:
    print >>out,"Number of copies of all sequences is: %s" %(copies_base)
    for seq in unique_sequences:
      copies_in_unique[seq]=1 # unique set has 1 of this one sequence
    return copies_in_unique,copies_base,unique_sequence_dict
  else:
    for cf in copies_list[1:]:
      if cf//copies_base != cf/copies_base:  # not integral
        all_ok=False
        break
    if all_ok:
      print >>out,"Copies are all multiples of %s...taking " %(
          copies_base)
      for seq in unique_sequences:
        copies_found=best_chain_copies_dict[seq]
        copies_in_unique[seq]=copies_found//copies_base
      return copies_in_unique,copies_base,unique_sequence_dict
    else:
      print >>out,"Copies are not all multiples of %s...taking all" %(
          copies_base)
      for seq in unique_sequences:
        copies_in_unique[seq]=best_chain_copies_dict[seq]
      return copies_in_unique,copies_base,unique_sequence_dict

def get_matching_ids(unique_seq=None,sequences=None,
      chains=None,unique_sequence_dict=None):

    matching_ids=[]
    matching_chains=[]
    for s1,c1 in zip(sequences,chains):
      if unique_sequence_dict[s1]==unique_seq:
          matching_ids.append(c1.id)
          matching_chains.append(c1)
    return matching_ids,matching_chains

def get_sorted_matching_chains(
   chains=None,
   target_centroid_list=None):
  sort_list=[]
  for chain in chains:
    if target_centroid_list:
        xx=flex.vec3_double()
        xx.append(chain.atoms().extract_xyz().mean())
        dist=xx.min_distance_between_any_pair(target_centroid_list)
    else:
        dist=0.
    sort_list.append([dist,chain])
  sort_list.sort()
  sorted_chains=[]
  sorted_distances=[]
  for dist,chain in sort_list:
    sorted_chains.append(chain)
    sorted_distances.append(dist)
  return sorted_chains,sorted_distances

def extract_unique_part_of_hierarchy(ph,target_ph=None,
    allow_mismatch_in_number_of_copies=True,
    allow_extensions=False,
    min_similarity=1.0,out=sys.stdout):

  # Container for unique chains:

  new_hierarchy=iotbx.pdb.input(
    source_info="Model",lines=flex.split_lines("")).construct_hierarchy()
  mm=iotbx.pdb.hierarchy.model()
  new_hierarchy.append_model(mm)

  # Target location:

  if target_ph:
    target_centroid_list=flex.vec3_double()
    for model in target_ph.models()[:1]:
      target_centroid_list.append(model.atoms().extract_xyz().mean())
  else:
    target_centroid_list=None

  # Get unique set of sequences
  # Also save all the chains associated with each one

  sequences=[]
  chains=[]
  for model in ph.models()[:1]:
    for chain in model.chains():
      try:
        seq=chain.as_padded_sequence()  # has XXX for missing residues
        seq=seq.replace("X","")
        sequences.append(seq)
        chains.append(chain)
      except Exception, e:
        pass
  copies_in_unique,base_copies,unique_sequence_dict=\
        extract_unique_part_of_sequences(
    sequence_list=sequences,
    allow_mismatch_in_number_of_copies=allow_mismatch_in_number_of_copies,
    allow_extensions=allow_extensions,
    min_similarity=min_similarity,out=out)

  sequences_matching_unique_dict={}
  for seq in sequences:
    unique_seq=unique_sequence_dict[seq]
    if not unique_seq in sequences_matching_unique_dict.keys():
      sequences_matching_unique_dict[unique_seq]=[]
    sequences_matching_unique_dict[unique_seq].append(seq)

  # Now we are going to return the unique set...if choice of which copies,
  #  take those closest to the target (in that order)

  unique_chains=[]

  print >>out,"Unique set of sequences. Copies of the unique set: %s" %(
      base_copies)
  print >>out,"Copies in unique set ID  sequence"
  for unique_seq in copies_in_unique.keys():
    matching_ids,matching_chains=get_matching_ids(
      unique_seq=unique_seq,sequences=sequences,chains=chains,
      unique_sequence_dict=unique_sequence_dict)
    print >>out," %s    %s    %s " %(
      copies_in_unique[unique_seq]," ".join(matching_ids),unique_seq)
    sorted_matching_chains,sorted_matching_distances=get_sorted_matching_chains(
      chains=matching_chains,
      target_centroid_list=target_centroid_list)
    for chain,dist in zip(
        sorted_matching_chains[:copies_in_unique[unique_seq]],
        sorted_matching_distances[:copies_in_unique[unique_seq]]):
      mm.append_chain(chain.detached_copy())
      print >>out, "Adding chain %s: %s (%s): %7.2f" %(
         chain.id,unique_seq,str(chain.atoms().extract_xyz()[0]),
         dist)
  return new_hierarchy

def run_test_unique_part_of_target_only(params=None,
       out=sys.stdout,
       ncs_obj=None,
       target_hierarchy=None,
       chain_hierarchy=None,
       target_file=None, # model
       chain_file=None, # query
       crystal_symmetry=None,
       max_dist=None,
       quiet=None,
       verbose=None,
       use_crystal_symmetry=None,
       chain_type=None,
       target_length_from_matching_chains=None,
       distance_per_site=None,
       min_similarity=None):
  if params.control.verbose:
    local_out=out
  else:
    local_out=null_out()
  rv_list=[]
  file_list=[]
  best_rv=None
  best_t=None
  best_percent_close=None
  best_rv_reverse=None
  for t in [True,False]:
    rv_reverse=None
    for reverse in [True,False]:
      local_params=deepcopy(params)
      local_params.input_files.test_unique_part_of_target_only=False
      local_params.input_files.unique_part_of_target_only=t
      if reverse:
       local_chain_hierarchy=target_hierarchy
       local_target_hierarchy=chain_hierarchy
      else:
       local_chain_hierarchy=chain_hierarchy
       local_target_hierarchy=target_hierarchy

      rv=run(params=local_params,out=local_out,
          ncs_obj=ncs_obj,
          target_hierarchy=local_target_hierarchy,
          chain_hierarchy=local_chain_hierarchy,
          crystal_symmetry=crystal_symmetry,
          max_dist=max_dist,
          quiet=quiet,
          verbose=verbose,
          use_crystal_symmetry=use_crystal_symmetry,
          chain_type=chain_type,
          target_length_from_matching_chains=target_length_from_matching_chains,
          distance_per_site=distance_per_site,
          min_similarity=min_similarity,
          )
      if reverse: # save number far from target
        rv_reverse=rv
      else:
        percent_close=rv.get_close_to_target_percent('close')
        print >>out,"Percent close with unique_part_of_target_only=%s: %7.1f" %(
          t,percent_close)
        if best_percent_close is None or percent_close>best_percent_close:
          best_percent_close=percent_close
          best_rv=rv
          best_t=t
  print >>out,"\nOriginal residues in target: %s  In query: %s" %(
    best_rv.total_target,best_rv.total_chain)
  print >>out,"Used residues in target:  %s  In query: %s" %(
    best_rv.used_target,best_rv.used_chain)
  rv_list=[best_rv]
  if best_t:
    file_list=['Unique_target']
  else:
    file_list=['Entire_target']
  write_summary(params=params,file_list=file_list,rv_list=rv_list, out=out)

def run_all(params=None,
       out=sys.stdout,
       ncs_obj=None,
       target_hierarchy=None,
       chain_hierarchy=None,
       target_file=None, # model
       chain_file=None, # query
       crystal_symmetry=None,
       max_dist=None,
       quiet=None,
       verbose=None,
       use_crystal_symmetry=None,
       chain_type=None,
       target_length_from_matching_chains=None,
       distance_per_site=None,
       min_similarity=None):

  if params.control.verbose:
    local_out=out
  else:
    local_out=null_out()
  rv_list=[]
  file_list=[]
  for query_file in os.listdir(params.input_files.query_dir):
    file_name=os.path.join(params.input_files.query_dir,query_file)
    if not os.path.isfile(file_name): continue
    local_params=deepcopy(params)
    local_params.input_files.query_dir=None
    local_params.input_files.pdb_in.append(file_name)
    try:
      rv=run(params=local_params,out=local_out,
        ncs_obj=ncs_obj,
        target_hierarchy=target_hierarchy,
        chain_hierarchy=chain_hierarchy,
        crystal_symmetry=crystal_symmetry,
        max_dist=max_dist,
        quiet=quiet,
        verbose=verbose,
        use_crystal_symmetry=use_crystal_symmetry,
        chain_type=chain_type,
        target_length_from_matching_chains=target_length_from_matching_chains,
        distance_per_site=distance_per_site,
        min_similarity=min_similarity,
      )

    except Exception,e:
      if str(e).find("CifParserError"):
        print >>out,"NOTE: skipping %s as it is not a valid model file" %(
           file_name)
        continue # it was not a valid PDB file...just skip it
      else:
        raise Sorry(str(e))
    rv_list.append(rv)
    file_list.append(file_name)

  write_summary(params=params,file_list=file_list,rv_list=rv_list, out=out)

def write_summary(params=None,file_list=None,rv_list=None,
    max_dist=None,write_header=True,out=sys.stdout):

  if params and max_dist is None:
     max_dist=params.comparison.max_dist
  if max_dist is None: max_dist=0.

  if write_header:
    print >>out,"\nCLOSE is within %4.1f A. FAR is greater than this." %(
      max_dist)
    print >>out,"\nCA SCORE is fraction in close CA / rmsd of these CA."
    print >>out,"\nSEQ SCORE is fraction (close and matching target sequence).\n"

    print >>out,"\n"
    print >>out,"               ----ALL RESIDUES---  CLOSE RESIDUES ONLY    %"
    print >>out,\
              "     MODEL     --CLOSE-    --FAR-- FORWARD REVERSE MIXED"+\
              " FOUND  CA                  SEQ"
    print >>out,"               RMSD   N      N       N       N      N  "+\
              "        SCORE  SEQ MATCH(%)  SCORE"+"\n"

  results_dict={}
  score_list=[]
  for rv,full_f in zip(rv_list,file_list):
    results_dict[full_f]=rv
    (rmsd,n)=rv.get_values('close')
    target_length=rv.get_target_length('close')
    score=n/(max(1,target_length)*max(0.1,rmsd))
    if rv.ratio_unique_to_total_target:
         score=score/rv.ratio_unique_to_total_target
    score_list.append([score,full_f])
  score_list.sort()
  score_list.reverse()
  for score,full_f in score_list:
    rv=results_dict[full_f]
    percent_close=rv.get_close_to_target_percent('close')

    seq_score=rv.get_match_percent('close')*percent_close/10000
    file_name=os.path.split(full_f)[-1]
    close_rmsd,close_n=rv.get_values('close')
    if not close_rmsd: close_rmsd=0
    far_away_rmsd,far_away_n=rv.get_values('far_away')
    forward_rmsd,forward_n=rv.get_values('forward')
    reverse_rmsd,reverse_n=rv.get_values('reverse')
    unaligned_rmsd,unaligned_n=rv.get_values('unaligned')
    match_percent=rv.get_match_percent('close')
    print >>out,"%14s %4.2f %4d   %4d   %4d    %4d    %4d  %5.1f %6.2f   %5.1f      %6.2f" %(file_name,close_rmsd,close_n,far_away_n,forward_n,
         reverse_n,unaligned_n,percent_close,score,match_percent,seq_score)

def get_target_length(target_chain_ids=None,hierarchy=None,
     target_length_from_matching_chains=None):
  total_length=0  # just counts residues
  for model in hierarchy.models()[:1]:
    for chain in model.chains():
      if (not target_length_from_matching_chains) or chain.id in target_chain_ids:
        for conformer in chain.conformers()[:1]:
          total_length+=len(conformer.residues())
  return total_length

def select_segments_that_match(params=None,
   chain_hierarchy=None,
   target_hierarchy=None,
   out=sys.stdout,
   ncs_obj=None,
   target_file=None, # model
   chain_file=None, # query
   crystal_symmetry=None,
   max_dist=None,
   quiet=None,
   verbose=None,
   use_crystal_symmetry=None,
   chain_type=None,
   target_length_from_matching_chains=None,
   distance_per_site=None,
   min_similarity=None):


  # Identify all the segments in chain_hierarchy that match target_hierarchy
  #  and write them out
  from mmtbx.secondary_structure.find_ss_from_ca import split_model,model_info,\
    merge_hierarchies_from_models
  chain_model=model_info(hierarchy=chain_hierarchy)
  if params.crystal_info.chain_type=="PROTEIN":
    distance_cutoff=5.
  else:
    distance_cutoff=15.
  chain_models=split_model(model=chain_model,distance_cutoff=distance_cutoff)
  print >>out,"Analyzing %s segments and identifying " %(len(chain_models)) +\
      " those with "+\
     "chain_type=%s and match percentage between %.1f %% and %.1f %% " %(
    params.crystal_info.chain_type,
    params.comparison.minimum_percent_match_to_select,
    params.comparison.maximum_percent_match_to_select)
  local_params=deepcopy(params)
  local_params.output_files.match_pdb_file=None # required
  models_to_keep=[]
  write_header=True
  for cm in chain_models:  # one segment
    rv_list=[]
    file_list=[]
    rv=run(
      params=local_params,
      ncs_obj=ncs_obj,
      target_hierarchy=target_hierarchy,
      quiet=True,
      chain_hierarchy=cm.hierarchy,out=null_out(),
        crystal_symmetry=crystal_symmetry,
        max_dist=max_dist,
        verbose=verbose,
        use_crystal_symmetry=use_crystal_symmetry,
        chain_type=chain_type,
        target_length_from_matching_chains=target_length_from_matching_chains,
        distance_per_site=distance_per_site,
        min_similarity=min_similarity,
      )

    rv_list.append(rv)
    file_list.append(params.crystal_info.chain_type)
    close_rmsd,close_n=rv.get_values('close')
    far_away_rmsd,far_away_n=rv.get_values('far_away')
    if close_n+far_away_n<1: continue # wrong chain type or other failure

    percent_matched=100.*close_n/max(1,close_n+far_away_n)
    if percent_matched < params.comparison.minimum_percent_match_to_select:
      continue
    if percent_matched > params.comparison.maximum_percent_match_to_select:
      continue

    write_summary(params=params,file_list=file_list,rv_list=rv_list,
      write_header=write_header,out=out)
    write_header=False
    models_to_keep.append(cm)

  new_model=merge_hierarchies_from_models(models=models_to_keep,resid_offset=5)
  ff=open(params.output_files.match_pdb_file,'w')
  print >>ff,new_model.hierarchy.as_pdb_string()
  ff.close()
  print >>out,"Wrote %s %s chains with %s residues to %s" %(
    len(models_to_keep),params.crystal_info.chain_type,
    new_model.hierarchy.overall_counts().n_residues,
    params.output_files.match_pdb_file)
  return new_model

def get_ncs_obj(file_name,out=sys.stdout):
  from mmtbx.ncs.ncs import ncs
  ncs_object=ncs()
  ncs_object.read_ncs(file_name=file_name,log=out)
  return ncs_object

def apply_ncs_to_hierarchy(ncs_obj=None,
        hierarchy=None,out=sys.stdout):
  if not ncs_obj or ncs_obj.max_operators()<2:
    return hierarchy
  try:
    from phenix.command_line.apply_ncs import apply_ncs as apply_ncs_to_atoms
  except exception, e:
    print "Need phenix for applying NCS"
    return hierarchy

  print >>out, "Applying NCS now..."
  from phenix.autosol.get_pdb_inp import get_pdb_hierarchy
  identity_copy=ncs_obj.identity_op_id_in_first_group()+1

  args=['pdb_out=None','match_copy=%s' %(identity_copy),
       'params_out=None' ]
  args.append("use_space_group_symmetry=False")
  an=apply_ncs_to_atoms(
      args,hierarchy=hierarchy,
      ncs_object=ncs_obj,
      out=out)
  new_hierarchy=get_pdb_hierarchy(text=an.output_text)
  return new_hierarchy

def run(args=None,
   ncs_obj=None,
   target_unique_hierarchy=None,
   target_hierarchy=None,
   chain_hierarchy=None,
   target_file=None, # model
   chain_file=None, # query
   crystal_symmetry=None,
   max_dist=None,
   quiet=None,
   verbose=None,
   use_crystal_symmetry=None,
   chain_type=None,
   params=None,
   target_length_from_matching_chains=None,
   distance_per_site=None,
   min_similarity=None,
   out=sys.stdout):
  if not args: args=[]
  if not params:
    params=get_params(args,out=out)
  if params.input_files.pdb_in:
    print >>out,"Using %s as target" %(params.input_files.pdb_in[0])
  elif chain_file or chain_hierarchy:
    pass # it is fine
  else:
    raise Sorry("Need target model (pdb_in)")
  if params.input_files.unique_query_only and \
     params.input_files.unique_part_of_target_only:
    print >>out,"Warning: You have specified unique_query_only and" +\
       " unique_part_of_target_only. \nThis is not normally appropriate "
  if params.input_files.unique_target_pdb_in and \
         params.input_files.unique_query_only:
    print >>out,"Using %s as target for unique chains" %(
       params.input_files.unique_target_pdb_in)
  if params.input_files.query_dir and \
      os.path.isdir(params.input_files.query_dir) and \
      not params.output_files.match_pdb_file:
    print >>out,"\nUsing all files in %s as queries\n" %(
       params.input_files.query_dir)
    return run_all(params=params,out=out)


  if not ncs_obj and params.input_files.ncs_file:
    ncs_obj=get_ncs_obj(params.input_files.ncs_file,out=out)
    print >>out,"NCS with %s operators read from %s" %(ncs_obj.max_operators(),
       params.input_files.ncs_file)
    if ncs_obj.max_operators()<2:
      print >>out,"Skipping NCS (no operators)"
      ncs_obj=None

  if verbose is None:
    verbose=params.control.verbose
  if quiet is None:
    quiet=params.control.quiet
  if chain_type is None:
    chain_type=params.crystal_info.chain_type
  if use_crystal_symmetry is None:
    use_crystal_symmetry=params.crystal_info.use_crystal_symmetry
  params.crystal_info.use_crystal_symmetry=use_crystal_symmetry
  if max_dist is None:
    max_dist=params.comparison.max_dist
  if distance_per_site is None:
    distance_per_site=params.comparison.distance_per_site
  if min_similarity is None:
    min_similarity=params.comparison.min_similarity
  if target_length_from_matching_chains is None:
    target_length_from_matching_chains=\
       params.comparison.target_length_from_matching_chains

  if verbose:
    local_out=out
  else:
    local_out=null_out()

  if not target_file and len(params.input_files.pdb_in)>0:
     target_file=params.input_files.pdb_in[0]  # model
  if not chain_file and len(params.input_files.pdb_in)>1:
     chain_file=params.input_files.pdb_in[1] # query

  # get the hierarchies
  if not chain_hierarchy or not target_hierarchy:
    assert chain_file and target_file
    pdb_inp=get_pdb_inp(file_name=chain_file)
    if params.input_files.unique_target_pdb_in:
      target_unique_hierarchy=get_pdb_inp(
        file_name=params.input_files.unique_target_pdb_in).construct_hierarchy()
    if not crystal_symmetry:
      crystal_symmetry=pdb_inp.crystal_symmetry_from_cryst1()
    chain_hierarchy=pdb_inp.construct_hierarchy()
    target_pdb_inp=get_pdb_inp(file_name=target_file)
    if not crystal_symmetry or not crystal_symmetry.unit_cell():
      crystal_symmetry=target_pdb_inp.crystal_symmetry_from_cryst1()
    target_hierarchy=target_pdb_inp.construct_hierarchy()
    # remove hetero atoms as they are not relevant
    chain_hierarchy=apply_atom_selection('not hetero',chain_hierarchy)
    target_hierarchy=apply_atom_selection('not hetero',target_hierarchy)

  total_target=target_hierarchy.overall_counts().n_residues
  total_chain=chain_hierarchy.overall_counts().n_residues

  if params.input_files.test_unique_part_of_target_only or \
      (params.input_files.test_unique_part_of_target_only is None and
        params.input_files.ncs_file):
    print>>out,"\nTesting unique_part_of_target_only as True and False and "
    print >>out,"reporting results for whichever gives higher fraction matched."
    return run_test_unique_part_of_target_only(params=params,out=out,
          ncs_obj=ncs_obj,
          target_hierarchy=target_hierarchy,
          chain_hierarchy=chain_hierarchy,
          crystal_symmetry=crystal_symmetry,
          max_dist=max_dist,
          quiet=quiet,
          verbose=verbose,
          use_crystal_symmetry=use_crystal_symmetry,
          chain_type=chain_type,
          target_length_from_matching_chains=target_length_from_matching_chains,
          distance_per_site=distance_per_site,
          min_similarity=min_similarity)


  # Take unique part of query if requested
  if target_unique_hierarchy:
    target_ph=target_unique_hierarchy
  else:
    target_ph=chain_hierarchy

  if params.input_files.unique_query_only or ncs_obj:
    print >>out,"\nExtracting unique part of query\n"
    chain_hierarchy=extract_unique_part_of_hierarchy(
      chain_hierarchy,target_ph=target_ph,
      allow_extensions=params.input_files.allow_extensions,
      min_similarity=min_similarity,out=local_out)
    print >>out,"Residues in unique part of query hierarchy: %s" %(
     chain_hierarchy.overall_counts().n_residues)
    if ncs_obj and not params.input_files.unique_query_only:
      # apply NCS to unique part of query
      print >>out,"Applying NCS to unique part of query"
      chain_hierarchy=apply_ncs_to_hierarchy(ncs_obj=ncs_obj,
        hierarchy=chain_hierarchy,out=out)
      print >>out,"Residues in full query hierarchy: %s" %(
        chain_hierarchy.overall_counts().n_residues)

  if params.input_files.unique_part_of_target_only:
    print >>out,"\nUsing only unique part of target \n"
    print >>out,"Residues in input target hierarchy: %s" %(
     target_hierarchy.overall_counts().n_residues)
    target_hierarchy=extract_unique_part_of_hierarchy(
      target_hierarchy,target_ph=target_ph,
      allow_extensions=params.input_files.allow_extensions,
      min_similarity=min_similarity,out=local_out)
    print >>out,"Residues in unique part of target hierarchy: %s" %(
     target_hierarchy.overall_counts().n_residues)


  if params.output_files.match_pdb_file and \
    params.comparison.minimum_percent_match_to_select is not None and \
    params.comparison.maximum_percent_match_to_select is not None:
      return select_segments_that_match(params=params,
       chain_hierarchy=chain_hierarchy,
       target_hierarchy=target_hierarchy,out=out)

  if params.input_files.unique_part_of_target_only:
    unique_part_of_target_hierarchy=extract_unique_part_of_hierarchy(
        target_hierarchy,
        min_similarity=min_similarity,
        allow_extensions=params.input_files.allow_extensions,
        target_ph=target_unique_hierarchy,out=local_out)
    ratio_unique_to_total_target=\
       unique_part_of_target_hierarchy.overall_counts().n_residues/  \
       max(1,target_hierarchy.overall_counts().n_residues)
    print >>out,"Counting unique residues in target when calculating "+\
      "\npercentage built (fraction=%.2f)" %(ratio_unique_to_total_target)
  else:
    ratio_unique_to_total_target=1.
    print >>out,"Counting all residues in target when calculating "+\
      "percentage built"

  used_target=target_hierarchy.overall_counts().n_residues
  used_chain=chain_hierarchy.overall_counts().n_residues

  if params.crystal_info.use_crystal_symmetry is None: # set default
    if crystal_symmetry and crystal_symmetry.space_group() and \
       (not crystal_symmetry.space_group().type().number() in [0,1]):
      params.crystal_info.use_crystal_symmetry=True
    else:
      params.crystal_info.use_crystal_symmetry=False
      crystal_symmetry=None
  elif params.crystal_info.use_crystal_symmetry==False:
      crystal_symmetry=None

  if not crystal_symmetry or not crystal_symmetry.unit_cell():
    crystal_symmetry=get_pdb_inp(
        text="CRYST1 1000.000 1000.000 1000.000  90.00  90.00  90.00 P 1"
        ).crystal_symmetry_from_cryst1()
    print >>out,"\nCrystal symmetry will not be used in comparison.\n"
    if use_crystal_symmetry:
        raise Sorry("Please set use_crystal_symmetry"+
           "=False (no crystal symmetry supplied)")
  else:
    print >>out,"\nCrystal symmetry will be used in comparison.\n"
    print>>out, "Space group: %s" %(crystal_symmetry.space_group().info()), \
         "Unit cell: %7.2f %7.2f %7.2f  %7.2f %7.2f %7.2f \n" %(
        crystal_symmetry.unit_cell().parameters())
    use_crystal_symmetry=True
  if not quiet:
    print >>out,"Looking for chain similarity for "+\
      "%s (%d residues) in the model %s (%d residues)" %(
     chain_file,chain_hierarchy.overall_counts().n_residues,
     target_file,target_hierarchy.overall_counts().n_residues)
    if verbose:
      print >>out,"Chain type is: %s" %(chain_type)
  if crystal_symmetry is None or crystal_symmetry.unit_cell() is None:
    raise Sorry("Need crystal symmetry in at least one input file")
  # get the CA residues
  if chain_type in ["RNA","DNA"]:
    atom_selection="name P"
    if not distance_per_site:
      distance_per_site=8.
  else:
    atom_selection="name ca and (not element Ca)"
    if not distance_per_site:
      distance_per_site=3.8
  chain_ca=apply_atom_selection(atom_selection,chain_hierarchy)
  chain_ca_lines=select_atom_lines(chain_ca)
  target_ca=apply_atom_selection(atom_selection,target_hierarchy)
  target_xyz_lines=select_atom_lines(target_ca)
  chain_xyz_cart=chain_ca.atoms().extract_xyz()
  target_xyz_cart=target_ca.atoms().extract_xyz()

  if target_xyz_cart.size()<1:
    print "No suitable atoms in target"
    return rmsd_values()
  if chain_xyz_cart.size()<1:
    print "No suitable atoms in query"
    return rmsd_values()

  # for each xyz in chain, figure out closest atom in target and dist
  best_i=None
  best_i_dd=None
  best_pair=None
  pair_list=[]
  from scitbx.array_family import flex
  chain_xyz_fract=crystal_symmetry.unit_cell().fractionalize(chain_xyz_cart)
  target_xyz_fract=crystal_symmetry.unit_cell().fractionalize(target_xyz_cart)
  far_away_match_list=[]
  far_away_match_rmsd_list=flex.double()
  if use_crystal_symmetry:
    working_crystal_symmetry=crystal_symmetry
  else:
    working_crystal_symmetry=None
  used_j_list=[]
  for i in xrange(chain_xyz_fract.size()):
    best_j=None
    best_dd=None
    distance=None
    if working_crystal_symmetry:
      info=get_best_match(
        flex.vec3_double([chain_xyz_fract[i]]),target_xyz_fract,
        crystal_symmetry=working_crystal_symmetry,
        distance_per_site=distance_per_site,used_j_list=used_j_list)
      if info:
        distance=info.dist()
    else:
      info=get_best_match(
        flex.vec3_double([chain_xyz_cart[i]]),target_xyz_cart,
        used_j_list=used_j_list)
      distance=info.distance
    if info and (best_dd is None or distance<best_dd):
        best_dd=distance
        best_j=info.j
    if best_dd > max_dist:
      far_away_match_list.append(i)
      far_away_match_rmsd_list.append(best_dd**2)
      if (not quiet) and verbose:
        print >>out,"%s" %(chain_ca_lines[i])
      continue
    if best_i is None or best_dd<best_i_dd:
      best_i=i
      best_i_dd=best_dd
      best_pair=[i,best_j]
      used_j_list.append(best_j)
    pair_list.append([i,best_j,best_dd])
  n_forward=0
  n_reverse=0
  forward_match_list=[]
  reverse_match_list=[]
  unaligned_match_list=[]
  close_match_list=[]
  forward_match_rmsd_list=flex.double()
  reverse_match_rmsd_list=flex.double()
  unaligned_match_rmsd_list=flex.double()
  close_match_rmsd_list=flex.double()
  last_i=None
  last_j=None
  for [i,j,dd],[next_i,next_j,next_dd] in zip(
      pair_list,pair_list[1:]+[[None,None,None]]):
    if i is None or j is None: continue
    found=False
    if last_i is None: # first time
      if next_i==i+1: # starting a segment
        if next_j==j+1:
          n_forward+=1
          forward_match_list.append([i,j])
          close_match_list.append([i,j])
          forward_match_rmsd_list.append(dd**2)
          close_match_rmsd_list.append(dd**2)
          found=True
        elif next_j==j-1:
          n_reverse+=1
          reverse_match_list.append([i,j])
          close_match_list.append([i,j])
          reverse_match_rmsd_list.append(dd**2)
          close_match_rmsd_list.append(dd**2)
          found=True
    else: # not the first time
      if i==last_i+1: # continuing a segment
        if j==last_j+1:
          n_forward+=1
          forward_match_list.append([i,j])
          close_match_list.append([i,j])
          forward_match_rmsd_list.append(dd**2)
          close_match_rmsd_list.append(dd**2)
          found=True
        elif j==last_j-1:
          n_reverse+=1
          reverse_match_list.append([i,j])
          close_match_list.append([i,j])
          reverse_match_rmsd_list.append(dd**2)
          close_match_rmsd_list.append(dd**2)
          found=True
    if not found:
      last_i=None
      last_j=None
      unaligned_match_list.append([i,j])
      close_match_list.append([i,j])
      unaligned_match_rmsd_list.append(dd**2)
      close_match_rmsd_list.append(dd**2)
    else:
      last_i=i
      last_j=j

  if n_forward==n_reverse==0:
    direction='none'
  elif n_forward>= n_reverse:
    direction='forward'
  else:
    direction='reverse'
  if (not quiet) and verbose:
    print >>out,"%s %d  %d  N: %d" %(
     direction,n_forward,n_reverse,chain_xyz_fract.size())

  rv=rmsd_values()

  id='forward'
  if forward_match_rmsd_list.size():
    rmsd=forward_match_rmsd_list.min_max_mean().mean**0.5
  else:
    rmsd=None
  n=forward_match_rmsd_list.size()
  rv.add_rmsd(id=id,rmsd=rmsd,n=n)

  id='reverse'
  if reverse_match_rmsd_list.size():
    rmsd=reverse_match_rmsd_list.min_max_mean().mean**0.5
  else:
    rmsd=None
  n=reverse_match_rmsd_list.size()
  rv.add_rmsd(id=id,rmsd=rmsd,n=n)

  id='unaligned'
  if unaligned_match_rmsd_list.size():
    rmsd=unaligned_match_rmsd_list.min_max_mean().mean**0.5
  else:
    rmsd=None
  n=unaligned_match_rmsd_list.size()
  rv.add_rmsd(id=id,rmsd=rmsd,n=n)
  id='close'
  if close_match_rmsd_list.size():
      rmsd=close_match_rmsd_list.min_max_mean().mean**0.5
  else:
      rmsd=None
  n=close_match_rmsd_list.size()
  rv.add_rmsd(id=id,rmsd=rmsd,n=n)

  id='far_away'
  if far_away_match_rmsd_list.size():
      rmsd=far_away_match_rmsd_list.min_max_mean().mean**0.5
  else:
      rmsd=None
  n=far_away_match_rmsd_list.size()
  rv.add_rmsd(id=id,rmsd=rmsd,n=n)

  if not quiet:
    if verbose:
      print >>out,"Total CA: %d  Too far to match: %d " %(
        chain_xyz_fract.size(),len(far_away_match_list))

    rmsd,n=rv.get_values(id='forward')
    if n:
      print >>out,\
          "\nResidues matching in forward direction:   %4d  RMSD: %6.2f" %(
         n,rmsd)
      if verbose:
        for i,j in forward_match_list:
          print >>out,"ID:%d:%d  RESIDUES:  \n%s\n%s" %( i,j, chain_ca_lines[i],
           target_xyz_lines[j])

    rmsd,n=rv.get_values(id='reverse')
    if n:
      print >>out,\
         "Residues matching in reverse direction:   %4d  RMSD: %6.2f" %(
         n,rmsd)
      if verbose:
        for i,j in reverse_match_list:
          print >>out,"ID:%d:%d  RESIDUES:  \n%s\n%s" %(
           i,j, chain_ca_lines[i],
           target_xyz_lines[j])

    rmsd,n=rv.get_values(id='unaligned')
    if n:
      print >>out,\
        "Residues near but not matching one-to-one:%4d  RMSD: %6.2f" %(
         n,rmsd)
      if verbose:
        for i,j in unaligned_match_list:
          print >>out,"ID:%d:%d  RESIDUES:  \n%s\n%s" %(i,j, chain_ca_lines[i],
            target_xyz_lines[j])

    rmsd,n=rv.get_values(id='close')
    if n:
      lines_chain_ca=[]
      lines_target_xyz=[]
      for i,j in close_match_list:
        lines_chain_ca.append(chain_ca_lines[i])
        lines_target_xyz.append(target_xyz_lines[j])
      seq_chain_ca=get_seq_from_lines(lines_chain_ca)
      seq_target_xyz=get_seq_from_lines(lines_target_xyz)
      target_chain_ids=get_chains_from_lines(lines_target_xyz)
      target_length=get_target_length(target_chain_ids=target_chain_ids,
        hierarchy=target_ca,
        target_length_from_matching_chains=target_length_from_matching_chains)
      rv.add_target_length(id='close',target_length=target_length)

      if verbose:
         print "SEQ1:",seq_chain_ca,len(lines_chain_ca)
         print "SEQ2:",seq_target_xyz,len(lines_target_xyz)

      match_n,match_percent=get_match_percent(seq_chain_ca,seq_target_xyz)
      rv.add_match_percent(id='close',match_percent=match_percent)

      percent_close=rv.get_close_to_target_percent('close')
      if ratio_unique_to_total_target:
        percent_close=percent_close/ratio_unique_to_total_target


      print >>out,\
        "\nAll residues near target: "+\
         "%4d  RMSD: %6.2f Seq match (%%):%5.1f  %% Found: %5.1f" %(
         n,rmsd,match_percent,percent_close)
      if verbose:
        for i,j in close_match_list:
          print >>out,"ID:%d:%d  RESIDUES:  \n%s\n%s" %(i,j, chain_ca_lines[i],
            target_xyz_lines[j])

    rmsd,n=rv.get_values(id='far_away')
    if n:
      print >>out,\
        "Residues far from target: %4d  RMSD: %6.2f" %(
         n,rmsd)
      if verbose:
        for i in far_away_match_list:
          print >>out,"ID:%d  RESIDUES:  \n%s" %(i,chain_ca_lines[i])

  rv.n_forward=n_forward
  rv.n_reverse=n_reverse
  rv.n=len(pair_list)
  rv.max_dist=params.comparison.max_dist
  rv.ratio_unique_to_total_target=ratio_unique_to_total_target
  rv.total_target=total_target
  rv.total_chain=total_chain
  rv.used_target=used_target
  rv.used_chain=used_chain
  return rv

if __name__=="__main__":
  args=sys.argv[1:]
  run(args=args,out=sys.stdout)
