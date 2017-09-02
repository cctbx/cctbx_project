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

    unique_only = True
      .type = bool
      .help = Use only unique chains in query
      .short_caption = Unique only

    unique_target_pdb_in = None
      .type = path
      .help = Target model to choose which element is selected with \
           unique_only. NOTE: must be specified by keyword.
      .short_caption = Target model

    query_dir = None
      .type = path
      .help = directory containing query PDB files (any number)
      .short_caption = Query directory
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
      .caption = Maximum distance between atoms to be considered close

    distance_per_site = None
      .type = float
      .short_caption = Maximum distance spanned by a pair of residues
      .help =Maximum distance spanned by a pair of residues.  Set by \
            default as 3.8 A for protein and 8 A for RNA

    target_length_from_matching_chains = False
      .type = bool
      .short_caption = Use matching chains to get length
      .caption = Use length of chains in target that are matched to \
                define full target length (as opposed to all unique \
                chains in target).
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
      return 100.*n/max(1.,target_length)

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

def apply_atom_selection(atom_selection,hierarchy=None):
  asc=hierarchy.atom_selection_cache()
  sel = asc.selection(string = atom_selection)
  return hierarchy.select(sel)

def select_atom_lines(hierarchy):
  lines=[]
  for atom in hierarchy.atoms():
    lines.append(atom.format_atom_record())
  return lines

def get_best_match(xyz1,xyz2,crystal_symmetry=None,
    distance_per_site=None):
  if crystal_symmetry:
    assert distance_per_site is not None
    return best_match(
      xyz1,xyz2,
      crystal_symmetry=crystal_symmetry,
      distance_per_site=distance_per_site)
  else: # do it without symmetry
    (distance,i,j)=xyz1.min_distance_between_any_pair_with_id(xyz2)
    from libtbx import group_args
    return group_args(i=i,j=j,distance=distance)

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
    if a==b: match_n+=1
  match_percent=100.*match_n/len(seq1)
  return match_n,match_percent

def extract_unique_part_of_hierarchy(ph,target_ph=None,out=sys.stdout):
  new_hierarchy=iotbx.pdb.input(
    source_info="Model",lines=flex.split_lines("")).construct_hierarchy()
  mm=iotbx.pdb.hierarchy.model()
  new_hierarchy.append_model(mm)

  if target_ph:
    target_centroid_list=flex.vec3_double()
    for model in target_ph.models()[:1]:
      target_centroid_list.append(model.atoms().extract_xyz().mean())
  else:
    target_centroid_list=None
  unique_sequences=[]
  best_chain_dict={}
  best_chain_dist_dict={}
  for model in ph.models()[:1]:
    for chain in model.chains():
      try:
        seq=chain.as_padded_sequence()  # has XXX for missing residues
      except Exception, e:
        seq="XXX"
      if not seq in unique_sequences:
        unique_sequences.append(seq)
      if target_centroid_list:
        xx=flex.vec3_double()
        xx.append(chain.atoms().extract_xyz().mean())
        dist=xx.min_distance_between_any_pair(target_centroid_list)
      else:
        dist=0.
      best_dist=best_chain_dist_dict.get(seq)
      if best_dist is None or dist<best_dist:
        best_chain_dist_dict[seq]=dist
        best_chain_dict[seq]=chain
 

  for seq in best_chain_dist_dict.keys():
    chain=best_chain_dict[seq]
    if not chain:
      print >>out,"Mising chain for sequence %s" %(seq)
    else:
      mm.append_chain(chain.detached_copy())
      print >>out, "Adding chain %s: %s (%s): %7.2f" %(
         chain.id,seq.replace("X",""),str(chain.atoms().extract_xyz()[0]),
         best_chain_dist_dict[seq])


  return new_hierarchy

def run_all(params=None,out=sys.stdout):
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
      rv=run(params=local_params,out=local_out)
    except Exception,e:
      if str(e).find("CifParserError"):
        print >>out,"NOTE: skipping %s as it is not a valid model file" %(
           file_name)
        continue # it was not a valid PDB file...just skip it
      else:
        raise Sorry(str(e))
    rv_list.append(rv)
    file_list.append(file_name)

  print >>out,"\nCLOSE is within %4.1f A. FAR is greater than this." %(
    params.comparison.max_dist)
  print >>out,"\nCA SCORE is fraction in close CA / rmsd of these CA."
  print >>out,"\nSEQ SCORE is fraction (close and matching target sequence).\n"

  print >>out,"\n"
  print >>out,"               ----ALL RESIDUES----     CLOSE RESIDUES ONLY    %"
  print >>out,\
              "     MODEL     --CLOSE-    ---FAR--    FORWARD REVERSE MIXED"+\
              " FOUND   CA                   SEQ"
  print >>out,"               RMSD   N    RMSD   N       N       N      N  "+\
              "        SCORE  SEQ MATCH(%)  SCORE"+"\n"

  results_dict={}
  score_list=[]
  for rv,full_f in zip(rv_list,file_list):
    results_dict[full_f]=rv
    (rmsd,n)=rv.get_values('close')
    target_length=rv.get_target_length('close')
    score=n/(max(1,target_length)*max(0.1,rmsd))
    score_list.append([score,full_f])
    print "ZZ",full_f,rmsd,n,target_length,score
  score_list.sort()
  score_list.reverse()
  for score,full_f in score_list:
    rv=results_dict[full_f]
    percent_close=rv.get_close_to_target_percent('close')
    seq_score=rv.get_match_percent('close')*percent_close/10000
    file_name=os.path.split(full_f)[-1]
    close_rmsd,close_n=rv.get_values('close')
    far_away_rmsd,far_away_n=rv.get_values('far_away')
    forward_rmsd,forward_n=rv.get_values('forward')
    reverse_rmsd,reverse_n=rv.get_values('reverse')
    unaligned_rmsd,unaligned_n=rv.get_values('unaligned')
    match_percent=rv.get_match_percent('close')
    print >>out,"%14s %4.2f %4d   %4.1f %4d   %4d    %4d    %4d  %5.1f %6.2f   %5.1f      %6.2f" %(file_name,close_rmsd,close_n,far_away_rmsd,far_away_n,forward_n,
         reverse_n,unaligned_n,percent_close,score,match_percent,seq_score)

def get_target_length(target_chain_ids=None,hierarchy=None,
     target_length_from_matching_chains=None):
  total_length=0  # just counts residues
  for model in hierarchy.models()[:1]:
    for chain in model.chains():
      if (not target_length_from_matching_chains) or chain.id in target_chain_ids:
        total_length+=len(chain.residues())
  return total_length

def run(args=None,
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
  if params.input_files.unique_target_pdb_in and params.input_files.unique_only:
    print >>out,"Using %s as target for unique chains" %(
       params.input_files.unique_target_pdb_in)
  if params.input_files.query_dir and \
      os.path.isdir(params.input_files.query_dir):
    print >>out,"\nUsing all files in %s as queries\n" %(
       params.input_files.query_dir)
    return run_all(params=params,out=out)

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
    pdb_inp=get_pdb_inp(file_name=chain_file  )
    if params.input_files.unique_target_pdb_in:
      target_unique_hierarchy=get_pdb_inp(
        file_name=params.input_files.unique_target_pdb_in).construct_hierarchy()
    else:
      target_unique_hierarchy=None
    if not crystal_symmetry:
      crystal_symmetry=pdb_inp.crystal_symmetry_from_cryst1()
    chain_hierarchy=pdb_inp.construct_hierarchy()
    if params.input_files.unique_only:
      print >>out,"\nUsing only unique part of query\n"
      chain_hierarchy=extract_unique_part_of_hierarchy(
        chain_hierarchy,target_ph=target_unique_hierarchy,out=local_out)
    target_pdb_inp=get_pdb_inp(file_name=target_file)
    if not crystal_symmetry or not crystal_symmetry.unit_cell():
      crystal_symmetry=target_pdb_inp.crystal_symmetry_from_cryst1()
    target_hierarchy=target_pdb_inp.construct_hierarchy()

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
  for i in xrange(chain_xyz_fract.size()):
    best_j=None
    best_dd=None
    distance=None
    if working_crystal_symmetry:
      info=get_best_match(
        flex.vec3_double([chain_xyz_fract[i]]),target_xyz_fract,
        crystal_symmetry=working_crystal_symmetry,
        distance_per_site=distance_per_site)
      if info:
        distance=info.dist()
    else:
      info=get_best_match(
        flex.vec3_double([chain_xyz_cart[i]]),target_xyz_cart)
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
    pair_list.append([i,best_j,best_dd])
  n_forward=0
  n_reverse=0
  forward_match_list=[]
  reverse_match_list=[]
  forward_match_rmsd_list=flex.double()
  reverse_match_rmsd_list=flex.double()
  unaligned_match_list=[]
  unaligned_match_rmsd_list=flex.double()
  close_match_rmsd_list=flex.double()
  close_match_list=[]
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
  if forward_match_rmsd_list.size():
      id='forward'
      rmsd=forward_match_rmsd_list.min_max_mean().mean**0.5
      n=forward_match_rmsd_list.size()
      rv.add_rmsd(id=id,rmsd=rmsd,n=n)
  if reverse_match_rmsd_list.size():
      id='reverse'
      rmsd=reverse_match_rmsd_list.min_max_mean().mean**0.5
      n=reverse_match_rmsd_list.size()
      rv.add_rmsd(id=id,rmsd=rmsd,n=n)
  if unaligned_match_rmsd_list.size():
      id='unaligned'
      rmsd=unaligned_match_rmsd_list.min_max_mean().mean**0.5
      n=unaligned_match_rmsd_list.size()
      rv.add_rmsd(id=id,rmsd=rmsd,n=n)
  if close_match_rmsd_list.size():
      id='close'
      rmsd=close_match_rmsd_list.min_max_mean().mean**0.5
      n=close_match_rmsd_list.size()
      rv.add_rmsd(id=id,rmsd=rmsd,n=n)

  if far_away_match_rmsd_list.size():
      id='far_away'
      rmsd=far_away_match_rmsd_list.min_max_mean().mean**0.5
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
  return rv

if __name__=="__main__":
  args=sys.argv[1:]
  run(args=args,out=sys.stdout)
