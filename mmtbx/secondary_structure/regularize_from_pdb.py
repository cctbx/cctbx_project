from __future__ import absolute_import, division, print_function

# regularize_from_pdb.py
# a tool to replace segments of a structure with similar segment with good
# geometry from the PDB

import math
from iotbx.pdb import resseq_encode
import iotbx.phil
import os,sys
from libtbx.utils import Sorry,null_out
from scitbx.matrix import col
from scitbx.array_family import flex
from mmtbx.secondary_structure.find_ss_from_ca import \
   find_secondary_structure, \
   find_helix,find_beta_strand,find_other_structure,helix,strand,other,\
   get_first_resno,get_last_resno,get_sequence,get_chain_id,get_atom_list,\
   apply_atom_selection,model_info,split_model,merge_hierarchies_from_models, \
   get_pdb_hierarchy
from six.moves import zip
from six.moves import range

master_phil = iotbx.phil.parse("""

  input_files {
    map_coeffs_file = None
      .type = path
      .help = File with map coefficients
      .short_caption = Map coefficients
      .style = bold file_type:hkl input_file process_hkl child:fobs:data_labels\
        child:space_group:space_group child:unit_cell:unit_cell

    map_coeffs_labels = None
      .type = str
      .input_size = 160
      .help = Optional label specifying which columns of of map coefficients \
          to use
      .short_caption = Map coeffs label
      .style = bold renderer:draw_fobs_label_widget

    seq_file = None
      .type = path
      .help = Sequence file
      .short_caption = Sequence file

    pdb_in = None
      .type = path
      .help = Input PDB file
      .short_caption = Input PDB file

    secondary_structure_input = None
      .type = bool
      .help = Not used
      .style = hidden

    force_secondary_structure_input = None
      .type = bool
      .help = Not used
      .style = hidden

  }
  output_files {

    pdb_out = placed.pdb
      .type = path
      .help = Output PDB file with placed segments only
      .short_caption = Output PDB file with placed segments only

    pdb_records_file = None
      .type = path
      .help = Not used
      .style = hidden

  }
  crystal_info {
     resolution = None
       .type = float
       .help = High-resolution limit. Data will be truncated at this\
               resolution. If a map is supplied, it will be Fourier \
               filtered at this resolution. Required if input is a map and \
                only_original_map is not set.
       .short_caption = High-resolution limit
       .style = resolution
     space_group = None
       .type = space_group
       .short_caption = Space Group
       .help = Space group (normally read from the data file)
     unit_cell = None
       .type = unit_cell
       .short_caption = Unit Cell
       .help = Unit Cell (normally read from the data file)

     sequence = None
       .type = str
       .short_caption = Sequence
       .help = Sequence as text string. Normally supply a sequence file instead

  }
  find_ss_structure { # Note overwrites values in find_ss_structure.py

     ss_by_chain = None
       .type = bool
       .help = Find secondary structure only within individual chains. \
               Alternative is to allow H-bonds between chains. Can be \
               much slower with ss_by_chain=False. If your model is complete \
               use ss_by_chain=True. If your model is many fragments, use \
               ss_by_chain=False.  Not used in regularize_pdb
       .short_caption = Secondary structure by chain
       .style = hidden

     auto_choose_ca_vs_ca_n_o = True
       .type = bool
       .help = Automatically identify whether chains are mostly CA or mostly \
                contain CA/N/O atoms (requires min_ca_n_o_completeness).
       .short_caption = Auto choose CA vs CA/N/O

     min_ca_n_o_completeness = 0.95
       .type = float
       .help = Minimum completeness of CA/N/O atoms to not use CA-only
       .short_caption = Minimum completeness of CA/N/O

     max_rmsd = 1
       .type = float
       .help = Maximum rmsd to consider two chains with identical sequences \
               as the same for ss identification. \
               Not used in regularize_pdb
       .short_caption = Maximum rmsd
       .style = hidden
     use_representative_chains = True
       .type = bool
       .help = Use a representative of all chains with the same sequence. \
               Alternative is to examine each chain individually. Can be \
               much slower with use_representative_of_chain=False if there \
               are many symmetry copies. Ignored unless ss_by_chain is True. \
               Not used in regularize_pdb
       .style = hidden
     max_representative_chains = 100
       .type = float
       .help = Maximum number of representative chains\
               Not used in regularize_pdb
       .short_caption = Maximum representative chains
       .style = hidden

     find_alpha = True
       .type = bool
       .help = Find alpha helices
       .short_caption = Find alpha helices

     helices_are_alpha = True
       .type = bool
       .help = Find alpha helices and not three_ten or pi
       .short_caption = Helices are alpha

     find_three_ten = False
       .type = bool
       .help = Find three_ten helices
       .short_caption = Find three_ten helices

     find_pi = False
       .type = bool
       .help = Find pi helices
       .short_caption = Find pi helices

     find_beta = True
       .type = bool
       .help = Find beta structure
       .short_caption = Find beta structure

     find_other = True
       .type = bool
       .help = Find other structure
       .short_caption = Find other structure

     exclude_alpha_in_beta  = False
       .type = bool
       .help = Exclude regions already identified as alpha from three_ten, pi,\
               and beta
       .short_caption = Exclude alpha from beta

     make_unique = False
       .type = bool
       .help = Assign each residue to a unique type of structure
       .short_caption = Assign residues to unique structure

     cut_up_segments = True
       .type = bool
       .help = Cut up segments (make short segments of secondary structure)
       .short_caption = Cut up segments

     extend_segments = False
       .type = bool
       .help = Try to extend segments in both directions one residue at a time
       .short_caption = Extend segments

     write_helix_sheet_records = False
       .type = bool
       .help = Write HELIX and SHEET records
       .short_caption = Write HELIX/SHEET records

     set_up_helices_sheets = False
       .type = bool
       .help = Set up HELIX and SHEET records
       .short_caption = Set up HELIX/SHEET records

     include_single_strands = False
       .type = bool
       .help = Write SHEET records that contain a single strand
       .short_caption = Write single strands

     remove_missing_atom_annotation = False
       .type = bool
       .help = Remove annotation that refers to atoms that are not present
       .short_caption = Remove missing atom annotation

     max_h_bond_length = 3.5
       .type = float
       .help = Maximum H-bond length to include in secondary structure
       .short_caption = Maximum H-bond length

    search_secondary_structure = True
      .type = bool
      .help = Search for secondary structure in input model. \
              (Alternative is to just use secondary structure from \
              secondary_structure_input.)
      .short_caption = Find secondary structure

    combine_annotations = True
      .type = bool
      .help = Combine annotations if an input annotation is provided
      .short_caption = Combine annotations

    require_h_bonds = False
      .type = bool
      .help = Remove all secondary structure records that have fewer than \
              minimum_h_bonds good hydrogen bonds
      .short_caption = Require H-bonds

    minimum_h_bonds = 1
      .type = int
      .help = Minimum number of good hydrogen bonds to keep secondary \
              structure element (helix/sheet) if require_h_bonds is set
      .short_caption = Minimum number of H bonds

    maximum_poor_h_bonds = None
      .type = int
      .help = Maximum number of poor hydrogen bonds to keep secondary \
              structure element (helix/sheet) if require_h_bonds is set. \
              Note: None means ignore this test, 0 means allow no poor H-bonds.
      .short_caption = Maximum number of poor H bonds

  }

  extract_segments_from_pdb {
    extract = *None alpha beta other
       .type = choice
       .help = Extract all segments (alpha/beta) from a PDB file. \
               Used to create libraries of segments
       .short_caption = Extract segments from PDB
  }
  alpha {
     include scope mmtbx.secondary_structure.secondary_structure_params.alpha_params

     library = 'helices.pdb'
       .type = path
       .help = Helix library
       .short_caption = Helix library

     replace = True
       .type = bool
       .help = Replace alpha helices
       .short_caption = Replace alpha helices

     maximum_rmsd = 1.5
       .type = float
       .help = Maximum rmsd of CA atoms to template for helices
       .short_caption = Maximum rmsd of CA atoms to template for helices

     index_delta_x = None
       .type = float
       .help = Index gridding in x. \
                 Segments are indexed by 1-N distance and 1-N/2 distance in \
                 bins of index_delta_x and index_delta_y .  \
                 1 A is usually about right for both.
       .short_caption = Index gridding in x

     index_delta_y = None
       .type = float
       .help = Index gridding in y. \
                 Segments are indexed by 1-N distance and 1-N/2 distance in \
                 bins of index_delta_x and index_delta_y .  \
                 1 A is usually about right for both.

  }
  beta {
     include scope mmtbx.secondary_structure.secondary_structure_params.beta_params

     library = 'strands.pdb'
       .type = path
       .help = Strand library
       .short_caption = Strand library

     replace = True
       .type = bool
       .help = Replace beta structure
       .short_caption = Replace beta structure

     maximum_rmsd = 1.5
       .type = float
       .help = Maximum rmsd of CA atoms to template for beta
       .short_caption = Maximum rmsd of CA atoms to template for beta

     index_delta_x = None
       .type = float
       .help = Index gridding in x. \
                 Segments are indexed by 1-N distance and 1-N/2 distance in \
                 bins of index_delta_x and index_delta_y .  \
                 1 A is usually about right for both.
       .short_caption = Index gridding in x

     index_delta_y = None
       .type = float
       .help = Index gridding in y. \
                 Segments are indexed by 1-N distance and 1-N/2 distance in \
                 bins of index_delta_x and index_delta_y .  \
                 1 A is usually about right for both.

  }
  other {
     include scope mmtbx.secondary_structure.secondary_structure_params.other_params

     library = 'lib8.pdb'
       .type = path
       .help = Other library
       .short_caption = Other library

     replace = True
       .type = bool
       .help = Replace other structure
       .short_caption = Replace other structure

     maximum_rmsd = 1.5
       .type = float
       .help = Maximum rmsd of CA atoms to template for other
       .short_caption = Maximum rmsd of CA atoms to template for other

     index_delta_x = 1.
       .type = float
       .help = Index gridding in x. \
                 Segments are indexed by 1-N distance and 1-N/2 distance in \
                 bins of index_delta_x and index_delta_y .  \
                 1 A is usually about right for both.
       .short_caption = Index gridding in x

     index_delta_y = 1.
       .type = float
       .help = Index gridding in y. \
                 Segments are indexed by 1-N distance and 1-N/2 distance in \
                 bins of index_delta_x and index_delta_y .  \
                 1 A is usually about right for both.
       .short_caption = Index gridding in y

  }
  regularization {

     include_side_chain = False
       .type = bool
       .help = include side chain (beyond CB) in replacement
       .short_caption = include side chains beyond CB

     maximum_overlap = 2
       .type = int
       .help = Maximum overlap of fragments in assembly
       .short_caption = Maximum overlap

     maximum_junction_rmsd = 1.0
       .type = float
       .help = Maximum rmsd of main-chain atoms at residues where fragments \
                are joined
       .short_caption = Maximum rmsd at junctions

     good_enough_ratio = 0.5
       .type = float
       .help = Accept a segment if rmsd to target is good_enough_ratio times \
            the maximum_rmsd or better (and do not look further). Note that \
            segments are sorted on frequency of occurrence so those at the \
            top of the list are more likely to occur than those at the end.
       .short_caption = Good enough ratio for maximum_rmsd

  }
  control {
      resolve_size = None
        .type = str
        .help = "Size of resolve to use. "
        .short_caption = Size of RESOLVE to use
      verbose = False
        .type = bool
        .help = Verbose output
        .short_caption = Verbose output
  }
""", process_includes=True)
master_params = master_phil

def get_and_split_model(pdb_hierarchy=None,
     pdb_in=None,get_start_end_length=None,out=sys.stdout):
    if not pdb_hierarchy:
      print("Reading model from %s" %(pdb_in), file=out)
      if not pdb_in:
        return []
      if not os.path.isfile (pdb_in):
        raise Sorry("Missing the file %s" %(pdb_in))
      pdb_hierarchy=get_pdb_hierarchy(text=open(pdb_in).read())
    model=model_info(
      hierarchy=pdb_hierarchy,
      info={})
    models=split_model(model=model,verbose=False)
    print("Split model into %d chains" %(len(models)), file=out)

    for model in models:
      model.hierarchy.remove_alt_confs(always_keep_one_conformer=False)
      if get_start_end_length:
        model.info['model_start_resno']=get_first_resno(model.hierarchy)
        model.info['model_end_resno']=get_last_resno(model.hierarchy)
        id=model.info['chain_number']
        ll=model.info['model_end_resno']-model.info['model_start_resno']+1
        model.info['length']=ll
    return models

class segment_library:

  def __init__(self,params=None,
     segment_class=None,out=sys.stdout):

    self.segment_class=segment_class # strand helix other
    self.index_length=None

    # get the segments as model_info objects
    if params.library is None:
      params.library=""
    elif not os.path.isfile(params.library):
      import libtbx.load_env
      name=os.path.join('cctbx_project','mmtbx',
         'secondary_structure','regularize_from_pdb_lib',params.library)
      full_name=libtbx.env.find_in_repositories(
         relative_path=name, test=os.path.exists)
      if not full_name:
        raise Sorry("Cannot find the library file %s" %(params.library))
      params.library=full_name

    self.models=get_and_split_model(pdb_in=params.library,out=out)

    # get information about the segments, if available
    info_file=params.library.replace(".pdb",".info")
    self.get_info_file(info_file=info_file,out=out)

    # set up the segments with information if available
    self.set_up_segments(params=params,out=out)

    # set up indexing so that we can find segments quickly
    if self.have_info:
      self.set_up_indexing()

  def get_info_file(self,info_file=None,out=sys.stdout):
    self.max_count=None
    if os.path.isfile(info_file):
      print("Reading other library info from %s" %(info_file), file=out)
      self.model_number_list=[]
      self.model_n_obs=[]
      self.model_rms=[]
      info_number_of_segments=None
      for line in open(info_file).readlines():
        spl=line.split()
        if info_number_of_segments is None:
          info_number_of_segments=int(spl[-1])
        else:
          self.model_number_list.append(int(spl[1]))
          nn=int(spl[4])
          self.model_n_obs.append(int(spl[4]))
          self.model_rms.append(float(spl[-1]))
          if self.max_count is None or nn>self.max_count:
            self.max_count=nn
      if len(self.model_number_list)!=info_number_of_segments:
       print("Model numbers: %d  Segments: %d " %(
         len(self.model_number_list),info_number_of_segments))
      assert len(self.model_number_list)==info_number_of_segments
      assert len(self.model_number_list)==len(self.models)
      self.have_info=True
    else:
      self.model_number_list=[None]*len(self.models)
      self.model_n_obs=[None]*len(self.models)
      self.model_rms=[None]*len(self.models)
      self.have_info=False

  def set_up_segments(self,params=None,out=sys.stdout):
    self.segments=[]
    frequency=None
    for model,count in zip(self.models,self.model_n_obs):
      if self.max_count:
        frequency=count/self.max_count
      self.segments.append(
        self.segment_class(
          params=params,hierarchy=model.hierarchy,frequency=frequency,out=out))

  def set_up_indexing(self):
    # set up indexing to find correct ones quickly
    index_dict={}
    for segment in self.segments:
      if self.index_length is None:
        self.index_length=segment.standard_length
      assert self.index_length==segment.standard_length
      index_list=self.get_index_list(segment=segment)
      for index in index_list:
        if not index in index_dict:
          index_dict[index]=[]
        index_dict[index].append(segment)
    self.index_dict=index_dict

  def get_index_from_x_y(self,x,y,delta_x=None,delta_y=None):
     index_x=int(0.5+x/delta_x)
     index_y=int(0.5+y/delta_y)
     return "%d:%d" %(index_x,index_y)

  def get_index_list(self,segment=None,sites=None,
      delta_x=1.,delta_y=1.,only_best=False):
     # get list of indices for this segment.  They correspond to the
     # values of dist(1,n) and dist(1,n/2) in this segment
     if sites is None:
       sites=segment.get_sites()
     if self.index_length!=len(sites):
       return None # the indexing only applies if we have the same length segment
     n=len(sites)-1
     n2=n//2
     dist1n=col(sites[n])-col(sites[0])
     dist1n=dist1n.length()
     dist1n2=col(sites[n2])-col(sites[0])
     dist1n2=dist1n2.length()
     if only_best:
       return self.get_index_from_x_y(
         dist1n,dist1n2,delta_x=delta_x,delta_y=delta_y)

     # Here when we are setting up the indices (only_best=False)
     index_list=[]
     for i in [-delta_x,0,delta_x]:
       for j in [-delta_y,0,delta_y]:
         index_list.append(
           self.get_index_from_x_y(
             dist1n+i,dist1n2+j,delta_x=delta_x,delta_y=delta_y))
     return index_list

  def get_xyz_in_order(self,xyz=None,original_order=None,
        target_order=None):
    #put xyz in the order required to match order to target_order
    if original_order==target_order:
      return xyz  # already was ok

    recover_target_order=range(len(target_order))
    recover_sort_list=[]
    for target,recover in zip(target_order,recover_target_order):
      recover_sort_list.append([target,recover])
    recover_sort_list.sort() # sorted on target_order, tag is place these go

    sort_list=[]
    for x,orig in zip(xyz,original_order):
      sort_list.append([orig,x])
    sort_list.sort()  # sorted on original order, tag is xyz


    second_sort_list=[]
    for [orig,x],[target,recover] in zip(sort_list,recover_sort_list):
      second_sort_list.append([recover,x,orig,target])
    second_sort_list.sort()

    new_xyz=flex.vec3_double()
    new_order=[]
    for [recover,x,orig,target] in second_sort_list:
      new_xyz.append(x)
      new_order.append(orig)

    assert new_order==target_order

    return new_xyz

  def get_main_chain_rmsd(self,lsq_fit=None,before_junction=None,
            start_res=None,end_res=None,  # in other_segment
            model_segment=None,other_segment=None):
    # find matching residues in model_segment and other_segment that have
    #  low rms for main_chain rmsd. Must have residues on up to before_junction
    #   and higher than before_junction
    # how to select matching residues...
    first_res_in_model_segment=model_segment.get_start_resno()
    last_res_in_model_segment=model_segment.get_end_resno()
    assert last_res_in_model_segment-first_res_in_model_segment==end_res-start_res
    placed_other_hierarchy=model_segment.apply_lsq_fit(
         lsq_fit_obj=lsq_fit,
          hierarchy=other_segment.hierarchy,
          start_res=start_res,
          end_res=end_res)

    rms_list=[]
    model_resseq_list=[]
    other_resseq_list=[]
    for i,j in zip(
       range(first_res_in_model_segment,last_res_in_model_segment+1),
       range(start_res,end_res+1)):
      atom_selection=\
         "resid %s through %s and (name n or name c or name ca or name o)" %(
          resseq_encode(i),resseq_encode(i))
      model_res_hierarchy=apply_atom_selection(
        atom_selection,hierarchy=model_segment.hierarchy)
      atom_selection=\
         "resid %s through %s and (name n or name c or name ca or name o)" %(
          resseq_encode(j),resseq_encode(j))
      other_res_hierarchy=apply_atom_selection(
        atom_selection,hierarchy=placed_other_hierarchy)
      xyz1=model_res_hierarchy.atoms().extract_xyz()
      xyz2=other_res_hierarchy.atoms().extract_xyz()
      assert xyz1.size()==xyz2.size()
      # make sure xyz2 is in same order as xyz1 based
      atom_list_1=get_atom_list(model_res_hierarchy)
      atom_list_2=get_atom_list(other_res_hierarchy)
      xyz2_sorted=self.get_xyz_in_order(xyz=xyz2,original_order=atom_list_2,
        target_order=atom_list_1)
      diffs=xyz1-xyz2_sorted
      rms_list.append(diffs.rms_length())
      model_resseq_list.append(i)
      other_resseq_list.append(j)

    # get rms by residue pair
    best_left_rms=None
    best_right_rms=None
    best_left_crossover=None
    best_right_crossover=None
    best_left_crossover_other=None
    best_right_crossover_other=None
    for rms,model_resseq,other_resseq in zip(
      rms_list,model_resseq_list,other_resseq_list):
      if model_resseq<=before_junction:
        if best_left_rms is None or best_left_rms>rms:
          best_left_rms=rms
          best_left_crossover=model_resseq
          best_left_crossover_other=other_resseq
      else:
        if best_right_rms is None or best_right_rms>rms:
          best_right_rms=rms
          best_right_crossover=model_resseq
          best_right_crossover_other=other_resseq

    return best_left_rms,best_right_rms,\
       best_left_crossover,best_right_crossover,\
       best_left_crossover_other,best_right_crossover_other


  def get_replacement_segment(self,model_segment=None,
       maximum_rmsd=None,good_enough_rmsd=None,
       score_by_main_chain_rmsd=None,before_junction=None,
       start_position=None,
       delta_residues=0,
       verbose=None,out=sys.stdout):
    # find a segment in segment_lib that matches model_segment.  Choose
    #  n residues in segment from segment lib where n=length of model_segment
    # delta_residues is change in length of model_segment to insert from
    #   segment_lib
    n=model_segment.length()
    best_segment=None
    best_lsq_fit=None
    best_rms=None
    best_start_res=None
    best_end_res=None
    best_id=None
    best_left_crossover=None
    best_right_crossover=None

    if delta_residues is None:
      delta_residues=0

    placed_hierarchy=None
    all_too_short=True

    nn=0
    found=False
    segments_to_use=None
    if self.have_info and not score_by_main_chain_rmsd:
      index=self.get_index_list(model_segment,only_best=True)
      segments_to_use=self.index_dict.get(index)

    if segments_to_use is None:
      segments_to_use=self.segments

    for segment in segments_to_use:
      if found: break
      nn+=1
      m=segment.length()
      extra_residues=m-n-delta_residues
      if extra_residues<0:
        continue # skip if too short
      if start_position is None or start_position>extra_residues:
         range_to_use=range(extra_residues+1)
      else:  # just use start_position
         range_to_use=[start_position]
      all_too_short=False
      # get transformation to map segment on to model_segment
      for offset in range_to_use:
        if found: break
        start_res=offset+segment.get_start_resno() # start res in seg library
        end_res=start_res+n+delta_residues-1

        lsq_fit=model_segment.get_lsq_fit(other=segment,
          start_res=start_res,end_res=end_res)
        if not lsq_fit: continue
        if lsq_fit.rms > maximum_rmsd: continue
        rms=lsq_fit.rms

        if score_by_main_chain_rmsd: # see if we can find a pair of residues
        #  that overlap acceptably...
          left_rms,right_rms,left_crossover,right_crossover,\
             left_crossover_other,right_crossover_other=\
             self.get_main_chain_rmsd(lsq_fit=lsq_fit,
            before_junction=before_junction,
            start_res=start_res,end_res=end_res,
            model_segment=model_segment,other_segment=segment)
          rms=max(left_rms,right_rms)
        else:
          left_crossover,right_crossover,\
             left_crossover_other,right_crossover_other=None,None,None,None
        if best_rms is None or rms<best_rms:
            best_rms=rms
            best_segment=segment
            best_lsq_fit=lsq_fit
            best_start_res=start_res
            best_end_res=end_res
            best_id=nn
            best_left_crossover=left_crossover
            best_right_crossover=right_crossover
            best_left_crossover_other=left_crossover_other
            best_right_crossover_other=right_crossover_other
            if best_rms<=good_enough_rmsd:
              found=True
    if all_too_short:
      if verbose:
        print("No templates long enough for segment: %s " %(
          model_segment.summary()), file=out)
    if best_segment:
      # apply rt to the hierarchy in segment
      placed_hierarchy=model_segment.apply_lsq_fit(lsq_fit_obj=best_lsq_fit,
          hierarchy=best_segment.hierarchy,
          start_res=best_start_res,end_res=best_end_res)
      frequency=best_segment.frequency
    else:
      return None,None,None,None,None,None,None
    if frequency is not None:
      if verbose:
        print("Used segment # %d (examined %d) with freq %5.3f and rmsd %5.2f" %(
          best_id,nn,frequency,best_rms))
      log_frequency=math.log(frequency)
    else:
      log_frequency=None
      if verbose:
        print("Used segment # %d (examined %d) and rmsd %5.2f" %(
          best_id,nn,best_rms))
    return placed_hierarchy,best_rms,log_frequency,\
          best_left_crossover,best_right_crossover,\
             best_left_crossover_other,best_right_crossover_other

class connected_group_segment:
  def __init__(self,segment=None,start_resno=None,end_resno=None):
    # start_resno is where we will start this segment in this connected group
    #  (can be from segment.info['target_start_resno'] to
    #      segment.info['target_end_resno'])
    self.segment=segment
    self.start_resno=start_resno
    self.end_resno=end_resno
    self.delta_length=\
       self.segment.info['length']-self.segment.info['target_length']
    self.first_resno_of_hierarchy=get_first_resno(self.segment.hierarchy)
    self.get_score()


  #  resseq numbers that identify where this segment is going to go relative
  #  to the sequence of the target chain

  def get_target_start_resno(self): # resseq in target corresponding to first
    # residue in the hierarchy for this segment
    return self.segment.info['target_start_resno']

  def get_target_end_resno(self): # resseq in target corresponding to last
    # residue in the hierarchy for this segment
    return self.segment.info['target_end_resno']

  def get_start_resno(self): # resseq in target we want to start at
    # This is the left end of the final sequence for this target
    return self.start_resno

  def get_end_resno(self): # resseq in target we want to end at
    # This is the right end of the final sequence for this target
    return self.end_resno

  def get_delta_length(self): # the difference between the length to be inserted
    # and the number of residues it will replace. Usually zero.
    return self.delta_length

  def get_length(self): # length of the part of the segment to be used.
    #  Does not include delta_length (it is the length of the segment that
    #  it would replace)
    return self.end_resno-self.start_resno+1

  #  offsets and resseq numbers that identify the residue numbers in the
  #  (arbitrary) hierarchy for this segment that we are using in this
  #  segment

  def get_first_resno_of_hierarchy(self): # first resseq in hierarchy
    # This is the first residue number in the hierarchy (the whole thing,
    #  not just the part we are going to use)
    return self.first_resno_of_hierarchy

  def get_start_offset(self):
    # how far in to the hierarchy in this cgs.segment we start:
    return self.start_resno-self.get_target_start_resno()

  def get_start_resno_of_hierarchy(self): # resseq in hierarchy to start at
    # this is the left end of the part of the hierarchy to use
    return self.get_first_resno_of_hierarchy()+self.get_start_offset()

  def get_end_resno_of_hierarchy(self): # resseq in hierarchy to end at
    # this is the right end of the part of the hierarchy to use
    # Note if there is an insertion then this may be a longer segment
    #  than the target
    return self.get_start_resno_of_hierarchy()+\
        self.get_length()-1+self.get_delta_length()

  def has_left_end(self):
    return self.segment.info.get('contains_left_end_of_chain')

  def has_right_end(self):
    return self.segment.info.get('contains_right_end_of_chain')

  def get_score(self,
     rms_weight=-0.2,
     use_frequency_score=True,
     freq_factor=0.02,
     length_factor=1.,
     rms_factor=1.):
    # For now, score + on length and - on length*rmsd from target structure
    length=self.end_resno-self.start_resno+self.delta_length
    # do not include end but do include extra residues inserted
    rms=self.segment.info['rms'] # deviation from target
    base_score=self.segment.info['base_score']
    if use_frequency_score:
      log_frequency=self.segment.info['log_frequency']
      if log_frequency is not None:
        self.score=base_score*length*(
          length_factor+rms_factor*math.exp(-rms)+freq_factor*log_frequency)
      else: # for helices etc use freq=1
        self.score=base_score*length*math.exp(-rms)
    else:
      self.score=base_score*length* (1+rms_weight*rms)   # old version
    return self.score

  def show_summary(self,out=None,return_text=True):
    if self.delta_length:
      extra_text="(%d)" %(self.end_resno+self.delta_length)
    else:
      extra_text="    "
    text="CS: %d:%d%s Hierarchy: %s:%s Score: %7.2f" %(
        self.start_resno,self.end_resno,extra_text,
        str(self.get_start_resno_of_hierarchy()),
        str(self.get_end_resno_of_hierarchy()),
        self.get_score())
    if return_text: return text
    if out is None: out=sys.stdout
    print(text, file=out)
    return ""


class connected_group:
  def __init__(self,segment=None,start_resno=None,end_resno=None):
    self.connected_group_segments=[]
    if segment:
      self.connected_group_segments.append(
       connected_group_segment(segment=segment,
       start_resno=start_resno,end_resno=end_resno))
    self.find_left_connection()
    self.find_right_connection()
    self.connection_rms_dict={}
    self.get_score()
    self.model_is_complete=None
    self.model_rms=None
    self.model_rms_n=None

  def customized_copy(self):
    new_copy=connected_group()
    for cgs in self.connected_group_segments:
      new_copy.connected_group_segments.append(cgs)
    new_copy.score=self.score
    new_copy.left_connection=self.left_connection
    new_copy.right_connection=self.right_connection
    return new_copy

  def as_model(self,renumber=True,is_left_end=None,is_right_end=None,
     offset_to_renumber_right_end_at_target=False,
     model_to_match=None,
     use_resname_from_model_to_match=None,
     verbose=None,
     out=sys.stdout):
    # if model_to_match is provided, use chain name and residue numbers
    # Also use residue name from that model if use_resname_from_model_to_match
    # Also keep only up to CB and discard CB if GLY

    models=[]
    diffs=flex.vec3_double()
    skipped=0
    sequence=[]
    sequences=[]  # list of sequences, one for each model
    if model_to_match:
      chain_id=get_chain_id(model_to_match.hierarchy)
    else:
      chain_id='A'
    overall_start_resno=None
    overall_delta_residues=0
    for cgs,is_first,is_last in zip(self.connected_group_segments,
      [True]+[False]*(len(self.connected_group_segments)-1),
      [False]*(len(self.connected_group_segments)-1)+[True],
       ):
      # 2015-04-09 cgs.start_resno and cgs.end_resno define what we use

      if is_first and is_left_end: # keep original start of this segment
        target_start_resno=cgs.segment.info['target_start_resno']
      else:
        target_start_resno=cgs.start_resno
      if overall_start_resno is None:
        overall_start_resno=target_start_resno

      if is_last and is_right_end: # keep original end of this segment
        target_end_resno=cgs.segment.info['target_end_resno']
      else:
        target_end_resno=cgs.end_resno

      # how far in to the hierarchy in this cgs.segment we start:
      start_offset=target_start_resno-cgs.segment.info['target_start_resno']

      # where we start (by resid) in the hierarchy
      start_resno=get_first_resno(cgs.segment.hierarchy)+start_offset
      length=target_end_resno-target_start_resno+1
      end_resno=start_resno+length-1+cgs.delta_length # insert may be different
      overall_delta_residues+=cgs.delta_length


      if not is_last:  # skip last residue as it is first of next one
        length-=1
        end_resno-=1
        target_end_resno-=1

      if length==0:
        continue

      atom_selection="resid %s through %s" %(resseq_encode(start_resno),
       resseq_encode(end_resno))
      h=apply_atom_selection(atom_selection,hierarchy=cgs.segment.hierarchy)

      if model_to_match: # get comparison and sequence
        atom_selection="name ca and resid %s through %s" %(
         resseq_encode(target_start_resno),
         resseq_encode(target_end_resno))
        original_ca=apply_atom_selection(atom_selection,
          hierarchy=model_to_match.hierarchy)
        original_xyz=original_ca.atoms().extract_xyz()
        original_sequence=get_sequence(original_ca)

        atom_selection="name ca "
        replacement_ca=apply_atom_selection(atom_selection,hierarchy=h)
        replacement_xyz=replacement_ca.atoms().extract_xyz()
        replacement_sequence=get_sequence(replacement_ca)

        if original_xyz.size()==replacement_xyz.size():
          diffs.extend(original_xyz-replacement_xyz)
          new_sequence=original_sequence
        else:
          skipped+=original_xyz.size()
          new_sequence=original_sequence+(['GLY']*max(0,
            len(replacement_sequence)-len(original_sequence))
            )[:len(replacement_sequence)]
        sequences.append(new_sequence)
        sequence+=new_sequence


      m=model_info(hierarchy=h,
       info={ 'target_start_residue':target_start_resno,
              'target_end_residue':target_end_resno,
              'delta_length':cgs.delta_length,
              'inserted_length':target_end_resno-target_start_resno+1+\
                 cgs.delta_length } )
      models.append(m)

    if diffs:
      print("RMSD (%d residues): %5.2f A" %(
          diffs.size(),diffs.rms_length()), end=' ', file=out)
      self.model_rms=diffs.rms_length()
      self.model_rms_n=diffs.size()
      if skipped:
        print("\n(Not including %d residues in insertions/deletions)" %(skipped), file=out)
      else: print(file=out)
      if verbose:
        print("Sequence: %s" %(" ".join(sequence)), file=out)


    if offset_to_renumber_right_end_at_target:
      overall_start_resno=overall_start_resno-overall_delta_residues
    new_model=merge_hierarchies_from_models(models=models,renumber=True,
      first_residue_number=overall_start_resno,sequences=sequences,
      chain_id=chain_id,trim_side_chains=True)
    return new_model

  def has_insertions_deletions(self):
    for cgs in self.connected_group_segments:
      if cgs.delta_length !=0:
         return True
    return False

  def overlaps_with(self,other=None):
    # either end of other overlaps with something between left and right
    #  means it overlaps
    if self.get_left_connection() <= other.get_left_connection() and \
       self.get_right_connection() >= other.get_left_connection():
      return True
    elif self.get_left_connection() <= other.get_right_connection() and \
       self.get_right_connection() >= other.get_right_connection():
      return True
    else:
      return False

  def contains(self,other=None):
    if self.get_left_connection() <= other.get_left_connection() and \
       self.get_right_connection() >= other.get_right_connection():
      return True
    else:
      return False

  def is_duplicate(self,other_list=[]):
    for other in other_list:
      if self.is_same_as(other=other):
        return True
    return False

  def is_same_as(self,other=None):
    if self.connected_group_segments==other.connected_group_segments:
      return True
    else:
      return False

  def shared_residue_range(self,other=None):
    start_res=max(self.get_left_connection(),other.get_left_connection())
    end_res=min(self.get_right_connection(),other.get_right_connection())
    if end_res>=start_res:
      return [start_res,end_res]
    else:
      return None

  def shared_segment(self,other=None):
    for cgs in self.connected_group_segments:
      for cgs1 in other.connected_group_segments:
        if cgs==cgs1: return cgs

  def get_cgs_containing_resno(self,resno):
    self_cgs=None
    for cgs in self.connected_group_segments:
      if resno>=cgs.start_resno and resno<=cgs.end_resno:
         self_cgs=cgs
         break
    return self_cgs

  def join_by_residue(self,other=None,
      shared_residue_range=shared_residue_range,
      maximum_junction_rmsd=None):
    # find rms for each possible crossover
    start_resno,end_resno=shared_residue_range
    best_resno=None
    best_rms=None

    for resno in range(start_resno,end_resno+1):
      self_cgs=self.get_cgs_containing_resno(resno)
      other_cgs=other.get_cgs_containing_resno(resno)
      rms=self.get_connection_rms(
        self_cgs,other_cgs,resno,allow_unit_length=True)
      if rms is not None and \
          rms < maximum_junction_rmsd and (best_rms is None or rms<best_rms):
        best_rms=rms
        best_resno=resno
    if best_resno is None:
      return  # failed to find a connection

    # figure out which one goes to the left
    if self.get_left_connection()<=other.get_left_connection():
      left_cg=self
      right_cg=other
    else:
      left_cg=other
      right_cg=self

    # make sure the new right one is not contained within the left:
    if left_cg.get_right_connection()>=right_cg.get_right_connection():
      return False

    # create new connected group starting with left left_cg, including part of
    #  the cgs that connects, and then the right one

    left_junction_cgs=left_cg.get_cgs_containing_resno(best_resno)
    right_junction_cgs=right_cg.get_cgs_containing_resno(best_resno)
    # modify these two so that together they cover the junction

    lcgs=connected_group_segment(segment=left_junction_cgs.segment,
        start_resno=left_junction_cgs.start_resno,end_resno=best_resno)
    lcgs.get_score()

    rcgs=connected_group_segment(segment=right_junction_cgs.segment,
        start_resno=best_resno,end_resno=right_junction_cgs.end_resno)
    rcgs.get_score()

    # now construct new connected group with left half of left_cg, the two
    # joining segments and the right half of right_cg
    new_cg=connected_group()

    # put in segments from the left
    for cgs in left_cg.connected_group_segments:
      if cgs==left_junction_cgs: break
      new_cg.connected_group_segments.append(cgs)
    # put in our new two segments making the connection
    if lcgs.get_length()>=2: # do not include segments of length 1
      new_cg.connected_group_segments.append(lcgs)
    if rcgs.get_length()>=2:
      new_cg.connected_group_segments.append(rcgs)

    # add on segments to the right
    found=False
    for cgs in right_cg.connected_group_segments:
      if cgs==right_junction_cgs:
        found=True
        continue
      if not found:
        continue
      new_cg.connected_group_segments.append(cgs)

    # make these this cg
    self.connected_group_segments=new_cg.connected_group_segments
    # sort and score
    self.sort_connected_group_segments()
    self.get_score()
    self.find_left_connection()
    self.find_right_connection()
    return True

  def get_left_right_cg(self,cg1,cg2):
    if cg1.get_left_connection()<=cg2.get_left_connection():
      left_cg=cg1
      right_cg=cg2
    else:
      left_cg=cg2
      right_cg=cg1
    return left_cg,right_cg

  def find_residue_range(self,trim_dict=None):
    # find longest segment that is not trimmed
    keys=list(trim_dict.keys())
    keys.sort()
    res_start=None
    res_end=None
    best_res_range=None
    for key in keys:
      if trim_dict[key]: # bad residue
        if res_start is None: # nothing to do
          pass
        else:  # end what we have
          if best_res_range is None or \
             res_end-res_start > best_res_range[1]-best_res_range[0]:
            best_res_range=[res_start,res_end]
            res_start=None
            res_end=None
      else:
        if res_start is not None: # continue
          res_end=key
        else: # start it
          res_start=key
          res_end=key
    if res_start is not None and (best_res_range is None or \
      res_end-res_start > best_res_range[1]-best_res_range[0]):
         best_res_range=[res_start,res_end]
    return best_res_range


  def trim_residues(self,trim_dict=None,min_length=2):
    # trim out all residues (ranges) that are marked in trim_dict
    # find longest segment that is not trimmed
    # if result is less than min_length, drop it
    res_range=self.find_residue_range(trim_dict)
    if res_range is None or res_range[1]-res_range[0]+1<min_length:
      return False# remove the connected group
    elif res_range[0]==self.get_left_connection() and \
        res_range[1]==self.get_right_connection():
      return True # fine as is
    else:
      return self.keep_residue_range(start_res=res_range[0],end_res=res_range[1])

  def keep_residue_range(self,start_res=None,end_res=None):
    # remove everything except for this range
    # Any segment less than min_length is to be removed (typically 1 long)

    # put in segments from the left starting with segment that contains
    # leftmost residue (start_res)
    first_segment_to_use=self.get_cgs_containing_resno(start_res)
    last_segment_to_use=self.get_cgs_containing_resno(end_res)

    new_cg=connected_group()
    if first_segment_to_use==last_segment_to_use:
      # just take part of this one
      cgs=connected_group_segment(segment=first_segment_to_use.segment,
              start_resno=start_res,end_resno=end_res)
      cgs.get_score()
      new_cg.connected_group_segments.append(cgs)

    else:
      started=False
      for cgs in self.connected_group_segments:
        if not started:
          if cgs==first_segment_to_use:
            started=True
            # put in the right part of this segment
            lcgs=connected_group_segment(segment=cgs.segment,
              start_resno=start_res,end_resno=cgs.end_resno)
            lcgs.get_score()
            if lcgs.get_length()>1: # never keep shorter than 2
              new_cg.connected_group_segments.append(lcgs)
          else: # do nothing
            pass
        else: # already started
          if cgs==last_segment_to_use:
            # put in the right part of this segment
            rcgs=connected_group_segment(segment=cgs.segment,
              start_resno=cgs.start_resno,end_resno=end_res)
            rcgs.get_score()
            if rcgs.get_length()>1: # never keep shorter than 2
              new_cg.connected_group_segments.append(rcgs)
            break # all done
          else:  # put it in
            new_cg.connected_group_segments.append(cgs)

    # make these this cg
    self.connected_group_segments=new_cg.connected_group_segments
    # sort and score
    self.sort_connected_group_segments()
    self.get_score()
    self.find_left_connection()
    self.find_right_connection()
    if len(self.connected_group_segments)>0:
      return True
    else:
      return False


  def join(self,other=None):
    shared_cgs=self.shared_segment(other=other)
    new_segment_list=[]
    if self.get_left_connection()<other.get_left_connection() and \
       self.get_right_connection()<other.get_right_connection():
      # join self on left
      left_cg=self
      right_cg=other
    elif self.get_left_connection()>other.get_left_connection() and \
       self.get_right_connection()>other.get_right_connection():
      # join on right
      left_cg=other
      right_cg=self
    else:  # skip (probably the same)
      return

    for cgs in left_cg.connected_group_segments:
      if cgs==shared_cgs: break
      new_segment_list.append(cgs)
    found=False
    for cgs in right_cg.connected_group_segments:
      if cgs==shared_cgs: found=True
      if found: new_segment_list.append(cgs)

    self.connected_group_segments=new_segment_list
    self.sort_connected_group_segments()
    self.get_score()
    self.find_left_connection()
    self.find_right_connection()
    return True

  def add_connected_group_to_right(self,cg):
    if self.get_left_connection()!=cg.get_right_connection() and \
       self.get_right_connection()!=cg.get_left_connection():
      raise Sorry("No connection between these connected groups:"+\
        self.show_summary()+":"+cg.show_summary())
    for cgs in cg.connected_group_segments:
      self.connected_group_segments.append(cgs)
    self.sort_connected_group_segments()
    self.get_score()
    self.find_left_connection()
    self.find_right_connection()

  def sort_connected_group_segments(self):
    sort_list=[]
    for cgs in self.connected_group_segments:
      sort_list.append([cgs.start_resno,cgs])
    sort_list.sort()
    self.connected_group_segments=[]
    for id,cgs in sort_list:
      self.connected_group_segments.append(cgs)

  def get_left_connection(self):
    return self.left_connection

  def get_right_connection(self):
    return self.right_connection

  def find_left_connection(self):
    lc=None
    for cg in self.connected_group_segments:
      if lc is None or cg.start_resno< lc:
        lc=cg.start_resno
    self.left_connection=lc

  def find_right_connection(self):
    lc=None
    for cg in self.connected_group_segments:
      if lc is None or cg.end_resno> lc:
        lc=cg.end_resno
    self.right_connection=lc

  def get_left_segment_hierarchy(self):
    return self.connected_group_segments[0].segment.hierarchy

  def get_right_segment_hierarchy(self):
    return self.connected_group_segments[-1].segment.hierarchy

  def is_complete(self):
    return self.model_is_complete

  def get_rms(self):
    return self.model_rms

  def get_rms_n(self):
    return self.model_rms_n

  def get_junction_rms(self):
    junction_rms=0.
    junction_rms_n=0.
    if len(self.connected_group_segments)>1:
      for s1,s2 in zip(
        self.connected_group_segments[:-1],self.connected_group_segments[1:]):
        junction_rms+=self.get_connection_rms(s1,s2)**2
        junction_rms_n+=1.
      junction_rms=(junction_rms/junction_rms_n)**0.5
    return junction_rms

  def get_junction_rms_n(self):
    return max(0,len(self.connected_group_segments)-1)

  def get_score(self,segment_weight=1.,connection_weight=-0.5,
     end_score=500.,complete_score=500):
    # score based on individual segments, minus penalty for connections, plus
    # bonus for getting ends and no gaps
    self.model_is_complete=False
    self.score=0.
    have_left_end=False
    have_right_end=False
    # score from individual segments
    for s in self.connected_group_segments:
      self.score+=s.get_score()*segment_weight
      if s.has_left_end(): have_left_end=True
      if s.has_right_end(): have_right_end=True
    if have_left_end:  # to ensure that we cover the ends
      self.score+=end_score
    if have_right_end:
      self.score+=end_score

    if have_left_end and have_right_end:
      expected_length=self.get_right_connection()-self.get_left_connection()+1
      actual_length=self.get_length()
      if expected_length==actual_length:
        self.score+=complete_score
        self.model_is_complete=True

    # score based on connections
    if len(self.connected_group_segments)>1:
      for s1,s2 in zip(
        self.connected_group_segments[:-1],self.connected_group_segments[1:]):
        self.score+=connection_weight*self.get_connection_rms(s1,s2)

    return self.score

  def get_length(self):  # does not include delta_resno (just the length it
    # would replace
    n=0
    for s in self.connected_group_segments:
      n+=s.get_length()
    # subtract off 1 for all but one (connection overlap)
    n=n-len(self.connected_group_segments)+1
    return n

  def get_connection_rms(self,s1,s2,
      resno=None,allow_unit_length=False):
    # get rms of atoms in connecting residue of segment 1 and segment 2

    # get main chain of last residue of s1 and first residue of s2 and compare
    if resno is None:
      s1_connection_resno=s1.end_resno
    else:
      s1_connection_resno=resno
    start_offset_s1=s1.start_resno-s1.segment.info['target_start_resno']
    last_resno_s1=get_first_resno(s1.segment.hierarchy)+\
       start_offset_s1+s1_connection_resno-s1.start_resno
    last_resno_s1+=s1.delta_length
    # do we already have it:
    dd=self.connection_rms_dict.get(s1.segment)
    if dd:
      dd=dd.get(s2.segment)
      if dd:
        value=dd.get(last_resno_s1)
        if value is not None: return value

    if resno is None:
      s2_connection_resno=s2.start_resno
    else:
      s2_connection_resno=resno
    start_offset_s2=s2_connection_resno-s2.segment.info['target_start_resno']
    first_resno_s2=get_first_resno(s2.segment.hierarchy)+start_offset_s2

    # make sure we create something long enough (not length 1):
    first_resno_s1=get_first_resno(s1.segment.hierarchy)
    start_offset_s2=s2.start_resno-s2.segment.info['target_start_resno']
    last_resno_s2=get_first_resno(s2.segment.hierarchy)+\
        start_offset_s2+s2.end_resno-s2.start_resno

    if not allow_unit_length and (
       last_resno_s1-first_resno_s1 < 1 or \
       last_resno_s2-first_resno_s2 < 1):  # would create a connection length 1
      rms=None
    else:
      rms=self.get_rms_pair_of_residues_main_chain(
        h1=s1.segment.hierarchy,resno1=last_resno_s1,
        h2=s2.segment.hierarchy,resno2=first_resno_s2)

    # save it so we don't calculate twice...
    if not s1.segment in self.connection_rms_dict:
      self.connection_rms_dict[s1.segment]={}
    if not s2.segment in self.connection_rms_dict[s1.segment]:
      self.connection_rms_dict[s1.segment][s2.segment]={}
    self.connection_rms_dict[s1.segment][s2.segment][last_resno_s1]=rms

    return rms

  def get_rms_pair_of_residues_main_chain(self,h1=None,resno1=None,
     h2=None,resno2=None,value_if_missing=10.):
    xyz1=self.get_one_residue_main_chain(hierarchy=h1,resno=resno1)
    xyz2=self.get_one_residue_main_chain(hierarchy=h2,resno=resno2)
    if xyz1.size()!=xyz2.size(): return value_if_missing
    diffs=xyz1-xyz2
    return diffs.rms_length()

  def get_one_residue_main_chain(self,hierarchy=None,resno=None):
    atom_selection="resseq %s and (name ca or name c or name o or name n)"  %(
     resseq_encode(resno))

    sele=apply_atom_selection(atom_selection,hierarchy=hierarchy)
    return sele.atoms().extract_xyz()


  def show_summary(self,out=None,return_text=True):
    rc=self.get_right_connection()
    lc=self.get_left_connection()
    if rc is None or lc is None:
      rc=0
      lc=0
    text="CG: Segments:%d Score: %7.2f Left: %d Right: %d Length %d" %(
     len(self.connected_group_segments),self.score,lc,rc,rc-lc+1)
    if return_text: return text
    if out is None: out=sys.stdout
    print(text, file=out)
    return ""

  def show_comprehensive_summary(self,out=None,return_text=True):
    text=self.show_summary(out=out,return_text=return_text)
    for cs,cs1 in zip(
        self.connected_group_segments,self.connected_group_segments[1:]+[None]):
      text+="\n%s" %(cs.show_summary())
      text+=" RMS:%5.2f A " %(cs.segment.info['rms'])
      if cs.segment.info['log_frequency']:
         text+=" ln(freq):%5.3f " %(cs.segment.info['log_frequency'])


      if cs1:
        text+=" Junction:%5.2f A " %(self.get_connection_rms(cs,cs1))
    if return_text: return text
    if out is None: out=sys.stdout
    print(text, file=out)
    return ""

class replacement_segment_summary:
  def __init__(self,model=None):
    self.model=model
    self.id=model.info['chain_number']
    self.chain_id=get_chain_id(self.model.hierarchy)
    self.replacement_model=None
    self.model_is_complete=None
    self.connected_groups=None
    self.model_has_insertions_deletions=None

  def show_summary(self,out=sys.stdout):
    print("\nID: %d ChainID: '%s'  RMSD: %5.2f A  (n=%d) " %(
       self.get_segment_id(),
       self.get_chain_id(),self.get_rms(),self.get_rms_n(),
        )+"Junction RMSD: %5.2f A (n=%d)" %(
      self.get_junction_rms(),self.get_junction_rms_n(),
       ), file=out)
    print("Complete: %s  Insertions/deletions: %s" %(
       self.is_complete(),self.has_insertions_deletions()), file=out)
    print("Input model start: %d  end: %d  length: %d " %(
         self.input_first_resno(),
         self.input_last_resno(),
         self.input_number_of_residues(),)+\
       "\nReplacement start: %d  end: %d  length: %d" %(
         self.output_first_resno(),
         self.output_last_resno(),
         self.output_number_of_residues(),
         ), file=out)

  def get_segment_id(self):
    return self.id

  def get_chain_id(self):
    return self.chain_id
  def add_replacement_model(self,replacement_model):
    self.replacement_model=replacement_model

  def add_connected_groups(self,connected_groups):
    self.connected_groups=connected_groups

  def get_rms(self):
    if self.connected_groups:
      sum_rms=0.
      sum_rms_n=0.
      for cg in self.connected_groups:
        if cg.get_rms():
          sum_rms+=cg.get_rms()**2*cg.get_rms_n()
          sum_rms_n+=cg.get_rms_n()
      if sum_rms_n:
        return (sum_rms/sum_rms_n)**0.5

    return 0.

  def get_rms_n(self):
    sum_rms_n=0.
    if self.connected_groups:
      for cg in self.connected_groups:
        if cg.get_rms_n():
          sum_rms_n+=cg.get_rms_n()
    return sum_rms_n

  def get_junction_rms(self):
    if self.connected_groups:
      sum_junction_rms=0.
      sum_junction_rms_n=0.
      for cg in self.connected_groups:
        if cg.get_junction_rms():
          sum_junction_rms+=cg.get_junction_rms()**2*cg.get_junction_rms_n()
          sum_junction_rms_n+=cg.get_junction_rms_n()
      if sum_junction_rms_n:
        return (sum_junction_rms/sum_junction_rms_n)**0.5

    return 0.

  def get_junction_rms_n(self):
    sum_junction_rms_n=0.
    if self.connected_groups:
      for cg in self.connected_groups:
        if cg.get_junction_rms_n():
          sum_junction_rms_n+=cg.get_junction_rms_n()
    return sum_junction_rms_n

  def set_is_complete(self,is_complete):
    self.model_is_complete=is_complete

  def is_complete(self):
    return self.model_is_complete

  def set_has_insertions_deletions(self,has_insertions_deletions):
    self.model_has_insertions_deletions=has_insertions_deletions

  def has_insertions_deletions(self):
    return self.model_has_insertions_deletions

  def input_number_of_residues(self):
    return self.model.hierarchy.overall_counts().n_residues

  def output_number_of_residues(self):
    if self.replacement_model:
      return self.replacement_model.hierarchy.overall_counts().n_residues
    else:
      return 0

  def input_first_resno(self):
    return get_first_resno(self.model.hierarchy)

  def input_last_resno(self):
    return get_last_resno(self.model.hierarchy)

  def output_first_resno(self):
    if self.replacement_model:
      return get_first_resno(self.replacement_model.hierarchy)
    else:
      return 0

  def output_last_resno(self):
    if self.replacement_model:
      return get_last_resno(self.replacement_model.hierarchy)
    else:
      return 0

class replace_with_segments_from_pdb:

  def __init__(self,pdb_hierarchy=None,args=None,
       helices_are_alpha=None,resid_offset=None,
       out=sys.stdout):

    self.replacement_model=model_info() # empty object
    self.model_output_number_of_residues_by_segment={}
    self.model_segment_ids=[]

    # get parameters
    params=self.get_params(args,out=out)

    if helices_are_alpha:
      params.find_ss_structure.helices_are_alpha=True

    if params.extract_segments_from_pdb.extract in [None,'None']:
      params.extract_segments_from_pdb.extract=None
      # get libraries
      helix_lib,strand_lib,other_lib=self.get_libraries(
         params,out=out)
    else:
      helix_lib,strand_lib,other_lib=None,None,None,None

    # find segments of secondary structure
    models=self.find_ss_from_pdb(params,pdb_hierarchy=pdb_hierarchy,out=out)

    # see if we can replace any secondary structure
    all_replacement_models=self.replace_secondary_structure(params,models=models,
        helix_lib=helix_lib,strand_lib=strand_lib,other_lib=other_lib,
         out=out)
    replacement_model=merge_hierarchies_from_models(
       models=all_replacement_models, resid_offset=resid_offset)

    # write out results
    if params.output_files.pdb_out:
      print("\nWriting output PDB file to %s" %(
        params.output_files.pdb_out), file=out)
      f=open(params.output_files.pdb_out,'w')
      print(replacement_model.hierarchy.as_pdb_string(), file=f)
      f.close()
    self.replacement_model=replacement_model

  def replacement_hierarchy(self):
    return self.replacement_model.hierarchy

  def output_number_of_residues_by_segment(self):
    return self.model_output_number_of_residues_by_segment

  def get_params(self,args,out=sys.stdout):
    command_line = iotbx.phil.process_command_line_with_files(
      reflection_file_def="input_files.map_coeffs_file",
      map_file_def="input_files.map_file",
      pdb_file_def="input_files.pdb_in",
      args=args,
      master_phil=master_phil)
    params = command_line.work.extract()
    print("\nRegularize from pdb: Use PDB segment to replace parts of a model", file=out)
    master_phil.format(python_object=params).show(out=out)
    return params


  def find_ss_from_pdb(self,params,pdb_hierarchy=None,out=sys.stdout):
    # read in model and split into segments
    models=get_and_split_model(pdb_in=params.input_files.pdb_in,
      pdb_hierarchy=pdb_hierarchy,
      get_start_end_length=True,out=out)

    # Note what we have
    self.model_replacement_segment_summaries=[]
    for model in models:
      rss=replacement_segment_summary(model=model)
      self.model_replacement_segment_summaries.append(rss)

    # identify secondary structure
    fss=find_secondary_structure(params=params,models=models,
     ss_by_chain=False,  # required (actually not as ignored for model input)
     helices_are_alpha=params.find_ss_structure.helices_are_alpha,
       out=out)
    return fss.models

  def get_libraries(self,params,out=sys.stdout):
    # read in libraries
    print("Reading libraries...", file=out)
    if params.control.verbose:
      lout=out
    else:
      lout=null_out()
    helix_lib=segment_library(segment_class=helix,params=params.alpha,out=lout)
    strand_lib=segment_library(segment_class=strand,params=params.beta,out=lout)
    other_lib=segment_library(segment_class=other,params=params.other,out=lout)

    print("Libraries read with %d helices and %d strands and %d other\n" %(
      len(helix_lib.segments),len(strand_lib.segments),len(other_lib.segments)), file=out)
    return helix_lib,strand_lib,other_lib

  def get_ss_structure(self,params,model=None,out=sys.stdout):
    if params.control.verbose:
      print("\nLooking for secondary structure in model %d with %d residues" %(
        model.info['chain_number'],model.hierarchy.overall_counts().n_residues)+\
        " starting at %d" %(get_first_resno(model.hierarchy)), file=out)

    if params.alpha.find_alpha:
      find_alpha=find_helix(params=params.alpha,model=model,
       extract_segments_from_pdb=params.extract_segments_from_pdb.extract,
         verbose=params.control.verbose,out=out)
      if params.control.verbose:
        find_alpha.show_summary(out=out)
    else:
      find_alpha=None

    if params.beta.find_beta:
      find_beta=find_beta_strand(params=params.beta,model=model,
         extract_segments_from_pdb=params.extract_segments_from_pdb.extract,
         verbose=params.control.verbose,out=out)
      if params.control.verbose:
        find_beta.show_summary(out=out)
    else:
      find_beta=None

    if params.other.find_other:
      find_other=find_other_structure(params=params.other,model=model,
         find_alpha=find_alpha,find_beta=find_beta,
         extract_segments_from_pdb=\
            params.extract_segments_from_pdb.extract,
         verbose=params.control.verbose,out=out)
      if params.control.verbose:
        find_other.show_summary(out=out)
    else:
      find_other=None

    return find_alpha,find_beta,find_other

  def sort_connected_groups(self,connected_groups,sort_by_start=False,
       small_number=-0.0001):
    sort_list=[]
    for cg in connected_groups:
      if sort_by_start:
        sort_list.append(
          [cg.get_left_connection()+cg.get_score()*small_number,cg])
      else:
        sort_list.append([cg.get_score(),cg])
    sort_list.sort(key = lambda x: x[0])
    if not sort_by_start: # high to low in score, low to high in sort_by_start
      sort_list.reverse()
    connected_groups=[]
    for sc,cg in sort_list:
      connected_groups.append(cg)
    return connected_groups

  def any_cg_has_insertions_deletions(self,connected_groups):
    for cg in connected_groups:
      if cg.has_insertions_deletions():  return True
    return False

  def remove_element(self,list_elements=[],element=None):
    if not element in list_elements:
      return list_elements
    else:
      new_list_elements=[]
      for x in list_elements:
         if x!=element: new_list_elements.append(x)
      return new_list_elements

  def combine_groups_by_residue(self,
      connected_groups,maximum_junction_rmsd=None,
      verbose=None,out=sys.stdout):
    connected_groups=self.sort_connected_groups(connected_groups)
    # sorted on score
    found=True
    while found:
      used_groups=[]
      found=False
      for cg in connected_groups:
        if found: break
        if cg in used_groups: continue
        used_groups.append(cg)
        for cg1 in connected_groups[1:]:
          if found: break
          if cg1 in used_groups: continue
          shared_residue_range=cg.shared_residue_range(cg1)
          if shared_residue_range:
            left_cg,right_cg=cg.get_left_right_cg(cg,cg1)
            original_right_cg_score=right_cg.get_score()

            found=left_cg.join_by_residue(other=right_cg,
              shared_residue_range=shared_residue_range,
              maximum_junction_rmsd=maximum_junction_rmsd)
            if found:
              combined_score=left_cg.get_score()
              used_groups=self.remove_element(
                 list_elements=used_groups,element=left_cg)
              new_groups=[]
              for poss_cg in connected_groups:
                if poss_cg != right_cg or original_right_cg_score>combined_score:
                  new_groups.append(poss_cg)
              connected_groups=new_groups
              connected_groups=self.sort_connected_groups(connected_groups)
              break

    connected_groups=self.sort_connected_groups(
       connected_groups,sort_by_start=True)
    if verbose:
      print("\nNumber of groups after combining by shared residue: %d"%(
         len(connected_groups)), file=out)
      print("\nConnected groups after combining by shared residue:", file=out)
      for cg in connected_groups:
        print(cg.show_summary(out=out))

    # remove overlapping now
    connected_groups=self.sort_and_remove_overlapping_groups(connected_groups,
      trim=True,contained_only=False)

    return connected_groups

  def combine_groups(self,connected_groups,verbose=None,out=sys.stdout):
    connected_groups=self.sort_connected_groups(connected_groups)
    # sorted on score
    found=True
    while found:
      used_groups=[]
      found=False
      new_cg=None
      for cg in connected_groups:
        if cg in used_groups: continue
        used_groups.append(cg)
        if found: break
        for cg1 in connected_groups[1:]:
          if cg1 in used_groups: continue
          if found: break
          if cg.shared_segment(cg1):
            found=cg.join(cg1)
            if found:
              connected_groups=self.sort_and_remove_overlapping_groups(
                connected_groups)
              connected_groups=self.sort_connected_groups(connected_groups)
    connected_groups=self.sort_connected_groups(
       connected_groups,sort_by_start=True)

    if verbose:
      print("\nNumber of groups after combining by shared segment: %d"%(
         len(connected_groups)), file=out)
      print("\nConnected groups after combining by shared segment:", file=out)
      for cg in connected_groups:
        print(cg.show_summary(out=out))

    connected_groups=self.sort_and_remove_overlapping_groups(
       connected_groups,out=out)

    return connected_groups

  def sort_and_remove_overlapping_groups(self,
       connected_groups,contained_only=True,
       trim=False,verbose=None,out=sys.stdout):

    connected_groups=self.sort_connected_groups(
        connected_groups) # Sorted on score

    if trim:
      connected_groups=self.trim_overlapping_groups(connected_groups)
    else:
      connected_groups=self.remove_overlapping_groups(
         connected_groups,contained_only=contained_only)

    # now sort on start_position
    connected_groups=self.sort_connected_groups(
       connected_groups,sort_by_start=True)

    if verbose:
      print("\nNumber of groups after removing overlapping: %d"%(
         len(connected_groups)), file=out)
      print("\nConnected groups after removing overlapping:", file=out)
      for cg in connected_groups:
        print(cg.show_summary(out=out))

    return connected_groups

  def trim_overlapping_groups(self,connected_groups):

    cg_dict={}
    for cg in connected_groups:
      cg_dict[cg]={}
      for i in range(cg.get_left_connection(),cg.get_right_connection()+1):
        cg_dict[cg][i]=None  # present

    new_groups=[]
    # mark all residues covered by a better segment and keep remaining good ones
    for i in range(len(connected_groups)):
      cg=connected_groups[i]
      keys=list(cg_dict[cg].keys())
      keys.sort()
      ok=cg.trim_residues(trim_dict=cg_dict[cg])
      if ok:
        new_groups.append(cg)
        for j in range(i+1,len(connected_groups)):
          cg_j=connected_groups[j]
          for k in range(cg.get_left_connection(),cg.get_right_connection()+1):
            if k in cg_dict[cg_j].keys():
              cg_dict[cg_j][k]=True


    return new_groups

  def remove_overlapping_groups(self,connected_groups,contained_only=False):
     removed_list=[None]*len(connected_groups)
     for i in range(len(connected_groups)):
       for j in range(i+1,len(connected_groups)):
         if removed_list[j]:
            continue # already removed
         if contained_only:
           if connected_groups[i].contains(connected_groups[j]):
             removed_list[j]=True
         elif connected_groups[i].overlaps_with(connected_groups[j]):
            removed_list[j]=True
     new_list=[]
     for cg,rm in zip(connected_groups,removed_list):
       if not rm:
         new_list.append(cg)
     return new_list


  def get_left_and_right_connections_for_segment(self,rs,minimum_length=2,
      maximum_overlap=1):
    target_start=rs.info['target_start_resno']
    target_end=rs.info['target_end_resno']
    length=rs.info['length']
    minimum_overlap=rs.info['minimum_overlap']
    target_length=rs.info['target_length']
    connection_pairs=[]
    if target_length!=length:  # allows variable lengths
      connection_pairs.append([target_start,target_end])
    else:
      start_of_left_connections=target_start+minimum_overlap
      end_of_left_connections=target_start+maximum_overlap
      start_of_right_connections=target_end-maximum_overlap
      end_of_right_connections=target_end-minimum_overlap
      if end_of_left_connections>=end_of_right_connections:
        end_of_left_connections=end_of_right_connections-1
      if start_of_right_connections<=start_of_left_connections:
        start_of_right_connections=start_of_left_connections+1

      for i in range(start_of_left_connections,end_of_left_connections+1):
        for j in range(start_of_right_connections,end_of_right_connections+1):
          if i>=j: continue
          connection_pairs.append([i,j])

    return connection_pairs

  def get_model_from_connected_groups(self,connected_groups,
       model_to_match=None,renumber=None,verbose=None,out=sys.stdout):
    models=[]
    for cg,is_left_end,is_right_end in zip(
       connected_groups,
       [True]+[False]*(len(connected_groups)-1),
       [False]*(len(connected_groups)-1)+[True]):
      models.append(
        cg.as_model(is_left_end=is_left_end,is_right_end=is_right_end,
        model_to_match=model_to_match,verbose=verbose,out=out))
    if not models:
       return None
    if renumber:
      new_model=merge_hierarchies_from_models(models=models,resid_offset=100,
        first_residue_number=get_first_resno(models[0].hierarchy))
    else:
      new_model=merge_hierarchies_from_models(models=models,renumber=False)
    return new_model

  def sort_connection_list(self,connection_list):
    connection_list.sort()
    connection_list.reverse()
    connections=[]
    best_score=0.0
    first=True
    for [score,connection] in connection_list:
      connections.append(connection)
      if first:
        best_score=score
        first=False
    return connections,best_score

  def get_all_connections(self,cg=None,used_connected_groups=None,
        possible_left_connection_dict=None,
        possible_right_connection_dict=None,
        maximum_junction_rmsd=None,look_ahead_score_only=False,
        make_right_connections=True,make_left_connections=True):
    found=False
    best_score=0.
    # make all right and left connections possible to this cg
    right_connection_list=[]
    right_connections=[]
    for possible_right_connection in possible_right_connection_dict.get(
           cg.get_right_connection(),[]):
      if not make_right_connections:continue
      if possible_right_connection in used_connected_groups:
        continue
      # make sure this connection does not duplicate one we have
      if cg.get_right_segment_hierarchy()==\
              possible_right_connection.get_left_segment_hierarchy():
        continue
      rms=cg.get_connection_rms(cg.connected_group_segments[-1],
             possible_right_connection.connected_group_segments[0])
      if rms is None or rms > maximum_junction_rmsd:
        continue
      score=possible_right_connection.get_score()
      if not look_ahead_score_only:
        best_score=self.get_all_connections(
          cg=possible_right_connection,
          used_connected_groups=used_connected_groups,
          possible_left_connection_dict=possible_left_connection_dict,
          possible_right_connection_dict=possible_right_connection_dict,
          maximum_junction_rmsd=maximum_junction_rmsd,make_left_connections=False,
          look_ahead_score_only=True)
        score+=best_score # adds on score for best right connection
      right_connection_list.append([score,possible_right_connection])
    right_connections,best_score=self.sort_connection_list(right_connection_list)
    if look_ahead_score_only:
      return best_score
    # Add on the best one including look-ahead score
    if right_connections:
      cg.add_connected_group_to_right(right_connections[0])
      used_connected_groups+=right_connections
      found=True

    # and now for left direction
    left_connection_list=[]
    left_connections=[]
    for possible_left_connection in possible_left_connection_dict.get(
           cg.get_left_connection(),[]):
      if not make_left_connections:continue
      if possible_left_connection in used_connected_groups:
        continue
      # make sure this connection does not duplicate one we have
      if cg.get_left_segment_hierarchy()==\
              possible_left_connection.get_right_segment_hierarchy():
        continue
      rms=cg.get_connection_rms(
             possible_left_connection.connected_group_segments[-1],
             cg.connected_group_segments[0],)
      if rms is None or rms > maximum_junction_rmsd:
        continue
      score=possible_left_connection.get_score()
      if not look_ahead_score_only:
        best_score=get_all_connections(
          cg=possible_left_connection,
          used_connected_groups=used_connected_groups,
          possible_left_connection_dict=possible_left_connection_dict,
          possible_right_connection_dict=possible_right_connection_dict,
          maximum_junction_rmsd=maximum_junction_rmsd,
          make_right_connections=False,
          look_ahead_score_only=True)
        score+=best_score # adds on score for best right connection
      left_connection_list.append(
        [score,possible_left_connection])
    left_connections,best_score=self.sort_connection_list(left_connection_list)
    if look_ahead_score_only:
      return best_score
    # Add on the best one
    if left_connections:
      new_cg=left_connections[0]
      new_cg.add_connected_group_to_right(cg)
      used_connected_groups+=left_connections
      cg=new_cg
      found=True
    return cg,used_connected_groups,found

  def link_groups(self,params=None,model=None,connected_groups=None,
      maximum_junction_rmsd=None,
      other_lib=None,out=sys.stdout):

    found=False
    if len(connected_groups) < 2: return connected_groups,found

    if params.control.verbose:
      print("\nLinking groups that are not yet connected", file=out)
    connected_groups=self.sort_connected_groups(connected_groups,
      sort_by_start=True) # get them in order
    new_connected_groups=[]
    good_enough_rmsd=params.regularization.good_enough_ratio*\
                 params.other.maximum_rmsd
    for cg1,cg2 in zip(connected_groups[:-1],connected_groups[1:]):
      if cg1.get_right_connection()+1==cg2.get_left_connection():

        # Create a segment with the junction in it. Then see if we can replace
        #   it with something that works and matches all the main-chain atoms
        #   at crossover points

        # see if we can find a segment from our other_lib that covers some of
        # the residues in this junction.  Try 1234 residues in from either end
        # create a hierarchy with cg1+cg2

        # if cg1 has an insert, start it with lower residue numbers to compensate
        cg1_model=cg1.as_model(offset_to_renumber_right_end_at_target=True)
        cg2_model=cg2.as_model()

        combined_model=merge_hierarchies_from_models(
          models=[cg1_model,cg2_model],
          resid_offset=1,
          first_residue_number=get_first_resno(cg1_model.hierarchy),
          renumber=True)

        sites=combined_model.hierarchy.extract_xray_structure().sites_cart()

        # now select a few possibilities with at least 1 res overlap on each end
        #  and no more than other.standard_length residues
        first_resno=get_first_resno(combined_model.hierarchy)
        last_resno=get_last_resno(combined_model.hierarchy)
        before_junction=get_last_resno(cg1_model.hierarchy)

        target_overlap=2
        i_start=max(first_resno,before_junction-target_overlap+1)
        i_end=min(last_resno-2*target_overlap+1,before_junction)
        if i_start>i_end: continue # nothing to do

        if other_lib.index_length is None:
          start_position=None
        else:
          start_position=max(0,1+(other_lib.index_length-2*target_overlap)//2)

        best_cg=None
        best_rms=None
        segment_list=[]
        for i in range(i_start,i_end+1):
          atom_selection=\
           "resid %s through %s and (name n or name c or name ca or name o)" %(
            resseq_encode(i),resseq_encode(i+2*target_overlap-1))
          segment_list.append(other(params=params.other,
            start_resno=i,
            hierarchy= apply_atom_selection(
            atom_selection,hierarchy=combined_model.hierarchy)))
        for segment in segment_list:
          # now find replacements for this segment
          rs,rms,log_frequency,lc,rc,lc_other,rc_other=\
             other_lib.get_replacement_segment(
              model_segment=segment,
              maximum_rmsd=params.other.maximum_rmsd,
              good_enough_rmsd=good_enough_rmsd,
              start_position=start_position,
              score_by_main_chain_rmsd=True,
              before_junction=before_junction,
              delta_residues=0,out=out)
          if rms is not None and rms < params.other.maximum_rmsd:
            if best_cg is None or best_rms> rms:
              # make a connected_segment object with residues from
              #lc_other to rc_other and link it
              contains_left_end_of_chain,contains_right_end_of_chain=\
                segment.contains_ends(
                  first_resno_of_chain=model.info['model_start_resno'],
                  last_resno_of_chain=model.info['model_end_resno'])
              replacement_segment=model_info(hierarchy=rs,
                info={'target_start_resno':segment.get_start_resno(),
                     'target_end_resno':segment.get_end_resno(),
                     'target_length':
                        segment.get_end_resno()-segment.get_start_resno()+1,
                     'length':rs.overall_counts().n_residues,
                     'minimum_overlap':params.other.minimum_overlap,
                     'rms':rms,
                     'base_score':segment.base_score,
                     'contains_right_end_of_chain':contains_right_end_of_chain,
                     'contains_left_end_of_chain':contains_left_end_of_chain,
                     'log_frequency':log_frequency},
                )
              best_cg=connected_group(segment=replacement_segment,
                 start_resno=lc,end_resno=rc)
              best_rms=rms
              if rms <= good_enough_rmsd:
                break # done
        if best_rms:
           best_cg.get_score()
           new_connected_groups.append(best_cg)
           found=True
    connected_groups+=new_connected_groups

    if params.control.verbose:
      print("\nGroups after adding linker groups:", file=out)
      for cg in connected_groups: print(cg.show_comprehensive_summary())
      for cg in connected_groups:
        print("\nGROUP:")
        print(cg.as_model(is_left_end=True).hierarchy.as_pdb_string())

    connected_groups=self.sort_and_remove_overlapping_groups(connected_groups)

    if params.control.verbose:
      print("\nGroups after removing overlapping:", file=out)
      for cg in connected_groups: print(cg.show_comprehensive_summary())
      for cg in connected_groups:
        print("\nGROUP:")
        print(cg.as_model(is_left_end=True).hierarchy.as_pdb_string())

    connected_groups=self.combine_groups_by_residue(connected_groups,
       maximum_junction_rmsd=maximum_junction_rmsd,
       verbose=params.control.verbose,out=out)

    return connected_groups,found

  def get_connected_groups_from_segments(self,params=None,
       replacement_segments=None,
       maximum_overlap=None,
       verbose=None,out=sys.stdout):
    # Set up connected groups
    connected_groups=[]
    possible_right_connection_dict={}
    possible_left_connection_dict={}
    for rs in replacement_segments:
      connection_pairs=self.get_left_and_right_connections_for_segment(rs,
        maximum_overlap=maximum_overlap)
      for lc,rc in connection_pairs:
          cg=connected_group(segment=rs,start_resno=lc,end_resno=rc)
          connected_groups.append(cg)
          if not lc in possible_right_connection_dict:
            possible_right_connection_dict[lc]=[]
          possible_right_connection_dict[lc].append(cg)
          if not rc in possible_left_connection_dict:
            possible_left_connection_dict[rc]=[]
          possible_left_connection_dict[rc].append(cg)
    # here possible_right_connection_dict[rc] is a list of connected_groups that
    #   can be linked to the right with their left end residue equal to rc

    connected_groups=self.sort_connected_groups(
       connected_groups,sort_by_start=True)

    if verbose:
      print("\nNumber of starting connected groups: %d"%(
          len(connected_groups)), file=out)
      print("\nConnected groups:", file=out)
      for cg in connected_groups:
        print(cg.show_summary(out=out))

    return connected_groups,possible_left_connection_dict,\
        possible_right_connection_dict

  def connect_groups(self,connected_groups=None,
     possible_left_connection_dict=None,
     possible_right_connection_dict=None,
     maximum_junction_rmsd=None,
     verbose=None,out=sys.stdout):

    if verbose:
      print("\nCreating new connected groups", file=out)

    used_connected_groups=[]
    final_connected_groups=[]
    for cg in connected_groups:  # work down the list
      if cg in used_connected_groups: continue
      used_connected_groups.append(cg)
      while 1:
        cg,used_connected_groups,found=self.get_all_connections(cg=cg,
          used_connected_groups=used_connected_groups,
          possible_left_connection_dict=possible_left_connection_dict,
          possible_right_connection_dict=possible_right_connection_dict,
          maximum_junction_rmsd=maximum_junction_rmsd)
        if not found: break
      final_connected_groups.append(cg)
    connected_groups=final_connected_groups
    connected_groups=self.sort_and_remove_overlapping_groups(connected_groups)
    connected_groups=self.sort_connected_groups(
       connected_groups,sort_by_start=True)

    if verbose:
      print("\nNumber of connected groups after connecting: %d"%(
        len(connected_groups)), file=out)
      print("\nConnected groups after connecting:", file=out)
      for cg in connected_groups:
        print(cg.show_summary(out=out))

    return connected_groups

  def get_connections(self,params,model=None,
     replacement_segments=None,maximum_overlap=None,
     maximum_junction_rmsd=None,other_lib=None,out=sys.stdout):
    # make longest possible contiguous segments from the connections available

    # A replacement_segment is a model_info object with a hierarchy that
    #  covers the whole segment and has info['target_start_resno'] etc that
    #  specifies what start residue
    #  we plan to use for this hierarchy (will not normally match the residue
    #  number in the hierarchy which is arbitrary)

    if params.control.verbose:
      print("\nReplacement segments:", file=out)
      for rs in replacement_segments:
        print(rs.show_summary(out=out))

    # Initial single connected groups and possible left and right connections
    connected_groups,possible_left_connection_dict,\
         possible_right_connection_dict=\
      self.get_connected_groups_from_segments(params=params,
        maximum_overlap=maximum_overlap,
        replacement_segments=replacement_segments,out=out)

    # link up the connected groups using left and right connections
    connected_groups=self.connect_groups(connected_groups=connected_groups,
        possible_left_connection_dict=possible_left_connection_dict,
        possible_right_connection_dict=possible_right_connection_dict,
        maximum_junction_rmsd=maximum_junction_rmsd,
        verbose=params.control.verbose,out=out)

    # combine groups by shared segments
    connected_groups=self.combine_groups(connected_groups,
      verbose=params.control.verbose,out=out)

    # combine groups by shared residue
    connected_groups=self.combine_groups_by_residue(connected_groups,
       maximum_junction_rmsd=maximum_junction_rmsd,
       verbose=params.control.verbose,out=out)

    n=0
    found=True
    while found:
      n+=1
      # combine groups by finding spanning segment
      connected_groups,found=self.link_groups(params=params,
        model=model,connected_groups=connected_groups,
        other_lib=other_lib,maximum_junction_rmsd=maximum_junction_rmsd,
        out=out)

      # final cleanup
      connected_groups=self.sort_and_remove_overlapping_groups(connected_groups,
        trim=True,contained_only=False)
      if n>1: break # makes infinite loop

    if params.control.verbose:
      print("\nFinal groups for this segment: ", file=out)
      for cg in connected_groups:
        print(cg.show_comprehensive_summary(), file=out)

    return connected_groups

  def assemble_segments(self,params,model=None,
        other_lib=None,
        replacement_segments=None,
        out=sys.stdout):
    # We have a set of segments with target_start_resno, target_end_resno,
    #    length (actual), target_length.  If length==target_length, any residue
    #    can be used as crossover, otherwise, only the end residues

    if params.control.verbose:
      print("\nAssembling segments for model. Start:"+\
         " %d length: %d Replacement segments: %d" %(
       get_first_resno(model.hierarchy),
       model.hierarchy.overall_counts().n_residues,
       len(replacement_segments)), file=out)
    connected_groups=self.get_connections(params,
      replacement_segments=replacement_segments,
      model=model,
      other_lib=other_lib,
      maximum_junction_rmsd=params.regularization.maximum_junction_rmsd,
      maximum_overlap=params.regularization.maximum_overlap,out=out)
    return connected_groups

  def replace_secondary_structure(self,params,models=None,
      helix_lib=None,strand_lib=None,other_lib=None,
       out=sys.stdout):
    all_replacement_models=[]
    completeness_of_all_replacement_models=[]
    insertions_deletions_of_all_replacement_models=[]
    for rss in self.model_replacement_segment_summaries:
      model=rss.model
      # get segments that might replace structure in this model (1 chain)
      replacement_segments=self.get_replacement_segments(
        params,model=model,helix_lib=helix_lib,strand_lib=strand_lib,
        other_lib=other_lib,out=out)

      if params.extract_segments_from_pdb.extract is not None:
        all_replacement_models+=replacement_segments
        insertions_deletions_of_all_replacement_models.append(True)
        completeness_of_all_replacement_models.append(True)
      else:
        print("\nReplacing segment %d (n=%d," %(
            model.info['chain_number'],
            model.hierarchy.overall_counts().n_residues)+\
            " %d - %d) ..." %(
           get_first_resno(model.hierarchy),get_last_resno(model.hierarchy)), file=out)

        connected_groups=self.assemble_segments(params,model=model,
          other_lib=other_lib,
          replacement_segments=replacement_segments,out=out)
        if connected_groups:
          is_complete=(
            len(connected_groups)==1 and connected_groups[0].is_complete())
          has_insertions_deletions=self.any_cg_has_insertions_deletions(
            connected_groups)

          replacement_model=self.get_model_from_connected_groups(
            connected_groups,
            renumber=has_insertions_deletions,
            model_to_match=model,verbose=params.control.verbose,out=out)

          rss.add_connected_groups(connected_groups)
          rss.add_replacement_model(replacement_model)
          rss.set_is_complete(is_complete)
          rss.set_has_insertions_deletions(
            self.any_cg_has_insertions_deletions(connected_groups))

          all_replacement_models.append(replacement_model)
          completeness_of_all_replacement_models.append(is_complete)
          insertions_deletions_of_all_replacement_models.append(
            has_insertions_deletions)
          id=model.info['chain_number']
          start_residue=get_first_resno(replacement_model.hierarchy)
          self.model_output_number_of_residues_by_segment[id]=\
            replacement_model.hierarchy.overall_counts().n_residues
          if params.control.verbose:
            print("Replacement model for segment %d with %d residues " %(
              model.info['chain_number'],
              replacement_model.hierarchy.overall_counts().n_residues) + \
             " from %d to %d:" %(
             get_first_resno(model.hierarchy),get_last_resno(model.hierarchy)), file=out)
        else:
          print("No replacement model found for this segment", file=out)
          all_replacement_models.append(None)
          completeness_of_all_replacement_models.append(None)
          insertions_deletions_of_all_replacement_models.append(None)

    return all_replacement_models

  def get_replacement_segments(self,params,model=None,
      helix_lib=None,strand_lib=None,other_lib=None,
      out=sys.stdout):

      replacement_segments=[]

      for replacement_type,segment_lib,\
          segment_params  in zip(
         ['alpha','beta','other'],
         [helix_lib,strand_lib,other_lib],
         [params.alpha,params.beta,params.other]):
        find_object=getattr(model,'find_%s' %(replacement_type))
        if not find_object or not segment_params.replace :continue
        for segment in find_object.get_segments():
          if params.extract_segments_from_pdb.extract is not None:
            if params.extract_segments_from_pdb.extract==replacement_type:
              replacement_segments.append(model_info(hierarchy=segment.hierarchy))
            else:
              pass # do nothing, just writing out one type of structure
          else: # usual
            delta_residues=segment.optimal_delta_length
            if delta_residues and params.control.verbose:
              print("\nReplacing segment from model length " + \
                "%d with library segment length %d"%(
                segment.length(),segment.optimal_delta_length+segment.length(),), file=out)
            elif params.control.verbose:
              print("\nReplacing segment from model length " +\
                "%d with library strand" %(segment.length()), file=out)
            rs,rms,log_frequency,lc,rc,lc_other,rc_other=\
               segment_lib.get_replacement_segment(
              model_segment=segment,
              maximum_rmsd=segment_params.maximum_rmsd,
              good_enough_rmsd=\
              params.regularization.good_enough_ratio*segment_params.maximum_rmsd,
              delta_residues=delta_residues,
              verbose=params.control.verbose,out=out)
            if rs: # Note length will be: segment.length()+delta_residues
              contains_left_end_of_chain,contains_right_end_of_chain=\
                segment.contains_ends(
                  first_resno_of_chain=model.info['model_start_resno'],
                  last_resno_of_chain=model.info['model_end_resno'])

              replacement_segments.append(model_info(hierarchy=rs,
                info={'target_start_resno':segment.get_start_resno(),
                     'target_end_resno':segment.get_end_resno(),
                     'target_length':
                        segment.get_end_resno()-segment.get_start_resno()+1,
                     'length':rs.overall_counts().n_residues,
                     'minimum_overlap':segment_params.minimum_overlap,
                     'rms':rms,
                     'base_score':segment.base_score,
                     'contains_right_end_of_chain':contains_right_end_of_chain,
                     'contains_left_end_of_chain':contains_left_end_of_chain,
                     'log_frequency':log_frequency},
                ))

      return replacement_segments


if __name__=="__main__":

  from mmtbx.secondary_structure.regularize_from_pdb import replace_with_segments_from_pdb
  r=replace_with_segments_from_pdb(args=sys.argv[1:],pdb_hierarchy=None,
    out=sys.stdout)

  print_output=False
  if print_output:
    print("Output model:\n%s" %(r.replacement_hierarchy().as_pdb_string()))


  print("\nSummary by segment:")
  for rss in r.model_replacement_segment_summaries:
    rss.show_summary()
