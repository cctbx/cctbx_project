from __future__ import division

# find_ss_from_ca.py
# a tool to find helices, strands, non-helices/strands in segments of
#  structure

from iotbx.pdb import resseq_encode
import iotbx.phil
import os,sys
from libtbx.utils import Sorry,null_out
from scitbx.matrix import col
from scitbx.math import superpose, matrix
from scitbx.array_family import flex
from copy import deepcopy

master_phil = iotbx.phil.parse("""

  input_files {
    pdb_in = None
      .type = path
      .help = Input PDB file
      .short_caption = Input PDB file
  }

  find_ss_structure {

     find_alpha = True
       .type = bool
       .help = Find alpha helices
       .short_caption = Find alpha helices

     find_beta = True
       .type = bool
       .help = Find beta structure
       .short_caption = Find beta structure

     find_other = False
       .type = bool
       .help = Find other structure
       .short_caption = Find other structure

     exclude_alpha_in_beta  = True
       .type = bool
       .help = Exclude regions already identified as alpha
       .short_caption = Exclude alpha regions from beta

     exclude_alpha_beta_in_other  = True
       .type = bool
       .help = Exclude regions already identified as alpha or beta
       .short_caption = Exclude alpha and beta regions from other

     make_unique = True
       .type = bool
       .help = Assign each residue to a unique type of structure
       .short_caption = Assign residues to unique structure

  }

  alpha {
    include scope mmtbx.secondary_structure.secondary_structure_params.alpha_params
  }

  beta {
    include scope  mmtbx.secondary_structure.secondary_structure_params.beta_params
  }

  other {
    include scope  mmtbx.secondary_structure.secondary_structure_params.other_params
  }


  extract_segments_from_pdb {
    extract = *None alpha beta other
       .type = choice
       .help = Extract all segments (alpha/beta) from a PDB file. \
               Used to create libraries of segments
       .short_caption = Extract segments from PDB
  }

  control {
      verbose = False
        .type = bool
        .help = Verbose output
        .short_caption = Verbose output
  }
""", process_includes=True)
master_params = master_phil

def apply_atom_selection(atom_selection,hierarchy=None):
  asc=hierarchy.atom_selection_cache()
  sel = asc.selection(string = atom_selection)
  return hierarchy.deep_copy().select(sel)  # deep copy is required

def get_pdb_hierarchy(text=None):
  return iotbx.pdb.input(
     source_info=None,lines=flex.split_lines(text)).construct_hierarchy()

def split_model(model=None,hierarchy=None,verbose=False,info=None,
     out=sys.stdout):

  model_list=[]
  if hierarchy:
    if not info: info={}
  else: # a model object with hierarchy and info
    hierarchy=model.hierarchy
    info=model.info
  for m in hierarchy.models():
    for chain in m.chains():
      new_hierarchy=iotbx.pdb.pdb_input(
         source_info="Model", lines="").construct_hierarchy()
      mm=iotbx.pdb.hierarchy.model()
      cc=iotbx.pdb.hierarchy.chain()
      cc.id=chain.id  # copy chain ID
      new_hierarchy.append_model(mm)
      mm.append_chain(cc)
      last_resseq=None
      for r in chain.residue_groups():
        if not (last_resseq is None or r.resseq_as_int()==last_resseq+1):
          # save and make new model
          new_model_info=model_info(hierarchy=new_hierarchy,info=deepcopy(info))
          model_list.append(new_model_info)
          new_model_info.info['chain_number']=len(model_list)

          # and make a new one
          new_hierarchy=iotbx.pdb.pdb_input(
             source_info="Model", lines="").construct_hierarchy()
          mm=iotbx.pdb.hierarchy.model()
          cc=iotbx.pdb.hierarchy.chain()
          cc.id=chain.id  # copy chain ID
          new_hierarchy.append_model(mm)
          mm.append_chain(cc)
          last_resseq=None
        # add on a residue...
        cc.append_residue_group(r.detached_copy())
        last_resseq=r.resseq_as_int()
      new_model_info=model_info(hierarchy=new_hierarchy,info=deepcopy(info))
      model_list.append(new_model_info)
      new_model_info.info['chain_number']=len(model_list)
  if verbose:
    print >>out,"Models after splitting:"
    for m in model_list:
      print >>out,"Chain: %d  Residues: %d" %(
        m.info.get('chain_number'),
        m.hierarchy.overall_counts().n_residues)
  return model_list

def merge_hierarchies_from_models(models=None,resid_offset=None,
    renumber=None,first_residue_number=None,
    sequences=None,chain_id=None,trim_side_chains=None):
  # assumes one chain from each model
  # if resid_offset, space by to next even n of this number of residues
  # otherwise if renumber, start at first_residue_number and sequence
  #  consecutively
  # If sequence or chain_id are supplied, use them
  # Trim off side chains (and CB for GLY) if trim_side_chains

  new_hierarchy=iotbx.pdb.pdb_input(
         source_info="Model", lines="").construct_hierarchy()
  mm=iotbx.pdb.hierarchy.model()
  new_hierarchy.append_model(mm)
  if renumber:  # just get chain once in advance
    cc=iotbx.pdb.hierarchy.chain()
    if chain_id:
      cc.id=chain_id
    mm.append_chain(cc)

  # set resid if necessary
  if resid_offset:
    if first_residue_number:
      resid=first_residue_number
    else:
      resid=1
  elif renumber:
    if first_residue_number is None:
      raise Sorry(
        "Need first_residue_number for merge_hierarchies_from_models if "+
         "renumber=True")
    resid=first_residue_number

  info={}
  from iotbx.pdb import resseq_encode

  if not sequences: sequences=len(models)*[None]
  for model,sequence in zip(models,sequences):
    if not model: continue
    if not info: info=model.info
    i=0
    for m in model.hierarchy.models():
      for chain in m.chains():
        if not renumber:
          cc=iotbx.pdb.hierarchy.chain()
          mm.append_chain(cc)
        for r in chain.residue_groups():
          rr=r.detached_copy()
          if resid_offset or renumber:
            rr.resseq=resseq_encode(resid)
            resid+=1
          if sequence:
            resname=sequence[i]
            i+=1
            for atom_group in rr.atom_groups():
              atom_group.resname=resname
          cc.append_residue_group(rr)
        if resid_offset and resid_offset>1:
          nn=resid_offset*((resid+resid_offset-1)//resid_offset)
          if nn<2: nn+=resid_offset

  if trim_side_chains:
    atom_selection=\
      "name ca or name c or name o or name n or (name cb and not resname gly)"
    new_hierarchy=apply_atom_selection(atom_selection,hierarchy=new_hierarchy)

  new_model=model_info(hierarchy=new_hierarchy,info=info)

  return new_model

def get_average_direction(diffs=None, i=None,j=None):
    if not diffs: return None
    if i is None and j is None:
      i=0
      j=len(diffs)-1
    average_direction=col(diffs[i])
    nn=1.
    for j in xrange(i+1,j):
      nn+=1.
      average_direction+=col(diffs[j])
    average_direction/=nn
    return average_direction.normalize()

def get_chain_id(hierarchy):
  if not hierarchy:
    return None
  for model in hierarchy.models():
    for chain in model.chains():
      return chain.id
  return None # nothing there

def get_sequence(hierarchy):
  if not hierarchy:
    return None
  sequence=[]
  for model in hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        for atom_group in rg.atom_groups():
          sequence.append(atom_group.resname)
          break
  return sequence

def get_atom_list(hierarchy):
  atom_list=[]
  if not hierarchy:
    return atom_list
  for model in hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        for atom_group in rg.atom_groups():
          for atom in atom_group.atoms():
            atom_list.append(atom.name)
  return atom_list

def get_first_resno(hierarchy):
  if not hierarchy:
    return None
  for model in hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        return rg.resseq_as_int()

def get_last_resno(hierarchy):
  if not hierarchy:
    return None
  last_resno=None
  for model in hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        last_resno=rg.resseq_as_int()
  return last_resno

class model_info: # mostly just a holder
  def __init__(self,hierarchy=None,id=0,info={},
      find_alpha=None,find_beta=None,find_other=None):
    self.hierarchy=hierarchy
    self.info=info
    self.id=id
    self.find_alpha=find_alpha
    self.find_beta=find_beta
    self.find_other=find_other

  def show_summary(self,out=sys.stdout):
    keys=self.info.keys()
    keys.sort()
    print >>out,"\nModel %d" %(self.id)
    for key in keys:
      print >>out,"%s: %s" %(key,str(self.info[key]))

  def first_residue(self):
    return get_first_resno(self.hierarchy)

  def last_residue(self):
    return get_last_resno(self.hierarchy)
  
  def length(self):
    return self.last_residue()-self.first_residue()+1

class segment:  # object for holding a helix or a strand or other

  def setup(self,sites=None,start_resno=None,hierarchy=None,
     segment_type='None',
     span=None,target_rise=None,
     rise_tolerance=None,dot_min=None,dot_min_single=None,
     minimum_length=None,
     buffer_residues=None,
     standard_length=None,
     is_n_terminus=None,
     is_c_terminus=None,
     frequency=None,
     base_score=None,
     optimal_delta_length=None,
     verbose=None,
     out=sys.stdout):

    self.out=out
    self.verbose=verbose
    self.orientation_points=None
    self.orientation_points_start=None
    self.orientation_points_end=None
    self.is_n_terminus=is_n_terminus
    self.is_c_terminus=is_c_terminus
    self.buffer_residues=buffer_residues
    self.standard_length=standard_length
    self.minimum_length=minimum_length
    self.segment_type=segment_type
    self.frequency=frequency
    self.base_score=base_score
    self.rise_tolerance=rise_tolerance
    self.dot_min=dot_min
    self.dot_min_single=dot_min_single
    self.diffs_single=None
    self.span=span # 3.5 for helices, 2 for strands, 0 for other
    self.target_rise=target_rise # 1.54 for helices, 3.3 strands
    self.optimal_delta_length=optimal_delta_length # target change 
        # in number of residues for this
        # segment (based on segment type and end-to-end distance)
    self.hierarchy=hierarchy # optional full hierarchy of this segment

    self.sites=None
    if sites:
      self.sites=sites
    elif self.hierarchy:
      self.get_sites_from_hierarchy()
    if start_resno is not None:
      self.start_resno=start_resno
    elif self.hierarchy:
      self.start_resno=get_first_resno(self.hierarchy)
      assert start_resno is None
    else:
      start_resno=1
      self.start_resno=start_resno

  def summary(self):
    return self.show_summary()

  def show_summary(self,return_text=False):
    text="Class: %6s  Length: %d  Start: %d  End: %d" %(
     self.segment_type,self.length(),self.get_start_resno(),
      self.get_start_resno()+self.length()-1)
    if return_text:
      return text
    else:
      print >>self.out, text

  def trim_ends(self,start_pos=None,end_pos=None):
    # trim back from start_pos to end_pos (not residue number, position in
    #   existing segment
    start_res=self.get_start_resno()+start_pos
    end_res=self.get_start_resno()+end_pos
    if start_pos==0 and end_pos==self.length()-1:
      return # nothing to do
    atom_selection="resseq %s:%s" %(resseq_encode(start_res),
       resseq_encode(end_res))

    self.hierarchy=apply_atom_selection(atom_selection,hierarchy=self.hierarchy)
    self.start_resno=self.start_resno+start_pos
    self.get_sites_from_hierarchy()
    if self.optimal_delta_length:  
      self.optimal_delta_length=0 # we no longer know how long it should be

  def contains_ends(self,first_resno_of_chain=None,last_resno_of_chain=None):
     if first_resno_of_chain>=self.get_start_resno() and \
        first_resno_of_chain<=self.get_end_resno():
       contains_left_end_of_chain=True
     else:
       contains_left_end_of_chain=False
     if last_resno_of_chain>=self.get_start_resno() and \
        last_resno_of_chain<=self.get_end_resno():
       contains_right_end_of_chain=True
     else:
       contains_right_end_of_chain=False
     return contains_left_end_of_chain,contains_right_end_of_chain

  def get_sites_from_hierarchy(self):
    atom_selection="name ca"
    sele=apply_atom_selection(atom_selection,hierarchy=self.hierarchy)

    # extract coordinates
    self.sites=sele.extract_xray_structure().sites_cart()


  def get_start_resno(self):
    return self.start_resno

  def get_end_resno(self):
    return self.start_resno+self.length()-1

  def length(self):
    return self.sites.size()

  def optimal_delta_length(self):
    return self.optimal_delta_length

  def get_norms(self):
    diffs,norms=self.get_diffs_and_norms()
    return self.norms

  def get_rise(self):
    # report mean distance from i to i+span-1
    if not self.span: return 0.
    return self.get_norms().min_max_mean().mean/self.span

  def get_cosine(self):
    self.mean_dot_single=None
    # report mean normalized dot product of diffs with mean direction
    diffs,norms=self.get_diffs_and_norms()
    if not diffs: return 0.
    average_direction=get_average_direction(diffs=diffs)
    cosines=flex.double()
    for j in xrange(len(diffs)):
      dot=col(diffs[j]).dot(average_direction)
      cosines.append(dot)
    mean_dot=cosines.min_max_mean().mean

    # now get mean cross product and dot with average_direction...
    diffs_single=self.get_diffs_single()
    cosines=flex.double()
    for j in xrange(len(diffs_single)):
      dot=col(diffs_single[j]).dot(average_direction)
      cosines.append(dot)
    self.mean_dot_single=cosines.min_max_mean().mean

    return mean_dot

  def get_diffs_single(self):
    if self.diffs_single:
      return self.diffs_single
    else:
      sites_offset_1=self.sites[1:]
      self.diffs_single=sites_offset_1-self.sites[:-1]
      norms=self.diffs_single.norms()
      self.diffs_single=self.diffs_single/norms
    return self.diffs_single

  def get_orientation_points(self,start_res=None,end_res=None): # for strand
    # get 2 points along strand axis and two point at a CA position at end
    #    These can be used to superimpose strands
    #   do not depend on the details of the strands
    n=self.length()
    if n < self.minimum_length: return None
    if start_res is None: start_res=self.get_start_resno()
    if end_res is None: end_res=start_res+n-1

    if 2*( (end_res-start_res+1)//2)==end_res-start_res+1: # even
      # take 2 residues centered at middle of chain
      i1=(start_res+end_res)//2 #
      center_middle=self.get_centroid(start_res=i1,end_res=i1+1)
    else:  # odd number of residues average just before and after
      im=(start_res+end_res)//2 # middle residue
      cm1=self.get_centroid(start_res=im-1,end_res=im)
      cm2=self.get_centroid(start_res=im,end_res=im+1)
      center_middle=0.5*(matrix.col(cm1)+matrix.col(cm2))

    center_1=self.get_centroid(start_res=start_res,end_res=start_res+1)
    center_2=self.get_centroid(start_res=end_res-1,end_res=end_res)
    ca_1=self.get_site(resno=start_res)
    ca_2=self.get_site(resno=end_res)
    orientation_points=flex.vec3_double((center_1,center_2,ca_1,ca_2,
      center_middle))
    return orientation_points

  def get_targets(self):
    target=self.target_rise*self.span
    tol=self.rise_tolerance
    dot_min=self.dot_min
    tol=tol*self.span # rise_tolerance * span
    minimum_length=self.minimum_length
    return target,tol,dot_min,minimum_length


  def is_ok(self):
    # check rise and cosine and dot product of single pairs with average
    #  direction and compare to target
    rise=self.get_rise()
    dot=self.get_cosine()
    dot_single=self.mean_dot_single # set by get_cosine()
    if self.minimum_length is not None and self.length()<self.minimum_length:
      return False
    elif self.dot_min is not None and dot < self.dot_min:
      return False
    elif self.target_rise is not None and self.rise_tolerance is not None and \
       (rise < self.target_rise-self.rise_tolerance or \
           rise > self.target_rise+self.rise_tolerance):
      return False
    elif self.dot_min_single is not None \
       and dot_single < self.dot_min_single:
      return False
    else:
      return True

  def get_lsq_fit(self,other=None,start_res=None,end_res=None):
    # superimpose other segment (sucha as a helix ) on this one using
    #  orientation points
    # start_res,end_res are start end in other to match
    # The lengths do not have to match

    # the orientation points for helix are center of
    #  N-term of helix, center of C-term of
    #  helix, and CA of first and last residues in the helix.

    self_orientation_points=self.get_orientation_points()
    other_orientation_points=other.get_orientation_points(
      start_res=start_res,end_res=end_res)

    if self_orientation_points.size()!=other_orientation_points.size():
      return None # did not match size

    if not self_orientation_points or not other_orientation_points:
      return None # could not make a fit

    # now superimpose other on to self...
    lsq_fit_obj=superpose.least_squares_fit(
       reference_sites=self_orientation_points,
       other_sites=other_orientation_points)

    # check it
    new_xyz=flex.vec3_double()
    for x in other_orientation_points:
      new_xyz.append(lsq_fit_obj.r * matrix.col(x) + lsq_fit_obj.t)
    delta=self_orientation_points-new_xyz
    lsq_fit_obj.rms=delta.norms().min_max_mean().mean

    return lsq_fit_obj

  def apply_lsq_fit(self,lsq_fit_obj=None,hierarchy=None,
     start_res=None,end_res=None):
    atom_selection="resseq %s:%s" %(resseq_encode(start_res),
       resseq_encode(end_res))
    sele=apply_atom_selection(atom_selection,hierarchy=hierarchy)

    asc=hierarchy.atom_selection_cache()
    sel = asc.selection(string = atom_selection)
    alt_hierarchy = hierarchy.deep_copy().select(sel)
    str1=sele.as_pdb_string()
    str2=alt_hierarchy.as_pdb_string()
    assert str1== str2
    sele=alt_hierarchy # ZZZ
    assert sele.overall_counts().n_residues  # should not be None
    from cStringIO import StringIO
    f=StringIO()
    for atom in sele.atoms_with_labels():
      new_xyz=lsq_fit_obj.r * matrix.col(atom.xyz) + lsq_fit_obj.t
      atom.set_xyz(new_xyz)
      print >>f,atom.format_atom_record()
    return get_pdb_hierarchy(text=f.getvalue())

  def get_site(self,resno=None):
    first_residue=self.get_start_resno()
    return self.sites[resno-first_residue]

  def get_sites(self,start_res=None,end_res=None):
    if start_res is None or end_res is None:
      return self.sites
    else:
      first_residue=self.get_start_resno()
      return self.sites[start_res-first_residue:end_res-first_residue+1]

  def get_centroid(self,start_res=None,end_res=None):
    # get centroid of residues from start_res to end_res
    sites=self.get_sites(start_res=start_res,end_res=end_res)
    return sites.mean()

class helix(segment): # Methods specific to helices

  def __init__(self,params=None,sites=None,start_resno=None,hierarchy=None,
     is_n_terminus=None,
     is_c_terminus=None,
     frequency=None,
     base_score=None,
     optimal_delta_length=None,
     verbose=None,
     out=sys.stdout):
    self.setup(sites=sites,start_resno=start_resno,hierarchy=hierarchy,
     segment_type='helix',
     minimum_length=params.minimum_length,
     buffer_residues=params.buffer_residues,
     standard_length=params.standard_length,
     frequency=frequency,
     base_score=params.base_score,
     is_n_terminus=is_n_terminus,
     is_c_terminus=is_c_terminus,
     span=params.span,
     target_rise=params.rise,
     rise_tolerance=params.rise_tolerance,
     dot_min=params.dot_min,
     dot_min_single=params.dot_min_single,
     optimal_delta_length=optimal_delta_length,
     verbose=verbose,
     out=out)

  def get_diffs_and_norms(self):
    assert self.span==3.5
    # report distances between i and average of i+3,i+4
    if not hasattr(self,'diffs'):
      sites_offset_3=self.sites[3:-1]
      sites_offset_4=self.sites[4:]
      average_offset=0.5*(sites_offset_3+sites_offset_4)
      self.diffs=average_offset-self.sites[:-4]
      self.norms=self.diffs.norms()
      self.diffs=self.diffs/self.diffs.norms()
    return self.diffs,self.norms


  def get_orientation_points(self,start_res=None,end_res=None):
    # get 3 points along helix axis, two CA at ends
    #   These can be used to superimpose helices and
    #   do not depend on the details of the helix
    n=self.length()
    if n < self.minimum_length: return None
    if start_res is None: start_res=self.get_start_resno()
    if end_res is None: end_res=start_res+n-1

    if 2*( (end_res-start_res+1)//2)==end_res-start_res+1: # even
      # take 4 residues centered at middle of chain
      i1=(start_res+end_res)//2-1 #
      center_middle=self.get_centroid(start_res=i1,end_res=i1+3)
    else:  # odd number of residues average just before and after
      im=(start_res+end_res)//2 # middle residue
      cm1=self.get_centroid(start_res=im-2,end_res=im+1)
      cm2=self.get_centroid(start_res=im-1,end_res=im+2)
      center_middle=0.5*(matrix.col(cm1)+matrix.col(cm2))

    center_1=self.get_centroid(start_res=start_res,end_res=start_res+3)
    center_2=self.get_centroid(start_res=end_res-3,end_res=end_res)
    ca_1=self.get_site(resno=start_res)
    ca_2=self.get_site(resno=end_res)
    orientation_points=flex.vec3_double(
      (center_1,center_2,ca_1,ca_2,center_middle))
    return orientation_points

class strand(segment):

  # Methods specific to strands

  def __init__(self,params=None,sites=None,start_resno=None,hierarchy=None,
     is_n_terminus=None,
     is_c_terminus=None,
     frequency=None,
     base_score=None,
     optimal_delta_length=None,
     verbose=None,
     out=sys.stdout):

    self.setup(sites=sites,start_resno=start_resno,hierarchy=hierarchy,
      segment_type='strand',
      span=params.span,
      minimum_length=params.minimum_length,
      buffer_residues=params.buffer_residues,
      standard_length=params.standard_length,
      frequency=frequency,
      base_score=params.base_score,
      is_n_terminus=is_n_terminus,
      is_c_terminus=is_c_terminus,
      target_rise=params.rise,
      rise_tolerance=params.rise_tolerance,
      dot_min=params.dot_min,
      dot_min_single=params.dot_min_single,
      optimal_delta_length=optimal_delta_length,
      verbose=verbose,
      out=out)

  def get_diffs_and_norms(self):
    # report distances between i and i+2
    assert self.span==2
    if not hasattr(self,'diffs'):
      sites_offset_2=self.sites[2:]
      self.diffs=sites_offset_2-self.sites[:-2]
      self.norms=self.diffs.norms()
      self.diffs=self.diffs/self.norms

    return self.diffs,self.norms

  def get_orientation_points(self,start_res=None,end_res=None): # for strand
    # get 2 points along strand axis and two point at a CA position at end
    #    These can be used to superimpose strands
    #   do not depend on the details of the strands
    n=self.length()
    if n < self.minimum_length: return None
    if start_res is None: start_res=self.get_start_resno()
    if end_res is None: end_res=start_res+n-1

    if 2*( (end_res-start_res+1)//2)==end_res-start_res+1: # even
      # take 2 residues centered at middle of chain
      i1=(start_res+end_res)//2 #
      center_middle=self.get_centroid(start_res=i1,end_res=i1+1)
    else:  # odd number of residues average just before and after
      im=(start_res+end_res)//2 # middle residue
      cm1=self.get_centroid(start_res=im-1,end_res=im)
      cm2=self.get_centroid(start_res=im,end_res=im+1)
      center_middle=0.5*(matrix.col(cm1)+matrix.col(cm2))

    center_1=self.get_centroid(start_res=start_res,end_res=start_res+1)
    center_2=self.get_centroid(start_res=end_res-1,end_res=end_res)
    ca_1=self.get_site(resno=start_res)
    ca_2=self.get_site(resno=end_res)
    orientation_points=flex.vec3_double((center_1,center_2,ca_1,ca_2,
      center_middle))
    return orientation_points

class other(segment):

  # Methods specific to other (not strand, not helix). Looks like strand mostly

  def __init__(self,params=None,sites=None,start_resno=None,hierarchy=None,
     is_n_terminus=None,
     is_c_terminus=None,
     frequency=None,
     base_score=None,
     optimal_delta_length=None,
     verbose=None,
     out=sys.stdout):

    self.setup(sites=sites,start_resno=start_resno,hierarchy=hierarchy,
      segment_type='other',
      span=params.span,
      minimum_length=params.minimum_length,
      buffer_residues=params.buffer_residues,
      standard_length=params.standard_length,
      base_score=params.base_score,
      is_n_terminus=is_n_terminus,
      is_c_terminus=is_c_terminus,
      target_rise=params.rise,
      rise_tolerance=params.rise_tolerance,
      dot_min=params.dot_min,
      dot_min_single=params.dot_min_single,
      frequency=frequency,
      optimal_delta_length=optimal_delta_length,
      verbose=verbose,
      out=out)

  def get_diffs_and_norms(self):
    self.diffs=None
    self.norms=None
    return self.diffs,self.norms

  def get_orientation_points(self,start_res=None,end_res=None): # for other
    # just use all points
    if not self.orientation_points or \
        not start_res==self.orientation_points_start or \
        not end_res==self.orientation_points_end: # calculate it
      n=self.length()
      if n < self.minimum_length:
         self.orientation_points=None
      else:
        first_residue=self.get_start_resno()
        self.orientation_points=flex.vec3_double(self.get_sites(
          start_res=start_res,end_res=end_res))
        self.orientation_points_start=start_res
        self.orientation_points_end=end_res
    return self.orientation_points

class find_segment: # class to look for a type of segment

  def setup(self,params=None,model=None,segment_type='helix',
      extract_segments_from_pdb=None,
      make_unique=None,
      verbose=None,
      out=sys.stdout):
    # Assumes model is just 1 chain of sequential residues
    #   obtained with split_model

    self.out=out
    self.extract_segments_from_pdb=extract_segments_from_pdb
    self.make_unique=make_unique
    self.verbose=verbose
    self.params=params
    self.model=model

    self.segment_type=segment_type
    if self.segment_type=='helix':
      self.segment_class=helix
      self.allow_insertions=params.allow_insertions
      self.allow_deletions=params.allow_deletions
    elif self.segment_type=='strand':
      self.segment_class=strand
      self.allow_insertions=params.allow_insertions
      self.allow_deletions=params.allow_deletions
    elif self.segment_type=='other':
      self.segment_class=other
      self.allow_insertions=params.allow_insertions
      self.allow_deletions=params.allow_deletions

    h=self.segment_class(params=params) # contains numbers we need
    self.target_rise=h.target_rise
    self.span=h.span
    self.buffer_residues=h.buffer_residues
    self.standard_length=h.standard_length

    # number of residues in segment compared to number of diffs
    # (4 more than diffs for helix; 2 more for strand)
    if self.span is not None:
      self.last_residue_offset=int(self.span+0.99)
    else:
      self.last_residue_offset=0

    # set start residue number if not set.
    self.start_resno=get_first_resno(self.model.hierarchy)

    self.segments=[]

    # select ca atoms
    sites=self.get_sites() # just a list of sites (CA atoms for entire chain)
    if not sites:
      return # nothing to do

    # get list of difference vectors i
    #  Then get list of segments in segment_dict
    diffs,norms,self.segment_dict=self.find_segments(params=params,sites=sites)

    # NOTE: separate find_segments method for each kind of segment
    # figure out how many residues really should be in these segments
    optimal_delta_length_dict,norm_dict=self.get_optimal_lengths(
       segment_dict=self.segment_dict,norms=norms)

    # create segment object (helix,strand) for each segment
    # If longer than the standard length, make a series of shorter segments

    keys=self.segment_dict.keys()
    keys.sort()

    # Specify which residues are in each segment
    #    (self.last_residue_offset past start)
    # and offset to start with first residue number

    # Possibly cut into pieces of length self.standard_length plus
    #   buffer_residues on each end for total length of
    #   standard_length+2*buffer_residues

    # If number of residues is to be changed, use just full length of the
    # segment (do not cut) and do not add buffer residues

    # NOTE: buffered on end with n_buffer residues unless
    #   at very ends of chain in which case set is_n_terminus or is_c_terminus

    start_end_list=[]
    # NOTE: self.segment_class is a substitution for helix/strand classes
    for i in keys:
      # this segment starts at i and ends at self.segment_dict[i]
      overall_start_res=i
      overall_end_res=self.segment_dict[i]+self.last_residue_offset
      overall_length=overall_end_res-overall_start_res+1
      optimal_delta_length=optimal_delta_length_dict[overall_start_res]
      if (optimal_delta_length > 0 and not self.allow_insertions) or \
         (optimal_delta_length < 0 and not self.allow_deletions):
        optimal_delta_length=0
      # add to self.segments
      # cut into pieces of size self.standard_length
      if overall_length<=self.standard_length or \
          self.extract_segments_from_pdb is not None: 
        start_end_list.append([overall_start_res,overall_end_res])
      else:
        for start_res_use in xrange(
           overall_start_res,
           overall_start_res+overall_length-self.standard_length+1):
          start_end_list.append(
             [start_res_use,start_res_use+self.standard_length-1])

    self.segment_dict={}
    for start_res_use,end_res_use in start_end_list:
        self.segment_dict[start_res_use]=end_res_use
        ok=self.extract_segment(params=params,
         start_res=start_res_use,end_res=end_res_use,sites=sites,
         optimal_delta_length=optimal_delta_length)
        if not ok:
          del self.segment_dict[start_res_use]

  def show_summary(self,out=None):
    if not out: out=self.out
    for h in self.segments:
      print >>out,\
       "Class: %6s  N: %d Start: %d End: %d " %(
          self.segment_type,h.get_end_resno()-h.get_start_resno()+1,
           h.get_start_resno(),h.get_end_resno(),
         ) +" Rise:%5.2f A Dot:%5.2f" %(
          h.get_rise(),h.get_cosine())

  def get_used_residues_list(self):
    # just return a list of used residues
    used_residues=[]
    for i in self.segment_dict.keys():
      for j in xrange(i,self.segment_dict[i]):
        if not j in used_residues:
          used_residues.append(j)
    return used_residues

  def set_used_residues(self,used_residues):
    for i in self.segment_dict.keys():
      for j in xrange(i,self.segment_dict[i]+1):
        used_residues[j]=True
    return used_residues

  def extract_segment(self,params=None,start_res=None,end_res=None,sites=None,
       optimal_delta_length=None):

      start_res_with_buffer=start_res  # -self.buffer_residues
      end_res_with_buffer=end_res   #  +self.buffer_residues

      if start_res_with_buffer<0: # at the n-terminus
        start_res_with_buffer=start_res
        is_n_terminus=True
      else:
        is_n_terminus=False
      if end_res_with_buffer>len(sites)-1: # at c-term
        end_res_with_buffer=end_res
        is_c_terminus=True
      else:
        is_c_terminus=False

      assert len(sites)>=end_res  # make sure we are in bounds

      # get the hierarchy if necessary
      if self.extract_segments_from_pdb:
        start_resno=start_res_with_buffer+self.start_resno
        end_resno=end_res_with_buffer+self.start_resno
        atom_selection="resseq %s:%s" %(resseq_encode(start_resno),
           resseq_encode(end_resno))
        hierarchy=apply_atom_selection(
           atom_selection,hierarchy=self.model.hierarchy)
      else:
        hierarchy=None

      h=self.segment_class(params=params,
        sites=sites[start_res_with_buffer:end_res_with_buffer+1],
        start_resno=start_res_with_buffer+self.start_resno,
        optimal_delta_length=optimal_delta_length,hierarchy=hierarchy,
        is_n_terminus=is_n_terminus,
        is_c_terminus=is_c_terminus,
          )
      if h.is_ok():
        self.segments.append(h)
        return True
      else:
        return False

  def get_sites(self):
    atom_selection="name ca"
    sele=apply_atom_selection(atom_selection,hierarchy=self.model.hierarchy)
    if not sele.overall_counts().n_residues:
      return []
    else:
      # extract coordinates
      sites=sele.extract_xray_structure().sites_cart()
      return sites

  def get_segments(self):
    return self.segments

  def find_segments(self,params=None,sites=None,hierarchy=None):
    # set up segment class for this kind of segment
    # for example, helix(sites=sites)
    h=self.segment_class(params=params,sites=sites,hierarchy=hierarchy)

    # get difference vectors i to i+2 (strands) or i to avg of i+3/i+4 (helix)
    diffs,norms=h.get_diffs_and_norms()  # difference vectors and lengths


    target,tol,dot_min,minimum_length=h.get_targets()

    # now find sequence of N residues where diffs are parallel and values of
    # diffs are all about the same and about span*rise = 5.4 A for helix

    n=len(diffs)
    segment_dict={}
    used_residues=[]

    if self.find_alpha:  # a find_alpha object is here if we want to exclude
      previously_used_residues=self.find_alpha.get_used_residues_list()
    else:
      previously_used_residues=[]

    for i in xrange(n):
      if i in used_residues or abs(norms[i]-target)>tol: continue
      if i in previously_used_residues: continue
      # i is start of segment
      segment_dict[i]=i  # lists end of segment
      used_residues.append(i)
      for j in xrange(i+1,n):
        if abs(norms[j]-target)>tol: break
        if j in previously_used_residues: break
        dot=col(diffs[j]).dot(col(diffs[j-1]))
        if dot < dot_min: break
        segment_dict[i]=j
        used_residues.append(j)
    # skip if only as long as a single span (requires 2 to be sure)
    for i in segment_dict.keys():
      segment_length=segment_dict[i]+1+self.last_residue_offset-i
      if segment_length<minimum_length:
        del segment_dict[i]

    # prune out anything not really in this segment type based on direction

    still_changing=True
    n_cycles=0
    max_cycles=2
    while still_changing and n_cycles < max_cycles:
      still_changing=False
      new_segment_dict={}
      keys=segment_dict.keys()
      keys.sort()
      for i in keys:
        average_direction=get_average_direction(
          diffs=diffs,i=i,j=segment_dict[i])
        segment_start=None
        for j in xrange(i,segment_dict[i]+1):
          dot=col(diffs[j]).dot(average_direction)
          if dot >= dot_min:
            if segment_start is None:
              segment_start=j
            new_segment_dict[segment_start]=j
          else:
            segment_start=None
            still_changing=True
      segment_dict=new_segment_dict
      for i in segment_dict.keys():
        segment_length=segment_dict[i]+1+self.last_residue_offset-i
        if segment_length<minimum_length:
          del segment_dict[i] # not long enough
      n_cycles+=1

    # merge any segments
    found=True
    n_cycles=0
    while found and n_cycles <=10:
      found=False
      n_cycles+=1
      keys=segment_dict.keys()
      keys.sort()  # sorted on start
      for i1,i2 in zip(keys[:-1],keys[1:]):
        if found: break
        if i2 <= segment_dict[i1]+self.last_residue_offset+1:
          # could merge. Test it
          h=self.segment_class(params=params,
             sites=sites[i1:segment_dict[i2]+self.last_residue_offset+1],
             start_resno=i1+self.start_resno,
             hierarchy=hierarchy)
          if h.is_ok():
           found=True
           segment_dict[i1]=segment_dict[i2]
           del segment_dict[i2]




    return diffs,norms,segment_dict

  def get_optimal_lengths(self,segment_dict=None,norms=None):

    keys=segment_dict.keys()
    keys.sort()

    # get average distance i to i+3/i+4 for this segment (to see if it is
    #   stretched out as can happen at lowres)
    # or for strand, distance i to i+2

    optimal_delta_length_dict={}
    norm_dict={}
    standard_rise=self.target_rise # 5.5 A for 3.5 residues
    for i in keys:
      norms_used=norms[i:segment_dict[i]+1]
      mean_norm=norms_used.min_max_mean().mean
      # guess number of residues it should be: standard rise
      segment_rise=mean_norm/self.span  # we used mean of i+3 and i+4 for helix
      segment_length=segment_dict[i]+1+self.last_residue_offset-i
      segment_distance=segment_rise*(segment_length-1)
      optimal_length=int(0.5+segment_distance/standard_rise)+1
      optimal_delta_length_dict[i]=optimal_length-segment_length
      norm_dict[i]=mean_norm
      if self.verbose:
        print >>self.out,\
       "Chain start: %d  end: %d  N: %d  Rise: %7.2f  Optimal length: %d" %(
         i+self.start_resno,i+segment_length-1+self.start_resno,
         segment_length,segment_rise,optimal_length)

    return optimal_delta_length_dict,norm_dict

class find_alpha_helix(find_segment):

  def __init__(self,params=None,model=None,verbose=None,
    extract_segments_from_pdb=None,
    make_unique=None,
    out=sys.stdout):
    self.find_alpha=None
    self.find_beta=None

    self.setup(params=params,model=model,segment_type='helix',
     extract_segments_from_pdb=extract_segments_from_pdb,
     make_unique=make_unique,
     verbose=verbose,out=out)

class find_beta_strand(find_segment):

  def __init__(self,params=None,model=None,verbose=None,
      extract_segments_from_pdb=None,
      make_unique=None,
      find_alpha=None,
      out=sys.stdout):

    self.find_beta=None
    self.find_alpha=find_alpha

    self.setup(params=params,
      model=model,segment_type='strand',
      extract_segments_from_pdb=extract_segments_from_pdb,
      make_unique=make_unique,
      verbose=verbose,out=out)

class find_other_structure(find_segment):

  def __init__(self,find_alpha=None,find_beta=None,
      params=None,model=None,
      extract_segments_from_pdb=None,
      make_unique=None,
      verbose=None,out=sys.stdout):

    self.find_alpha=find_alpha
    self.find_beta=find_beta

    self.setup(params=params,model=model,segment_type='other',
      extract_segments_from_pdb=extract_segments_from_pdb,
      make_unique=make_unique,
      verbose=verbose,out=out)

  def get_optimal_lengths(self,segment_dict=None,norms=None):
    # always zero
    keys=segment_dict.keys()
    optimal_delta_length_dict={}
    norm_dict={}
    for i in keys:
      optimal_delta_length_dict[i]=0
      norm_dict[i]=None
    return optimal_delta_length_dict,norm_dict

  def add_start_end_to_segment_dict(self,
       params,n=None,n_buf=None,segment_start=None,
       segment_end=None,segment_dict=None):
    if self.make_unique:
      nn=0
    else:
      nn=n_buf
    start_pos=max(0,segment_start-nn)
    end_pos=min(segment_end+nn,n-1)
    if start_pos==0 and end_pos<params.minimum_length-1:
      end_pos=min(params.minimum_length-1,end_pos)
    if end_pos==n-1 and start_pos > n-params.minimum_length:
      start_pos=max(0,n-params.minimum_length)
    segment_dict[start_pos]=end_pos
    return segment_dict

  def find_segments(self,params=None,sites=None,hierarchy=None):
    # find everything that is not alpha and not beta, put buffer_residues
    #  buffer on the end of each one.

    # set up segment class for this kind of segment SPECIFIC FOR OTHER
    # for example, helix(sites=sites)
    h=self.segment_class(params=params,sites=sites,hierarchy=hierarchy)

    # cross off used residues except for buffer of buffer_residues
    n_buf=params.buffer_residues
    n=len(sites)
    used_residues=n*[False]

    if self.find_alpha:  # a find_alpha object is here if we want to exclude
      used_residues=self.find_alpha.set_used_residues(used_residues)

    if self.find_beta:  # a find_beta object is here if we want to exclude
      used_residues=self.find_beta.set_used_residues(used_residues)
    segment_dict={}
    segment_start=None
    segment_end=None
    for i,used in zip(xrange(n),used_residues):
      if not used: # use it
        segment_end=i
        if segment_start is None:
          segment_start=i
      else:
        if segment_start is not None: # save it
          segment_dict=self.add_start_end_to_segment_dict(
             params,n=n,n_buf=n_buf,segment_start=segment_start,
             segment_end=segment_end,segment_dict=segment_dict)
        segment_start=None
        segment_end=None

    if segment_start is not None: # save it
      segment_dict=self.add_start_end_to_segment_dict(
             params,n=n,n_buf=n_buf,segment_start=segment_start,
             segment_end=segment_end,segment_dict=segment_dict)

    return None,None,segment_dict

class find_secondary_structure: # class to look for secondary structure

  def __init__(self,params=None,args=None,hierarchy=None,models=None,
      verbose=None,out=sys.stdout):

    if not args: args=[]

    if not params:  # get params from args if necessary
      params=self.get_params(args,out=out)

    if models:
      self.models=models
    else:
      if not hierarchy: # get a hierarchy if necessary
        if not params.input_files.pdb_in or \
            not os.path.isfile(params.input_files.pdb_in):
          raise Sorry("Missing file: %s" %(params.input_files.pdb_in))
        hierarchy=get_pdb_hierarchy(text=open(params.input_files.pdb_in).read())
      else: # make a copy as we are going to modify hierarchy
        hierarchy=hierarchy.deep_copy()
   
      # remove alt conformers
      from mmtbx.pdbtools import remove_alt_confs
      remove_alt_confs(hierarchy,always_keep_one_conformer=True)

      self.models=split_model(hierarchy=hierarchy)

    for model in self.models:
      print "\nModel %d  N: %d  Start: %d End; %d" %(
          model.info['chain_number'],
          model.length(),model.first_residue(),model.last_residue())
      model=self.find_ss_in_model(params=params,model=model,out=out)
      if model.find_alpha and verbose:
        model.find_alpha.show_summary()
      if model.find_beta and verbose:
        model.find_beta.show_summary()
      if model.find_other and verbose:
        model.find_other.show_summary()

    # now each model has a find_alpha, find_beta, find_other object

  def find_ss_in_model(self,params=None,model=None,out=sys.stdout):
    if params.find_ss_structure.find_alpha:
      model.find_alpha=find_alpha_helix(params=params.alpha,
        model=model,verbose=params.control.verbose,
        make_unique=params.find_ss_structure.make_unique,
        extract_segments_from_pdb=params.extract_segments_from_pdb.extract,
        out=out)

    self.beta_list_by_model=[]
    if params.find_ss_structure.find_beta:
      if params.find_ss_structure.exclude_alpha_in_beta:
        find_alpha=model.find_alpha
      else:
        find_alpha=None
      model.find_beta=find_beta_strand(params=params.beta,
        model=model,verbose=params.control.verbose,
        make_unique=params.find_ss_structure.make_unique,
        extract_segments_from_pdb=params.extract_segments_from_pdb.extract,
        find_alpha=find_alpha,
        out=out)

    if params.find_ss_structure.make_unique: # do this before finding other
      # get unique residues in alpha, beta if desired
      self.make_unique(model)

    if params.find_ss_structure.find_other:
      if params.find_ss_structure.exclude_alpha_beta_in_other:
        find_alpha=model.find_alpha
        find_beta=model.find_beta
      else:
        find_alpha=None
        find_beta=None
      model.find_other=find_other_structure(params=params.other,
        model=model,verbose=params.control.verbose,
        make_unique=params.find_ss_structure.make_unique,
        extract_segments_from_pdb=params.extract_segments_from_pdb.extract,
        find_alpha=find_alpha,find_beta=find_beta,
        out=out)
    return model

  def make_unique(self,model):
    # iteratively work down residues in this model starting with helix, then
    #  strand, then anything left
    # trim in from end if overlapping with previous, delete to nearest end
    # if a used residue in the middle

    is_used_list=model.length()*[None]
    first_res=model.first_residue()

    for find_group in ['find_alpha','find_beta','find_other']:
      fg=getattr(model,find_group)
      if not fg: continue
      new_segment_list=[]
      for s in fg.segments:
        start_resno=s.get_start_resno()
        first_pos=start_resno-first_res
        last_pos=s.get_end_resno()-first_res
	already_used=is_used_list[first_pos:last_pos+1]
        
        if not True in already_used:  # cross of used ones
           new_start=0
           new_end=len(already_used)-1
        else:
          # find range of ok positions...
          new_start,new_end=self.get_start_end(already_used=already_used)
          if new_start is None: # delete it
            continue
          s.trim_ends(start_pos=new_start,end_pos=new_end)
          if not s.is_ok():
            continue

        # mark used residues
        for i in xrange(first_pos+new_start,first_pos+new_end+1):
          is_used_list[i]=True
          
        #save it
        new_segment_list.append(s) # keep it if we get this far
      fg.segments=new_segment_list


  def get_start_end(self,already_used=None):
    # goes from first available to end of available (not necessarily optimal)
    new_start=None
    new_end=None 
    for i in xrange(len(already_used)):
      if already_used[i]: 
        if new_end is not None:
          return new_start,new_end
        else:
          pass # have not yet started
      else:
        if new_end is None:
          new_start=i
          new_end=i
        else:
          new_end=i
    return new_start,new_end

  def get_params(self,args,out=sys.stdout):
    command_line = iotbx.phil.process_command_line_with_files(
      args=args, master_phil=master_phil)
    params = command_line.work.extract()
    print >>out,"\nFind secondary structure in hierarchy"
    master_phil.format(python_object=params).show(out=out)
    return params


def tst_01():
  text="""
ATOM      1  N   GLY U  11     -20.099   8.864  10.230  1.00 52.86           N
ATOM      2  CA  GLY U  11     -19.400   7.571   9.992  1.00 52.86           C
ATOM      3  C   GLY U  11     -17.867   8.054   9.511  1.00 52.86           C
ATOM      4  O   GLY U  11     -18.033   8.077   8.293  1.00 52.86           O
ATOM      5  N   GLY U  12     -16.673   8.172  10.087  1.00 52.86           N
ATOM      6  CA  GLY U  12     -15.517   7.978   9.123  1.00 52.86           C
ATOM      7  C   GLY U  12     -15.338   6.468   9.675  1.00 52.86           C
ATOM      8  O   GLY U  12     -15.070   5.904  10.738  1.00 52.86           O
ATOM      9  N   GLY U  13     -15.479   5.785   8.520  1.00 52.86           N
ATOM     10  CA  GLY U  13     -14.900   4.879   9.412  1.00 52.86           C
ATOM     11  C   GLY U  13     -14.424   3.514   9.380  1.00 52.86           C
ATOM     12  O   GLY U  13     -13.660   2.944   8.602  1.00 52.86           O
ATOM     13  N   GLY U  14     -14.286   3.425  10.715  1.00 52.86           N
ATOM     14  CA  GLY U  14     -13.338   2.264  10.393  1.00 52.86           C
ATOM     15  C   GLY U  14     -12.091   1.361   9.704  1.00 52.86           C
ATOM     16  O   GLY U  14     -11.343   0.974   8.806  1.00 52.86           O
ATOM     17  N   GLY U  15     -12.859   0.538  10.412  1.00 52.86           N
ATOM     18  CA  GLY U  15     -12.030  -0.556  10.373  1.00 52.86           C
ATOM     19  C   GLY U  15     -12.456  -1.530   9.072  1.00 52.86           C
ATOM     20  O   GLY U  15     -13.665  -1.682   8.859  1.00 52.86           O
ATOM     21  N   GLY U  16     -11.909  -3.009   7.572  1.00 42.55           N
ATOM     22  CA  GLY U  16     -11.031  -3.180   6.040  1.00 42.55           C
ATOM     23  C   GLY U  16     -10.265  -3.837   5.718  1.00 42.55           C
ATOM     24  O   GLY U  16      -9.515  -4.299   6.577  1.00 42.55           O
ATOM     25  N   GLY U  17     -11.247  -3.944   5.165  1.00 52.86           N
ATOM     26  CA  GLY U  17     -10.273  -5.222   4.456  1.00 52.86           C
ATOM     27  C   GLY U  17      -9.735  -4.585   3.332  1.00 52.86           C
ATOM     28  O   GLY U  17     -10.280  -3.765   2.596  1.00 52.86           O
ATOM     29  N   GLY U  18      -8.496  -4.996   3.135  1.00 52.86           N
ATOM     30  CA  GLY U  18      -8.015  -4.822   2.186  1.00 52.86           C
ATOM     31  C   GLY U  18      -7.002  -5.187   1.126  1.00 52.86           C
ATOM     32  O   GLY U  18      -6.107  -5.649   1.834  1.00 52.86           O
ATOM     33  N   GLY U  19      -7.061  -5.398  -0.184  1.00 52.86           N
ATOM     34  CA  GLY U  19      -6.009  -6.212  -0.948  1.00 52.86           C
ATOM     35  C   GLY U  19      -4.828  -5.254  -1.458  1.00 52.86           C
ATOM     36  O   GLY U  19      -5.131  -4.184  -1.991  1.00 52.86           O
ATOM     37  N   GLY U  20      -3.582  -5.721  -1.443  1.00 52.86           N
ATOM     38  CA  GLY U  20      -2.729  -4.680  -2.268  1.00 52.86           C
ATOM     39  C   GLY U  20      -2.086  -4.815  -3.421  1.00 52.86           C
ATOM     40  O   GLY U  20      -1.108  -5.401  -3.871  1.00 52.86           O
ATOM     41  N   GLY U  21      -2.755  -3.911  -4.145  1.00 52.86           N
ATOM     42  CA  GLY U  21      -2.425  -3.693  -5.384  1.00 52.86           C
ATOM     43  C   GLY U  21      -0.789  -3.646  -5.564  1.00 52.86           C
ATOM     44  O   GLY U  21      -0.196  -4.136  -6.539  1.00 52.86           O
"""
  print "\nFinding beta strands"
  import iotbx.pdb
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(text)).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,verbose=True,out=sys.stdout)


def tst_02():
  text="""
ATOM      1  N   GLY A   1      42.375 -12.180  24.780  1.00 35.31
ATOM      2  CA  GLY A   1      43.603 -11.488  24.325  1.00 35.57
ATOM      3  C   GLY A   1      43.288 -10.171  23.615  1.00 34.64
ATOM      4  O   GLY A   1      42.111  -9.896  23.277  1.00 35.82
ATOM      5  N   ILE A   2      44.323  -9.391  23.299  1.00 32.23
ATOM      6  CA  ILE A   2      44.200  -8.183  22.475  1.00 27.55
ATOM      7  C   ILE A   2      43.750  -8.629  21.119  1.00 24.92
ATOM      8  O   ILE A   2      43.068  -7.904  20.409  1.00 23.73
ATOM      9  CB  ILE A   2      45.525  -7.320  22.425  1.00 30.10
ATOM     10  CG1 ILE A   2      45.924  -6.837  23.820  1.00 29.64
ATOM     11  CG2 ILE A   2      45.555  -6.173  21.386  1.00 30.54
ATOM     12  CD1 ILE A   2      44.837  -6.338  24.762  1.00 32.44
ATOM     13  N   GLY A   3      44.161  -9.867  20.749  1.00 22.69
ATOM     14  CA  GLY A   3      43.999 -10.264  19.329  1.00 21.05
ATOM     15  C   GLY A   3      42.433 -10.405  19.166  1.00 22.08
ATOM     16  O   GLY A   3      41.912 -10.061  18.096  1.00 22.86
ATOM     17  N   ALA A   4      41.862 -10.961  20.191  1.00 21.60
ATOM     18  CA  ALA A   4      40.378 -11.260  20.106  1.00 21.80
ATOM     19  C   ALA A   4      39.584  -9.950  20.087  1.00 21.67
ATOM     20  O   ALA A   4      38.676  -9.747  19.278  1.00 22.21
ATOM     21  CB  ALA A   4      40.061 -12.080  21.350  1.00 22.97
ATOM     22  N   VAL A   5      39.936  -9.001  20.956  1.00 22.13
ATOM     23  CA  VAL A   5      39.355  -7.658  21.083  1.00 19.34
ATOM     24  C   VAL A   5      39.536  -6.896  19.795  1.00 18.81
ATOM     25  O   VAL A   5      38.626  -6.314  19.126  1.00 17.22
ATOM     26  CB  VAL A   5      39.843  -6.933  22.338  1.00 20.40
ATOM     27  CG1 VAL A   5      39.237  -5.519  22.413  1.00 23.06
ATOM     28  CG2 VAL A   5      39.745  -7.587  23.653  1.00 21.67
ATOM     29  N   LEU A   6      40.752  -7.021  19.222  1.00 16.05
ATOM     30  CA  LEU A   6      41.062  -6.432  17.957  1.00 17.59
ATOM     31  C   LEU A   6      40.230  -6.935  16.870  1.00 20.01
ATOM     32  O   LEU A   6      39.649  -6.121  16.029  1.00 21.90
ATOM     33  CB  LEU A   6      42.627  -6.461  17.880  1.00 24.58
ATOM     34  CG  LEU A   6      43.125  -6.023  16.524  1.00 23.91
ATOM     35  CD1 LEU A   6      42.706  -4.584  16.210  1.00 27.44
ATOM     36  CD2 LEU A   6      44.669  -6.152  16.638  1.00 29.31
ATOM     37  N   LYS A   7      39.981  -8.229  16.721  1.00 19.83
ATOM     38  CA  LYS A   7      39.079  -8.646  15.636  1.00 22.55
ATOM     39  C   LYS A   7      37.648  -8.063  15.784  1.00 19.04
ATOM     40  O   LYS A   7      37.031  -7.839  14.731  1.00 21.18
ATOM     41  CB  LYS A   7      38.854 -10.176  15.616  1.00 27.62
ATOM     42  CG  LYS A   7      40.011 -10.993  15.144  1.00 40.15
ATOM     43  CD  LYS A   7      39.691 -12.487  15.325  1.00 47.84
ATOM     44  CE  LYS A   7      40.599 -13.394  14.493  1.00 53.11
ATOM     45  NZ  LYS A   7      39.966 -14.755  14.319  1.00 55.47
ATOM     46  N   VAL A   8      37.111  -7.988  16.981  1.00 19.69
ATOM     47  CA  VAL A   8      35.792  -7.369  17.211  1.00 20.52
ATOM     48  C   VAL A   8      35.776  -5.881  16.885  1.00 20.31
ATOM     49  O   VAL A   8      34.775  -5.402  16.257  1.00 20.12
ATOM     50  CB  VAL A   8      35.113  -7.600  18.562  1.00 23.09
ATOM     51  CG1 VAL A   8      34.774  -9.045  18.851  1.00 24.22
ATOM     52  CG2 VAL A   8      35.769  -6.970  19.726  1.00 27.95
ATOM     68  N   THR A  11      35.315  -5.870  13.113  1.00 21.12
ATOM     69  CA  THR A  11      34.132  -6.405  12.343  1.00 24.14
ATOM     70  C   THR A  11      32.874  -6.299  13.140  1.00 24.13
ATOM     71  O   THR A  11      31.823  -5.889  12.630  1.00 27.78
ATOM     72  CB  THR A  11      34.462  -7.925  11.935  1.00 27.78
ATOM     73  OG1 THR A  11      34.702  -8.591  13.199  1.00 31.22
ATOM     74  CG2 THR A  11      35.695  -7.959  11.047  1.00 30.31
ATOM     75  N   GLY A  12      32.918  -6.606  14.420  1.00 23.00
ATOM     76  CA  GLY A  12      31.584  -6.595  15.140  1.00 24.17
ATOM     77  C   GLY A  12      31.264  -5.168  15.584  1.00 24.45
ATOM     78  O   GLY A  12      30.060  -4.914  15.603  1.00 23.79
ATOM     79  N   LEU A  13      32.195  -4.276  15.948  1.00 22.05
ATOM     80  CA  LEU A  13      31.923  -2.919  16.364  1.00 23.24
ATOM     81  C   LEU A  13      31.251  -2.004  15.324  1.00 21.24
ATOM     82  O   LEU A  13      30.232  -1.308  15.679  1.00 23.02
ATOM     83  CB  LEU A  13      33.146  -2.183  17.008  1.00 25.55
ATOM     84  CG  LEU A  13      32.913  -1.523  18.351  1.00 27.01
ATOM     85  CD1 LEU A  13      33.999  -0.517  18.760  1.00 27.53
ATOM     86  CD2 LEU A  13      31.587  -0.842  18.531  1.00 24.54
ATOM     87  N   PRO A  14      31.689  -1.979  14.106  1.00 17.38
ATOM     88  CA  PRO A  14      31.026  -1.278  13.030  1.00 17.52
ATOM     89  C   PRO A  14      29.521  -1.670  12.857  1.00 18.98
ATOM     90  O   PRO A  14      28.657  -0.801  12.744  1.00 16.75
ATOM     91  CB  PRO A  14      31.845  -1.614  11.816  1.00 17.80
ATOM     92  CG  PRO A  14      33.118  -2.205  12.303  1.00 18.05
ATOM     93  CD  PRO A  14      33.098  -2.363  13.769  1.00 17.96
ATOM     94  N   ALA A  15      29.207  -2.952  12.868  1.00 16.39
ATOM     95  CA  ALA A  15      27.822  -3.418  12.724  1.00 17.10
ATOM     96  C   ALA A  15      27.023  -3.016  13.951  1.00 16.98
ATOM     97  O   ALA A  15      25.872  -2.551  13.769  1.00 16.78
ATOM     98  CB  ALA A  15      27.741  -4.906  12.502  1.00 19.58
ATOM     99  N   LEU A  16      27.570  -3.117  15.127  1.00 15.97
ATOM    100  CA  LEU A  16      26.958  -2.649  16.351  1.00 18.20
ATOM    101  C   LEU A  16      26.614  -1.169  16.344  1.00 20.28
ATOM    102  O   LEU A  16      25.599  -0.734  16.933  1.00 18.32
ATOM    103  CB  LEU A  16      27.811  -3.027  17.542  1.00 19.70
ATOM    104  CG  LEU A  16      27.384  -2.550  18.921  1.00 22.23
ATOM    105  CD1 LEU A  16      26.031  -3.234  19.257  1.00 27.80
ATOM    106  CD2 LEU A  16      28.445  -2.970  19.933  1.00 21.91
ATOM    107  N   ILE A  17      27.514  -0.365  15.791  1.00 20.97
ATOM    108  CA  ILE A  17      27.343   1.056  15.618  1.00 20.41
ATOM    109  C   ILE A  17      26.081   1.392  14.758  1.00 18.17
ATOM    110  O   ILE A  17      25.380   2.240  15.282  1.00 16.46
ATOM    111  CB  ILE A  17      28.579   1.847  15.132  1.00 21.10
ATOM    112  CG1 ILE A  17      29.586   1.858  16.352  1.00 25.66
ATOM    113  CG2 ILE A  17      28.268   3.288  14.691  1.00 22.04
ATOM    114  CD1 ILE A  17      30.856   2.696  16.161  1.00 27.00
ATOM    115  N   SER A  18      25.930   0.759  13.657  1.00 16.97
ATOM    116  CA  SER A  18      24.825   0.827  12.744  1.00 19.98
ATOM    117  C   SER A  18      23.499   0.405  13.438  1.00 18.89
ATOM    118  O   SER A  18      22.557   1.165  13.352  1.00 18.37
ATOM    119  CB  SER A  18      25.076   0.039  11.491  1.00 20.39
ATOM    120  OG  SER A  18      23.902   0.046  10.670  1.00 23.87
ATOM    121  N   TRP A  19      23.512  -0.661  14.161  1.00 17.71
ATOM    122  CA  TRP A  19      22.492  -1.085  15.081  1.00 15.72
ATOM    123  C   TRP A  19      22.083   0.004  16.012  1.00 18.02
ATOM    124  O   TRP A  19      20.820   0.244  16.160  1.00 16.93
ATOM    125  CB  TRP A  19      22.854  -2.410  15.767  1.00 15.59
ATOM    126  CG  TRP A  19      21.803  -2.993  16.678  1.00 17.94
ATOM    127  CD1 TRP A  19      20.917  -3.950  16.210  1.00 18.03
ATOM    128  CD2 TRP A  19      21.448  -2.745  18.041  1.00 16.03
ATOM    129  NE1 TRP A  19      20.060  -4.304  17.222  1.00 21.28
ATOM    130  CE2 TRP A  19      20.357  -3.624  18.372  1.00 20.45
ATOM    131  CE3 TRP A  19      21.879  -1.892  19.048  1.00 14.41
ATOM    132  CZ2 TRP A  19      19.784  -3.690  19.612  1.00 17.15
ATOM    133  CZ3 TRP A  19      21.292  -1.950  20.288  1.00 17.24
ATOM    134  CH2 TRP A  19      20.230  -2.805  20.601  1.00 15.13
ATOM    135  N   ILE A  20      22.930   0.594  16.823  1.00 14.82
ATOM    136  CA  ILE A  20      22.628   1.633  17.766  1.00 15.67
ATOM    137  C   ILE A  20      21.917   2.819  17.080  1.00 17.51
ATOM    138  O   ILE A  20      20.942   3.365  17.655  1.00 17.70
ATOM    139  CB  ILE A  20      23.902   2.192  18.499  1.00 16.25
ATOM    140  CG1 ILE A  20      24.481   0.986  19.363  1.00 15.50
ATOM    141  CG2 ILE A  20      23.599   3.421  19.390  1.00 14.54
ATOM    142  CD1 ILE A  20      26.033   1.304  19.637  1.00 18.18
ATOM    143  N   LYS A  21      22.464   3.177  15.957  1.00 16.61
ATOM    144  CA  LYS A  21      21.888   4.236  15.157  1.00 19.84
ATOM    145  C   LYS A  21      20.436   3.910  14.752  1.00 21.02
ATOM    146  O   LYS A  21      19.685   4.899  14.971  1.00 22.80
ATOM    147  CB  LYS A  21      22.699   4.646  13.935  1.00 16.73
ATOM    148  CG  LYS A  21      23.944   5.416  14.471  1.00 22.19
ATOM    149  CD  LYS A  21      24.919   5.669  13.300  1.00 25.86
ATOM    150  CE  LYS A  21      26.173   6.287  13.908  1.00 32.91
ATOM    151  NZ  LYS A  21      27.199   6.564  12.863  1.00 39.11
ATOM    152  N   ARG A  22      20.075   2.728  14.351  1.00 19.18
ATOM    153  CA  ARG A  22      18.740   2.273  14.020  1.00 20.38
ATOM    154  C   ARG A  22      17.807   2.365  15.177  1.00 22.46
ATOM    155  O   ARG A  22      16.648   2.899  15.096  1.00 23.24
ATOM    156  CB  ARG A  22      18.733   0.860  13.418  1.00 22.73
ATOM    157  CG  ARG A  22      19.309   0.956  12.004  1.00 22.29
ATOM    158  CD  ARG A  22      19.117  -0.300  11.254  1.00 25.78
ATOM    159  NE  ARG A  22      19.382  -1.562  11.991  1.00 28.94
ATOM    160  CZ  ARG A  22      20.624  -2.104  11.889  1.00 33.77
ATOM    161  NH1 ARG A  22      21.664  -1.554  11.252  1.00 31.21
ATOM    162  NH2 ARG A  22      20.870  -3.271  12.498  1.00 33.11
ATOM    163  N   LYS A  23      18.257   1.937  16.323  1.00 20.49
ATOM    164  CA  LYS A  23      17.500   1.928  17.550  1.00 22.62
ATOM    165  C   LYS A  23      17.216   3.361  18.100  1.00 25.47
ATOM    166  O   LYS A  23      16.204   3.554  18.811  1.00 24.62
ATOM    167  CB  LYS A  23      18.257   1.128  18.589  1.00 24.16
ATOM    168  CG  LYS A  23      17.979  -0.388  18.463  1.00 31.03
ATOM    169  CD  LYS A  23      16.858  -0.657  19.514  1.00 38.52
ATOM    170  CE  LYS A  23      16.197  -1.986  19.278  1.00 44.05
ATOM    171  NZ  LYS A  23      15.412  -2.493  20.477  1.00 48.30
ATOM    172  N   ARG A  24      18.155   4.268  17.844  1.00 23.99
ATOM    173  CA  ARG A  24      18.059   5.674  18.276  1.00 27.11
ATOM    174  C   ARG A  24      16.996   6.355  17.441  1.00 28.34
ATOM    175  O   ARG A  24      16.224   7.166  18.029  1.00 31.82
ATOM    176  CB  ARG A  24      19.446   6.379  18.188  1.00 22.24
ATOM    177  CG  ARG A  24      20.182   6.236  19.542  1.00 20.01
ATOM    178  CD  ARG A  24      21.577   6.742  19.433  1.00 25.00
ATOM    179  NE  ARG A  24      21.715   8.115  18.950  1.00 24.24
ATOM    180  CZ  ARG A  24      21.745   9.136  19.837  1.00 22.72
ATOM    181  NH1 ARG A  24      21.508   8.848  21.086  1.00 24.10
ATOM    182  NH2 ARG A  24      22.134  10.375  19.512  1.00 24.09
ATOM    183  N   GLN A  25      16.889   6.080  16.192  1.00 29.90
ATOM    184  CA  GLN A  25      15.836   6.730  15.339  1.00 37.50
ATOM    185  C   GLN A  25      14.463   6.162  15.554  1.00 39.92
ATOM    186  O   GLN A  25      13.435   6.799  15.227  1.00 42.39
ATOM    187  CB  GLN A  25      16.283   6.713  13.874  1.00 40.82
ATOM    188  CG  GLN A  25      17.607   7.505  13.800  1.00 46.91
ATOM    189  CD  GLN A  25      18.413   7.205  12.564  1.00 50.89
ATOM    190  OE1 GLN A  25      19.596   7.559  12.477  1.00 54.23
ATOM    191  NE2 GLN A  25      17.757   6.512  11.624  1.00 51.61
ATOM    192  N   GLN A  26      14.363   4.945  16.064  1.00 43.06
ATOM    193  CA  GLN A  26      13.132   4.360  16.583  1.00 46.66
ATOM    194  C   GLN A  26      12.618   5.071  17.837  1.00 48.00
ATOM    195  O   GLN A  26      11.409   5.344  17.597  1.00 50.85
ATOM    196  CB  GLN A  26      13.144   2.874  16.862  1.00 49.14
ATOM    197  CG  GLN A  26      13.166   1.972  15.625  1.00 54.49
ATOM    198  CD  GLN A  26      13.506   0.551  16.116  1.00 58.10
ATOM    199  OE1 GLN A  26      13.083   0.180  17.239  1.00 59.27
ATOM    200  NE2 GLN A  26      14.304  -0.117  15.290  1.00 57.61
"""
  print "\nFinding helices"
  import iotbx.pdb
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(text)).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,verbose=True,out=sys.stdout)

if __name__=="__main__":
  args=sys.argv[1:]
  if 'tst' in args:
    tst_01()
    tst_02()
  else:
    fss=find_secondary_structure(args=args,verbose=True,out=sys.stdout)
