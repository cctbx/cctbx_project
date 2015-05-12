from __future__ import division

# find_ss_from_ca.py
# a tool to find helices, strands, non-helices/strands in segments of
#  structure

from iotbx.pdb import resseq_encode
import iotbx.phil
import os,sys
from libtbx.utils import Sorry
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

  find_ss_structure {  # note values from regularize_from_pdb overwrite these

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

     write_helix_sheet_records = True
       .type = bool
       .help = Write HELIX and SHEET records
       .short_caption = Write HELIX/SHEET records

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

def get_first_residue(hierarchy):
  if not hierarchy:
    return None
  for model in hierarchy.models():
    for chain in model.chains():
      for conformer in chain.conformers():
        for residue in conformer.residues():
          return residue

def get_last_residue(hierarchy):
  if not hierarchy:
    return None
  last_residue=None
  for model in hierarchy.models():
    for chain in model.chains():
      for conformer in chain.conformers():
        for residue in conformer.residues():
          last_residue=residue
  return last_residue

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

  def pdb_records(self,last_id=0):
    from iotbx.pdb import secondary_structure

    records=[]
    k=last_id
    for s in self.segments:
      if not s.hierarchy: continue
      start=get_first_residue(s.hierarchy)
      end=get_last_residue(s.hierarchy)
      chain_id=get_chain_id(s.hierarchy)
      k=k+1
      record = secondary_structure.pdb_helix(
        serial=k,
        helix_id=k,
        start_resname=start.resname,
        start_chain_id=chain_id,
        start_resseq=start.resseq_as_int(),
        start_icode=start.icode,
        end_resname=end.resname,
        end_chain_id=chain_id,
        end_resseq=end.resseq_as_int(),
        end_icode=end.icode,
        helix_class=1, # alpha helix only
        comment="",
        length=s.length())
      records.append(record)
    return records

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

  def pdb_records(self,last_id=0):
    from iotbx.pdb import secondary_structure

    records=[]
    k=last_id
    for s in self.segments:
      if not s.hierarchy: continue
      start=get_first_residue(s.hierarchy)
      end=get_last_residue(s.hierarchy)
      chain_id=get_chain_id(s.hierarchy)
      k=k+1
      current_sheet = secondary_structure.pdb_sheet(
        sheet_id=k,
        n_strands=1, # XXX TODO figure out others and registration
        strands=[],
        registrations=[])
      first_strand = secondary_structure.pdb_strand(
        sheet_id=k,
        strand_id=1,
        start_resname=start.resname,
        start_chain_id=chain_id,
        start_resseq=start.resseq_as_int(),
        start_icode=start.icode,
        end_resname=end.resname,
        end_chain_id=chain_id,
        end_resseq=end.resseq_as_int(),
        end_icode=end.icode,
        sense=0)
      current_sheet.add_strand(first_strand)
      current_sheet.add_registration(None)
      records.append(current_sheet)
    return records

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

  def pdb_records(self,last_id=0):
    return []

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
      self.find_ss_in_model(params=params,model=model,out=out)

    if params.find_ss_structure.write_helix_sheet_records:
      self.find_sheets(out=out) # organize strands into sheets
      self.write_pdb_records(out=out)

  def show_summary(self,verbose=None,out=sys.stdout):

    for model in self.models:
      print >>out,"\nModel %d  N: %d  Start: %d End: %d" %(
          model.info['chain_number'],
          model.length(),model.first_residue(),model.last_residue())
      if verbose:
        if model.find_alpha:
          model.find_alpha.show_summary(out=out)
        if model.find_beta:
          model.find_beta.show_summary(out=out)
        if model.find_other:
          model.find_other.show_summary(out=out)
    if verbose:
      print >>out,"\nPDB RECORDS:"
      print >>out,self.all_pdb_records

  def find_sheets(self,out=sys.stdout,max_sheet_ca_ca_dist=6.):
    all_strands=[]
    for model in self.models:
      all_strands+=model.find_beta.segments
    if not all_strands: return
    print >>out,"\nFinding sheets from %d strands" %(len(all_strands))
    self.get_strand_pairs(all_strands=all_strands,tol=max_sheet_ca_ca_dist)
    # self.pair_dict is list of all the strands that each strand matches with
    # self.info_dict is information on a particular pair of strands:
    #  self.info_dict["%d:%d" %(i,j)]=
    #     [first_ca_1,last_ca_1,first_ca_2,last_ca_2,is_parallel]

    for i in self.pair_dict.keys():
      for j in self.pair_dict[i]:
        key="%d:%d" %(i,j)
        if key in self.info_dict.keys():
          [first_ca_1,last_ca_1,first_ca_2,last_ca_2,is_parallel]=self.info_dict[key]
          print "PAIR:  Seg %d:  %d:%d   Seg %d:  %d:%d  Parallel: %s" %(
            i,first_ca_1,last_ca_1,j,first_ca_2,last_ca_2,is_parallel)

    # Create sheets from paired strands.
    self.used_strands=[] # keep track of which ones we have assigned
    self.all_strands=all_strands

    # single_strands are those with no matching strands
    single_strands=self.get_strands_by_pairs(pairs=0)
    pair_strands=self.get_strands_by_pairs(pairs=1)
    triple_strands=self.get_strands_by_pairs(pairs=2)
    multiple_strands=self.get_strands_by_pairs(pairs=None)

    self.used_strands=[] # initialize again
    self.used_strands+=single_strands
    # find all sheets with edges (beginning with a paired strand)
    for i in pair_strands:
      if i in self.used_strands:continue
      self.used_strands.append(i)
      print "Building sheet starting with %d" %(i)
      strand_list=[i]
      current_strand=i
      while current_strand is not None:
        current_strand=self.get_available_strand(current_strand=current_strand)
        if current_strand is not None: strand_list.append(current_strand)
      print "Strand list:",strand_list

    print "single strands:",single_strands
    print "pair strands:",pair_strands
    print "triple strands:",triple_strands
    print "multiple strands:",multiple_strands

  def get_available_strand(self,current_strand=None):
    return None

  def get_strands_by_pairs(self,pairs=None):
    strand_list=[]
    while 1:  # get all single strands
      i=self.get_unused_strand(n=len(self.all_strands),
         used_strands=self.used_strands,pairs=pairs)
      if i is None: break
      self.used_strands.append(i)
      strand_list.append(i)
    return strand_list

  def get_unused_strand(self,n=None,used_strands=None,pairs=None):
    for i in xrange(n):
      if i in used_strands: continue
      if pairs is None or len(self.pair_dict.get(i,[]))==pairs:
        return i
    return None

  def get_strand_pairs(self,all_strands=None,tol=None):
    self.info_dict={}
    self.pair_dict={}
    for i in xrange(len(all_strands)):
      self.pair_dict[i]=[]
    for i in xrange(len(all_strands)):
      for j in xrange(i+1,len(all_strands)):
        self.ca1=None
        self.ca2=None
        if self.ca_pair_is_close(all_strands[i],all_strands[j],
            tol=tol):

          # figure out alignment and whether it really is ok
          first_last_1_and_2=self.align_strands(
            all_strands[i],all_strands[j],tol=tol)

          if first_last_1_and_2:
            # we have a match
            [first_ca_1,last_ca_1,first_ca_2,last_ca_2,is_parallel]=\
               first_last_1_and_2
            self.pair_dict[i].append(j)
            self.info_dict["%d:%d" %(i,j)]=\
               [first_ca_1,last_ca_1,first_ca_2,last_ca_2,is_parallel]
            # and make an entry for the other way around
            self.pair_dict[j].append(i)
            self.info_dict["%d:%d" %(j,i)]=\
               [first_ca_2,last_ca_2,first_ca_1,last_ca_1,is_parallel]

  def align_strands(self,s1,s2,tol=None):
    # figure out best alignment and directions. Require at least 2 residues
    sites1=s1.get_sites()
    sites2=s2.get_sites()
    sites2_reversed=sites2.deep_copy().reversed()

    # self.ca1 and self.ca2 are pos of closest residues from ca_pair_is_close
    best_offset=None
    best_reverse=None
    best_keep_1=None
    best_keep_2=None
    best_score=None
    for offset in [-1,0,1]: # just in case an offset from best guess is better
      if self.ca1+offset < 0 or self.ca1+offset > len(sites1)-1: continue
      dd_list,keep_1,keep_2=self.get_residue_pairs_in_sheet(sites1,sites2,
       center1=self.ca1+offset,center2=self.ca2,tol=tol)
      dd_list_reverse,keep_1_reverse,keep_2_reverse=\
         self.get_residue_pairs_in_sheet(sites1,sites2_reversed,
         center1=self.ca1+offset,center2=len(sites2)-self.ca2-1,tol=tol)

      if len(keep_1)<2 and len(keep_1_reverse)<2:
        continue
      elif len(keep_1_reverse)>len(keep_1):
        score=len(keep_1_reverse)-0.001*flex.double(dd_list_reverse).norm()
        use_reverse=True
        if best_score is None or score>best_score:
          best_keep_1=keep_1_reverse
          best_keep_2=keep_2_reverse
          best_reverse=True
          best_offset=offset
          best_score=score
      else:
        score=len(keep_1)-0.001*flex.double(dd_list).norm()
        use_reverse=False
        if best_score is None or score>best_score:
          best_keep_1=keep_1
          best_keep_2=keep_2
          best_reverse=False
          best_offset=offset
          best_score=score

    if not best_score:
      return None

    first_ca_1=best_keep_1[0]
    last_ca_1=best_keep_1[-1]
    if best_reverse: # reversed
      first_ca_2=len(sites2)-(best_keep_2[-1]+1)
      last_ca_2=len(sites2)-(best_keep_2[0]+1)
      is_parallel=False
    else:  # forward
      first_ca_2=best_keep_2[0]
      last_ca_2=best_keep_2[-1]
      is_parallel=True

    return [first_ca_1,last_ca_1,first_ca_2,last_ca_2,is_parallel]


  def get_residue_pairs_in_sheet(self,sites1,sites2,
     center1=None,center2=None,tol=None):
    # figure out pairs that are within tol. Center pairs are center1-center2
    dd=(col(sites1[center1])-col(sites2[center2])).norm_sq()
    keep1_list=[]
    keep2_list=[]
    dd_list=[]
    if dd<=tol**2: # plausible at least
      start_offset=max(-center1,-center2)
      end_offset=min(len(sites1)-(center1+1),len(sites2)-(center2+1))
      for offset in xrange(start_offset,end_offset+1):
        i1=center1+offset
        i2=center2+offset
        dd=(col(sites1[i1])-col(sites2[i2])).norm_sq()
        if dd <= tol**2:
          keep1_list.append(i1)
          keep2_list.append(i2)
          dd_list.append(dd**0.5)
    return dd_list,keep1_list,keep2_list

  def ca_pair_is_close(self,s1,s2,tol=None,
      dist_per_residue=3.5,jump=4):
    best_dist_sq=None
    self.ca1=None
    self.ca2=None
    sites1=s1.get_sites()
    sites2=s2.get_sites()
    while jump > 0:
      # see if it is even close
      for i in xrange(jump//2,len(sites1),jump):
        for j in xrange(jump//2,len(sites2),jump):
          dd=(col(sites1[i])-col(sites2[j])).norm_sq()
          if best_dist_sq is None or dd < best_dist_sq:
            best_dist_sq=dd
            self.ca1=i
            self.ca2=j
      if best_dist_sq is None or \
         best_dist_sq**0.5 > (jump+1)*dist_per_residue+tol: # no hope
           break
      if jump==1:
           break
      jump=max(1,jump//2)  # try finer search
    if best_dist_sq is not None and best_dist_sq <= tol**2:
      return True
    else:
      return False

  def write_pdb_records(self,out=sys.stdout):
    # save everything as pdb_helix or pdb_sheet objects
    self.pdb_helix_list=[]
    self.pdb_sheet_list=[]
    for model in self.models:
      self.pdb_helix_list+=model.find_alpha.pdb_records(
        last_id=len(self.pdb_helix_list))
      self.pdb_sheet_list+=model.find_beta.pdb_records(
        last_id=len(self.pdb_sheet_list))

    from cStringIO import StringIO
    all_pdb_records=StringIO()
    for helix in self.pdb_helix_list:
      print >>all_pdb_records,helix.as_pdb_str()
    for sheet in self.pdb_sheet_list:
      print >>all_pdb_records,sheet.as_pdb_str()
    self.all_pdb_records=all_pdb_records.getvalue()

    self.show_summary(verbose=True,out=out)

  def get_all_pdb_records(self):
    return self.all_pdb_records

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
      args=args,
      master_phil=master_phil,
      pdb_file_def="input_files.pdb_in")

    params = command_line.work.extract()
    print >>out,"\nFind secondary structure in hierarchy"
    master_phil.format(python_object=params).show(out=out)
    return params


if __name__=="__main__":
  args=sys.argv[1:]
  fss=find_secondary_structure(args=args,verbose=True,out=sys.stdout)
