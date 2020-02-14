from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.find_ss_from_ca

# find_ss_from_ca.py
# a tool to find helices, strands, non-helices/strands in segments of
#  structure
#


from libtbx import adopt_init_args
from iotbx.pdb import resseq_encode
import iotbx.phil
import os,sys
from libtbx.utils import Sorry
from scitbx.matrix import col
from scitbx.math import superpose, matrix
from scitbx.array_family import flex
from copy import deepcopy
from iotbx.pdb import secondary_structure
from six.moves import zip
from six.moves import range
from six.moves import cStringIO as StringIO

master_phil = iotbx.phil.parse("""

  input_files {
    pdb_in = None
      .type = path
      .help = Input PDB file
      .short_caption = Input PDB file

    secondary_structure_input = None
      .type = path
      .help = Optional input secondary-structure file (can be a PDB or just \
              text) with secondary-structure (HELIX/SHEET) records. If \
              supplied, this secondary structure is used as a starting point \
              and additional information is added if possible.  If \
              force_ss_in=True, then this is used exactly as input. \
      .short_caption = Input secondary structure

    force_secondary_structure_input = False
      .type = bool
      .help = Force use of input secondary_structure without changes (even \
              if H-bonds are not withing max_h_bond_length)
      .short_caption = Force input secondary structure

  }
  output_files {
    pdb_records_file = None
      .type = path
      .help = Output file with HELIX/SHEET records
      .short_caption = Output HELIX/SHEET records file
  }

  find_ss_structure {  # note values from regularize_from_pdb overwrite these
                       # Also values for ss_by_chain and
                       #   max_rmsd are set in
                       #   mmtbx/secondary_structure/__init__.py
     ss_by_chain = True
       .type = bool
       .help = Find secondary structure only within individual chains. \
               Alternative is to allow H-bonds between chains. Can be \
               much slower with ss_by_chain=False. If your model is complete \
               use ss_by_chain=True. If your model is many fragments, use \
               ss_by_chain=False.
       .short_caption = Secondary structure by chain
       .expert_level = 1

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
               as the same for ss identification
       .short_caption = Maximum rmsd
       .expert_level = 3
     use_representative_chains = True
       .type = bool
       .help = Use a representative of all chains with the same sequence. \
               Alternative is to examine each chain individually. Can be \
               much slower with use_representative_of_chain=False if there \
               are many symmetry copies. Ignored unless ss_by_chain is True.
       .short_caption = Use representative chains
       .expert_level = 3
     max_representative_chains = 100
       .type = float
       .help = Maximum number of representative chains
       .short_caption = Maximum representative chains
       .expert_level = 3
     find_alpha = True
       .type = bool
       .help = Find alpha helices
       .short_caption = Find alpha helices

     helices_are_alpha = False
       .type = bool
       .help = Find alpha helices and not three_ten or pi
       .short_caption = Helices are alpha

     find_three_ten = True
       .type = bool
       .help = Find three_ten helices
       .short_caption = Find three_ten helices

     find_pi = True
       .type = bool
       .help = Find pi helices
       .short_caption = Find pi helices

     find_beta = True
       .type = bool
       .help = Find beta structure
       .short_caption = Find beta structure

     find_other = False
       .type = bool
       .help = Find other structure
       .short_caption = Find other structure

     exclude_alpha_in_beta  = False
       .type = bool
       .help = Exclude regions already identified as alpha from three_ten, pi,\
               and beta
       .short_caption = Exclude alpha from beta

     make_unique = True
       .type = bool
       .help = Assign each residue to a unique type of structure
       .short_caption = Assign residues to unique structure

     cut_up_segments = False
       .type = bool
       .help = Cut up segments (make short segments of secondary structure)
       .short_caption = Cut up segments

     extend_segments = False
       .type = bool
       .help = Try to extend segments in both directions one residue at a time
       .short_caption = Extend segments

     set_up_helices_sheets = True
       .type = bool
       .help = Set up HELIX and SHEET records
       .short_caption = Set up HELIX/SHEET records

     write_helix_sheet_records = True
       .type = bool
       .help = Write HELIX and SHEET records
       .short_caption = Write HELIX/SHEET records

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
              structure if require_h_bonds is set
      .short_caption = Minimum number of H bonds

    maximum_poor_h_bonds = None
      .type = int
      .help = Maximum number of poor hydrogen bonds to keep secondary \
              structure element (helix/sheet) if require_h_bonds is set. \
              Note: None means ignore this test, 0 means allow no poor H-bonds.
      .short_caption = Maximum number of poor H bonds

  }

  alpha {
    include scope mmtbx.secondary_structure.secondary_structure_params.alpha_params
  }

  three_ten {
    include scope mmtbx.secondary_structure.secondary_structure_params.three_ten_params
  }

  pi {
    include scope mmtbx.secondary_structure.secondary_structure_params.pi_params
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

def is_close_to(r,last_r,distance_cutoff=None,use_default_distance_cutoff=True):
  if not r or not last_r:
    return None
  r_atom=None
  atom_type=None
  last_atom_type=None
  for atom in r.atoms():
    if atom.name.strip() in ["CA","P"]:
       r_atom=atom
       atom_type=atom.name.strip()
  last_r_atom=None
  for atom in last_r.atoms():
    if atom.name.strip() in ["CA","P"]:
       last_r_atom=atom
       last_atom_type=atom.name.strip()
  if not r_atom or not last_r_atom: return None
  if last_atom_type != atom_type: return None # change of type
  if distance_cutoff is None and use_default_distance_cutoff:
     if atom_type=="P":
        distance_cutoff=10.  # pretty sure these are not connected
     else:
        distance_cutoff=6.

  dd=col(r_atom.xyz)-col(last_r_atom.xyz)
  if dd.length()<distance_cutoff:
    return True
  else:
    return False

def make_four_char_unique_chain_id(id,used_chain_ids=None):
  # Make a unique chain id that looks like id but has extra characters as
  #  necessary.  Use first 2 chars as input part, last 2 chars of 4 as unique
  original_id=id
  id=id.strip()
  if len(id)<1:
    id="XX"
  elif len(id)<2:
    id="%sX" %(id)
  elif len(id)==2:
    pass # already good
  else:
    raise Sorry(
      "Cannot use longer-than-2 chain ids in make_four_char_unique_chain_id")
  new_chars=" abcdefghijklmnopqrstuvwxyz"
  new_chars+="abcdefghijklmnopqrstuvwxyz".upper()
  new_chars+="0123456789"
  for a in new_chars[1:]:
    for b in new_chars:
      new_id=id+a+b.strip()
      if not new_id in used_chain_ids:
        used_chain_ids.append(new_id)
        return new_id,used_chain_ids
  raise Sorry("Unable to generate a new unique chain ID for %s" %(original_id))

def split_model(model=None,hierarchy=None,verbose=False,info=None,
     only_first_model=None,distance_cutoff=None,
     use_default_distance_cutoff=True,
     out=sys.stdout):
  # XXX NOTE: this splits model at all icode residues (one model per residue)
  # The routine extract_segment below assumes that the residues in an individual
  #  model are sequential (no insertion codes)
  # if CA-CA or P-P distance is > distance-cutoff then split there
  model_list=[]
  if hierarchy:
    if not info: info={}
  elif hasattr(model,'hierarchy'): # a model object with hierarchy and info
    hierarchy=model.hierarchy
    info=model.info
  else:
    return []  # nothing here

  if not hierarchy or hierarchy.overall_counts().n_residues<1:
    return [] # nothing to do
  total_models=0
  base_chain_ids=[]
  for m in hierarchy.models():
    total_models+=1
    if total_models>1:
      raise Sorry("Sorry, find_ss_from_ca cannot use multi-model files. "+\
        "Please use phenix.pdbtools to select just one model from your file.")
    for chain in m.chains():
      new_hierarchy=iotbx.pdb.input(
         source_info="Model", lines=flex.split_lines("")).construct_hierarchy()
      mm=iotbx.pdb.hierarchy.model()
      cc=iotbx.pdb.hierarchy.chain()
      cc.id=chain.id  # copy chain ID
      new_hierarchy.append_model(mm)
      mm.append_chain(cc)

      last_resseq=None
      last_r=None
      is_linked=None
      for r in chain.residue_groups():
        if distance_cutoff is not None or use_default_distance_cutoff:
          is_linked=is_close_to(r,last_r,distance_cutoff=distance_cutoff,
             use_default_distance_cutoff=use_default_distance_cutoff)

        if (last_resseq is not None )  and (r.resseq_as_int()!=last_resseq+1
           or (
            (distance_cutoff is not None or use_default_distance_cutoff) and
            (not is_linked) )):

          # save and make new model
          new_model_info=model_info(hierarchy=new_hierarchy,info=deepcopy(info))
          model_list.append(new_model_info)
          new_model_info.info['chain_number']=len(model_list)

          # and make a new one
          new_hierarchy=iotbx.pdb.input(
             source_info="Model",
             lines=flex.split_lines("")).construct_hierarchy()
          mm=iotbx.pdb.hierarchy.model()
          cc=iotbx.pdb.hierarchy.chain()
          cc.id=chain.id  # copy chain ID
          new_hierarchy.append_model(mm)
          mm.append_chain(cc)
          last_resseq=None
          last_r=None
        # add on a residue...
        cc.append_residue_group(r.detached_copy())
        last_resseq=r.resseq_as_int()
        last_r=r
      new_hierarchy.reset_atom_i_seqs()
      new_model_info=model_info(hierarchy=new_hierarchy,info=deepcopy(info))
      model_list.append(new_model_info)
      new_model_info.info['chain_number']=len(model_list)
    if only_first_model:
      break

  if verbose:
    print("Models after splitting:", file=out)
    for m in model_list:
      print("Chain: %d  Residues: %d" %(
        m.info.get('chain_number'),
        m.hierarchy.overall_counts().n_residues), file=out)
  return model_list

def sort_models_and_sequences(models,sequences):
  # sort based on chain type if available
  sort_dict={}
  for m,s in zip(models,sequences):
    if not m or not m.hierarchy: continue
    chain_type=str(m.info.get('chain_type'))
    if not chain_type in sort_dict: sort_dict[chain_type]=[]
    sort_dict[chain_type].append([m,s])
  keys=list(sort_dict.keys())
  if len(keys)<2:
    return models,sequences # do nothing if only one chain type

  keys.sort()
  keys.reverse()
  new_models=[]
  new_sequences=[]
  for key in keys:
    for [m,s] in sort_dict[key]:
      new_models.append(m)
      new_sequences.append(s)
  return new_models,new_sequences

def merge_hierarchies_from_models(models=None,resid_offset=None,
    renumber=None,first_residue_number=None,
    sequences=None,chain_id=None,trim_side_chains=None,
    remove_ter_records=False,):
  # assumes one chain from each model
  # if resid_offset, space by to next even n of this number of residues
  # otherwise if renumber, start at first_residue_number and sequence
  #  consecutively
  # If sequence or chain_id are supplied, use them
  # Trim off side chains (and CB for GLY) if trim_side_chains
  # sort by chain_type if provided in one or more

  new_hierarchy=iotbx.pdb.input(
         source_info="Model",
             lines=flex.split_lines("")).construct_hierarchy()
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
  models,sequences=sort_models_and_sequences(models,sequences)
  for model,sequence in zip(models,sequences):
    if not model or not model.hierarchy: continue
    if not model.info:
        model.info={}
    if not info: info=model.info
    i=0
    for m in model.hierarchy.models():
      for chain in m.chains():
        if not renumber:
          cc=iotbx.pdb.hierarchy.chain()
          if chain_id:
            cc.id=chain_id
          elif chain.id: # take what is already there
            cc.id=chain.id
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
          if nn-resid<2: nn+=resid_offset
          resid=nn
  new_hierarchy.reset_atom_i_seqs()
  if trim_side_chains:
    atom_selection=\
      "name ca or name c or name o or name n or (name cb and not resname gly)"
    new_hierarchy=apply_atom_selection(atom_selection,hierarchy=new_hierarchy)

  if remove_ter_records:
    new_records=flex.split_lines("")
    for line in new_hierarchy.as_pdb_string().splitlines():
      if not line.startswith("TER"):
        new_records.append(line)

    new_hierarchy=iotbx.pdb.input(
         source_info="Model", lines=new_records).construct_hierarchy()

  new_model=model_info(hierarchy=new_hierarchy,info=info)

  return new_model

def get_average_direction(diffs=None, i=None,j=None):
    if not diffs: return None
    if i is None and j is None:
      i=0
      j=len(diffs)-1
    average_direction=col(diffs[i])
    nn=1.
    for j in range(i+1,j):
      nn+=1.
      average_direction+=col(diffs[j])
    average_direction/=nn
    return average_direction.normalize()

def get_chain_ids(hierarchy,unique_only=None):
  chain_ids=[]
  if not hierarchy:
    return chain_ids
  for model in hierarchy.models():
    for chain in model.chains():
      if (not unique_only) or (not chain.id in chain_ids):
        chain_ids.append(chain.id)
  return chain_ids

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

def get_atom_from_residue(residue=None,
  atom_name=None,allow_ca_only_model=False,
  skip_n_for_pro=False):
  if not residue:
    return None,None
  # just make sure that atom_name is there
  for atom in residue.atoms():
    if atom.name.replace(" ","")==atom_name.replace(" ",""):
      if skip_n_for_pro and atom.name.replace(" ","")=='N' and \
         residue.resname.replace(" ","").upper()=="PRO":
        return None,None
      else:
        return atom.name,atom.xyz
  if allow_ca_only_model:
    return atom_name,None
  else:
    return None,None

def get_indexed_residue(hierarchy,index=0):
  if not hierarchy:
    return None
  count=0
  for model in hierarchy.models():
    for chain in model.chains():
      for conformer in chain.conformers():
        for residue in conformer.residues():
          if count==index:
            return residue
          count+=1

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

def has_atom(hierarchy,name=None):
  if not hierarchy:
    return None
  for model in hierarchy.models():
    for chain in model.chains():
      for conformer in chain.conformers():
        for residue in conformer.residues():
          for atom in residue.atoms():
            if atom.name.replace(" ","")==name.replace(" ",""):
              return True
  return False

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

def get_all_resno(hierarchy):
  resno_list=[]
  if not hierarchy:
    return resno_list
  for model in hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        resno_list.append(rg.resseq_as_int())
  return resno_list

def get_middle_resno(hierarchy,first_resno=None,last_resno=None):
  if not hierarchy:
    return None
  if first_resno is None:
     first_resno=get_first_resno(hierarchy)
  if last_resno is None:
     last_resno=get_last_resno(hierarchy)
  target_resno=int(0.5+(first_resno+last_resno)/2)
  for model in hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        middle_resno=rg.resseq_as_int()
        if middle_resno >= target_resno:
          return middle_resno
  return last_resno



def verify_existence(hierarchy=None,
   prev_hierarchy=None,strand=None,registration=None,helix=None):
  ok=True

  if not hierarchy:
    raise Sorry("Need a model for user annotation input")

  for sh in [strand,helix]:
    if sh:
      start_atom=get_atom_index(hierarchy=hierarchy,
        atom_name=None,
        resname=sh.start_resname,
        chain_id=sh.start_chain_id,
        resseq=sh.start_resseq,
        icode=sh.start_icode)
      if start_atom is None:
        raise Sorry(
        "The starting residue in an annotation is missing in the input model:"+
        "%s %s %s %s" %(
         sh.start_resname,sh.start_chain_id,sh.start_resseq,sh.start_icode))

    if sh:
      end_atom=get_atom_index(hierarchy=hierarchy,
        atom_name=None,
        resname=sh.end_resname,
        chain_id=sh.end_chain_id,
        resseq=sh.end_resseq,
        icode=sh.end_icode)
      if end_atom is None:
        raise Sorry(
        "The ending residue in an annotation is missing in the input model:"+
        "%s %s %s %s" %(
         sh.end_resname,sh.end_chain_id,sh.end_resseq,sh.end_icode))

  if registration:
    cur_atom=get_atom_index(hierarchy=hierarchy,
      atom_name=registration.cur_atom,
      resname=registration.cur_resname,
      chain_id=registration.cur_chain_id,
      resseq=registration.cur_resseq,
      icode=registration.cur_icode)
    if cur_atom is None:
        raise Sorry(
        "A residue in an strand registration is missing in the input model:"+
        "%s %s %s %s %s" %(
        registration.cur_atom,registration.cur_resname,
        registration.cur_chain_id,registration.cur_resseq,
        registration.cur_icode))

  if registration and prev_hierarchy:
    prev_atom=get_atom_index(hierarchy=prev_hierarchy,
      atom_name=registration.prev_atom,
      resname=registration.prev_resname,
      chain_id=registration.prev_chain_id,
      resseq=registration.prev_resseq,
      icode=registration.prev_icode)
    if prev_atom is None:
      raise Sorry(
        "A residue in an strand registration is missing in the input model:"+
        "%s %s %s %s %s" %(
        registration.prev_atom,registration.prev_resname,
        registration.prev_chain_id,registration.prev_resseq,
        registration.prev_icode))

def get_atom_index(hierarchy=None,atom_name=None,resname=None,
 chain_id=None,resseq=None,icode=None):
  # NOTE: could speed up using atom cache and iselection
  # find the index in hierarchy of this atom
  if not hierarchy:
    return None
  count=-1
  for model in hierarchy.models():
    for chain in model.chains():
      for conformer in chain.conformers():
        for residue in conformer.residues():
          count+=1
          if not chain.id==chain_id: continue
          if not residue.resseq.replace(" ","")==str(resseq).replace(" ","") \
             or not residue.icode==icode:
            continue
          for atom in residue.atoms():
            if atom_name is None or atom.name==atom_name:
              return count
  return None

def have_n_or_o(models):
    have_n=False
    have_o=False
    for model in models:
      if model.has_n():
        have_n=True
        break
    for model in models:
      if model.has_o():
        have_o=True
        break
    if have_n and have_o:
      return True
    return False

def remove_bad_annotation(annotation,hierarchy=None,
     maximize_h_bonds=True,
     max_h_bond_length=None,
     maximum_length_difference=None,
     minimum_overlap=None,
     remove_overlaps=True,
     out=sys.stdout):
  # XXX perhaps move this to secondary_structure.py

  # remove parts of annotation that do not exist in the hierarchy or
  #  that overlap with other annotations
  new_annotation=annotation.split_sheets()

  deleted_something=False
  new_helices=[]
  for helix in new_annotation.helices:
    ph=apply_atom_selection(
      get_string_or_first_element_of_list(helix.as_atom_selections()),
         hierarchy=hierarchy)
    try:
        verify_existence(hierarchy=ph,helix=helix)
        new_helices.append(helix)
    except Exception as e: # the helix needs to be deleted
        print("%s\nDeleting this helix" %(str(e)), file=out)
        deleted_something=True
  if deleted_something:
    new_annotation.helices=new_helices

  new_sheets=[]
  for sheet in new_annotation.sheets:
    prev_hierarchy=None
    strands_ok=True
    registrations_ok=True
    for strand,registration in zip(sheet.strands,sheet.registrations):
      # verify that first and last atom selections in strand exist
      ph=apply_atom_selection(
         get_string_or_first_element_of_list(strand.as_atom_selections()),
         hierarchy=hierarchy)
      try:
        verify_existence(hierarchy=ph,prev_hierarchy=prev_hierarchy,
         strand=strand)
      except Exception as e: # the strand needs to be deleted
        print("%s\nDeleting this strand" %(str(e)), file=out)
        strands_ok=False
        deleted_something=True
      try:
        verify_existence(hierarchy=ph,prev_hierarchy=prev_hierarchy,
         registration=registration)
      except Exception as e: # the registration needs to be deleted
        print("Deleting registration\n%s" %(str(e)), file=out)
        registrations_ok=False
        deleted_something=True
      prev_hierarchy=ph
    if strands_ok:
      if not registrations_ok:
        sheet.registrations=[None,None]
      new_sheets.append(sheet)
  new_annotation.sheets=new_sheets

  if deleted_something:
    new_annotation.merge_sheets()
    if new_annotation.as_pdb_str():
      print("User annotation after removing bad parts:", file=out)
      print(new_annotation.as_pdb_str(), file=out)
    else:
      print("No annotation left after removing bad parts", file=out)
      return None

  if not remove_overlaps:
    return new_annotation

  # Now remove overlaps

  no_overlap_annotation=new_annotation.remove_overlapping_annotations(
    hierarchy=hierarchy,
    maximize_h_bonds=maximize_h_bonds,
    max_h_bond_length=max_h_bond_length,
    maximum_length_difference=maximum_length_difference,
    minimum_overlap=minimum_overlap)

  if not no_overlap_annotation.is_same_as(other=new_annotation):
    print("Edited annotation without overlaps:", file=out)
    print(no_overlap_annotation.as_pdb_str(), file=out)
    print(file=out)

  return no_overlap_annotation

def hierarchy_from_chain(chain):
  # create a hierarchy from a chain
  new_hierarchy=iotbx.pdb.input(
     source_info="Model", lines=flex.split_lines("")).construct_hierarchy()
  mm=iotbx.pdb.hierarchy.model()
  cc=chain.detached_copy()
  cc.id=chain.id  # copy chain ID
  new_hierarchy.append_model(mm)
  mm.append_chain(cc)
  return new_hierarchy


def get_string_or_first_element_of_list(something):
  if type(something)==type([1,2,3]):
    return something[0]
  else:
    return something

def sequence_from_hierarchy(hierarchy,chain_type='PROTEIN'):
  from iotbx.pdb import amino_acid_codes
  ott = amino_acid_codes.one_letter_given_three_letter

  sequence=""
  for model in hierarchy.models()[:1]:
    for chain in model.chains()[:1]:
      for rr in chain.residue_groups():
        for atom_group in rr.atom_groups()[:1]:
           sequence+=ott.get(atom_group.resname,"")
  return sequence

def sites_are_similar(sites1,sites2,max_rmsd=1):
  # return True if sites can be superimposed
  if max_rmsd is None:
    return True
  if max_rmsd < 0:
    return False
  if sites1.size() != sites2.size():
    return False
  if sites1.size() < 3:
    return False
  lsq_fit_obj=superpose.least_squares_fit(
       reference_sites=sites1,
       other_sites=sites2)
  sites2_fitted=lsq_fit_obj.other_sites_best_fit()
  rmsd = sites1.rms_difference(sites2_fitted)
  if rmsd <=max_rmsd:
    return True
  else:
    return False


def is_ca_only_hierarchy(hierarchy):
  if not hierarchy: return None
  asc=hierarchy.atom_selection_cache()
  atom_selection="not (name ca)"
  sel = asc.selection(string = atom_selection)
  if sel.count(True)==0:
    return True
  else:
    return False

def ca_n_and_o_always_present(hierarchy):
  if not hierarchy: return None
  total_residues=hierarchy.overall_counts().n_residues
  asc=hierarchy.atom_selection_cache()
  for n in ("ca","n","o"):
    atom_selection="name %s" %(n)
    sel = asc.selection(string = atom_selection)
    if sel.count(True) != total_residues:
      return False
  return True


def get_fraction_complete_backbone(hierarchy):
  if not hierarchy: return 0
  total_residues=hierarchy.overall_counts().n_residues
  asc=hierarchy.atom_selection_cache()
  minimum_complete=total_residues
  for n in ("ca","n","o"):
    atom_selection="name %s" %(n)
    sel = asc.selection(string = atom_selection)
    if sel.count(True) < minimum_complete:
      minimum_complete=sel.count(True)
  if minimum_complete==total_residues:
    return 1 # just to make sure it is exactly 1
  else:
    return minimum_complete/max(1,total_residues)

def choose_ca_or_complete_backbone(hierarchy, params=None):
  # Purpose: return a hierarchy that is either completely CA or has no CA-only
  #   residues
  fraction_complete_backbone=get_fraction_complete_backbone(hierarchy)
  if fraction_complete_backbone == 1:
    return hierarchy # nothing to do
  if fraction_complete_backbone < \
       params.find_ss_structure.min_ca_n_o_completeness:
    # just use CA-only
    return apply_atom_selection('name CA',hierarchy=hierarchy)
  else:  # remove CA-only residues
    hierarchy.remove_incomplete_main_chain_protein()
    return hierarchy

def sites_and_seq_from_hierarchy(hierarchy):
  atom_selection="name ca"
  sele=apply_atom_selection(atom_selection,hierarchy=hierarchy)
  if sele.overall_counts().n_residues==0:
    sites=flex.vec3_double()
    sequence=""
  else:
    sites=sele.extract_xray_structure().sites_cart()
    sequence=sequence_from_hierarchy(sele)
  start_resno=get_first_resno(sele)
  end_resno=get_last_resno(sele)
  return sites,sequence,start_resno,end_resno

class model_info: # mostly just a holder
  def __init__(self,hierarchy=None,id=0,info={},
      find_alpha=None,find_beta=None,find_other=None,
      find_three_ten=None,find_pi=None):
    self.hierarchy=hierarchy
    self.info=info
    self.id=id
    self.find_alpha=find_alpha
    self.find_three_ten=find_three_ten
    self.find_pi=find_pi
    self.find_beta=find_beta
    self.find_other=find_other

  def show_summary(self,out=sys.stdout):
    keys=list(self.info.keys())
    keys.sort()
    print("\nModel %d" %(self.id), file=out)
    for key in keys:
      print("%s: %s" %(key,str(self.info[key])), file=out)

  def has_n(self):
    return has_atom(self.hierarchy,name="N")

  def has_o(self):
    return has_atom(self.hierarchy,name="O")

  def first_residue(self):
    return get_first_resno(self.hierarchy)

  def last_residue(self):
    return get_last_resno(self.hierarchy)

  def length(self):
    return self.last_residue()-self.first_residue()+1

  def set_chain_type(self): # XXX specific for atom names N O2'/O2*
    # set the chain type if not already set
    if self.info and self.info.get('chain_type'): return
    if not self.info: self.info={}
    if has_atom(self.hierarchy,name="N") or has_atom(self.hierarchy,name="CA"):
      self.info['chain_type']="PROTEIN"
    else:
      if has_atom(self.hierarchy,name="O2'") or  \
          has_atom(self.hierarchy,name="O2*"):
        self.info['chain_type']="RNA"
      else:
        self.info['chain_type']="DNA"

class segment:  # object for holding a helix or a strand or other

  def setup(self,sites=None,start_resno=None,hierarchy=None,
     segment_class=None,
     segment_type='None',name=None,
     span=None,target_rise=None,residues_per_turn=None,
     rise_tolerance=None,dot_min=None,dot_min_single=None,
     target_i_ip3=None,tol_i_ip3=None,
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

    self.segment_class=segment_class
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
    self.target_i_ip3=target_i_ip3
    self.tol_i_ip3=tol_i_ip3
    self.dot_min=dot_min
    self.dot_min_single=dot_min_single
    self.diffs_single=None
    self.name=name # alpha helix, etc
    self.span=span # 3.5 for helices, 2 for strands, 0 for other
    self.target_rise=target_rise # 1.54 for helices, 3.3 strands
    self.residues_per_turn=residues_per_turn #3.6 alpha, 3 3-10, 4 pi, strand na
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
      print(text, file=self.out)

  def get_diffs_and_norms_3(self):
    # report distances between i and i+3
    sites_offset_3=self.sites[3:]
    diffs=sites_offset_3-self.sites[:-3]
    return diffs.norms()

  def get_min_diff_i_i3(self):
    norms=self.get_diffs_and_norms_3()
    if norms is None: return
    mm=norms.min_max_mean()
    return mm.min

  def trim_ends(self,start_pos=None,end_pos=None):
    # trim back from start_pos to end_pos (not residue number, position in
    #   existing segment
    start_res=self.get_start_resno()+start_pos
    end_res=self.get_start_resno()+end_pos
    if start_pos==0 and end_pos==self.length()-1:
      return # nothing to do
    atom_selection="resid %s through %s" %(resseq_encode(start_res),
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
    if sele.overall_counts().n_residues==0:
      self.sites=flex.vec3_double()
    else:
      # extract coordinates
      self.sites=sele.extract_xray_structure().sites_cart()

    if sele.overall_counts().n_residues==0:
      self.sites=flex.vec3_double()
    else:
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

  def segment_average_direction(self):
    diffs,norms=self.get_diffs_and_norms()
    if not diffs: return 0.
    return get_average_direction(diffs=diffs)

  def get_cosine(self):
    self.mean_dot_single=None
    # report mean normalized dot product of diffs with mean direction
    diffs,norms=self.get_diffs_and_norms()
    if not diffs: return 0.
    average_direction=get_average_direction(diffs=diffs)
    cosines=flex.double()
    for j in range(len(diffs)):
      dot=col(diffs[j]).dot(average_direction)
      cosines.append(dot)
    mean_dot=cosines.min_max_mean().mean

    # now get mean cross product and dot with average_direction...
    diffs_single=self.get_diffs_single()
    cosines=flex.double()
    for j in range(len(diffs_single)):
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
    target_i_ip3=self.target_i_ip3
    tol_i_ip3=self.tol_i_ip3
    tol=tol*self.span # rise_tolerance * span
    minimum_length=self.minimum_length
    return target,tol,dot_min,minimum_length,target_i_ip3,tol_i_ip3


  def is_ok(self,check_sub_segments=False,sub_segment_length=8):

    ##########################################
    # Option to check sub-segments instead of the whole thing (not used)
    if check_sub_segments and len(self.sites) > sub_segment_length:
      dd=max(1,sub_segment_length//2)
      start=0
      end=len(self.sites)-sub_segment_length
      for i in range(start,end+dd,dd):
        ii=min(i,end)
        h=self.segment_class(params=self.params,
          sites=self.sites[ii:ii+sub_segment_length])
        if not h.is_ok(check_sub_segments=False):
          return False
      return True
    ##########################################

    # check rise and cosine and dot product of single pairs with average
    #  direction and compare to target
    if (self.minimum_length is not None and self.length()<self.minimum_length) \
      or self.length() <1:
      return False
    rise=self.get_rise()
    dot=self.get_cosine()
    dot_single=self.mean_dot_single # set by get_cosine()
    min_diff_i_i3=self.get_min_diff_i_i3()
    if self.dot_min is not None and dot < self.dot_min:
      return False
    elif self.target_rise is not None and self.rise_tolerance is not None and \
       (rise < self.target_rise-self.rise_tolerance or \
           rise > self.target_rise+self.rise_tolerance):
      return False
    elif self.dot_min_single is not None \
       and dot_single < self.dot_min_single:
      return False
    elif self.target_i_ip3 is not None and self.tol_i_ip3 is not None and \
       min_diff_i_i3 < self.target_i_ip3-self.tol_i_ip3:
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

    if not self_orientation_points or not other_orientation_points:
      return None # none found

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
    atom_selection="resid %s through %s" %(resseq_encode(start_res),
       resseq_encode(end_res))
    sele=apply_atom_selection(atom_selection,hierarchy=hierarchy)

    asc=hierarchy.atom_selection_cache()
    sel = asc.selection(string = atom_selection)
    assert sele.overall_counts().n_residues  # should not be None
    f=StringIO()
    for atom in sele.atoms_with_labels():
      new_xyz=lsq_fit_obj.r * matrix.col(atom.xyz) + lsq_fit_obj.t
      atom.set_xyz(new_xyz)
      print(atom.format_atom_record(), file=f)
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
    self.params=params
    self.setup(sites=sites,start_resno=start_resno,hierarchy=hierarchy,
     segment_type='helix',
     segment_class=helix,
     minimum_length=params.minimum_length,
     buffer_residues=params.buffer_residues,
     standard_length=params.standard_length,
     frequency=frequency,
     base_score=params.base_score,
     is_n_terminus=is_n_terminus,
     is_c_terminus=is_c_terminus,
     name=params.name,
     span=params.span,
     target_rise=params.rise,
     residues_per_turn=params.residues_per_turn,
     rise_tolerance=params.rise_tolerance,
     target_i_ip3=params.target_i_ip3,
     tol_i_ip3=params.tol_i_ip3,
     dot_min=params.dot_min,
     dot_min_single=params.dot_min_single,
     optimal_delta_length=optimal_delta_length,
     verbose=verbose,
     out=out)

  def get_diffs_and_norms(self):
    # report distances between i and average of i+3,i+4
    if not hasattr(self,'diffs'):
      sites_offset_3=self.sites[3:-1]
      sites_offset_4=self.sites[4:]
      if self.residues_per_turn<=3:
        average_offset=sites_offset_3
      elif self.residues_per_turn>=4:
        average_offset=sites_offset_4
      else: # alpha helix, take average
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
    self.params=params

    self.setup(sites=sites,start_resno=start_resno,hierarchy=hierarchy,
      segment_type='strand',
      segment_class=strand,
      name=params.name,
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
      target_i_ip3=params.target_i_ip3,
      tol_i_ip3=params.tol_i_ip3,
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

    self.params=params
    self.setup(sites=sites,start_resno=start_resno,hierarchy=hierarchy,
      segment_type='other',
      segment_class=other,
      name=params.name,
      span=params.span,
      minimum_length=params.minimum_length,
      buffer_residues=params.buffer_residues,
      standard_length=params.standard_length,
      base_score=params.base_score,
      is_n_terminus=is_n_terminus,
      is_c_terminus=is_c_terminus,
      target_rise=params.rise,
      rise_tolerance=params.rise_tolerance,
      target_i_ip3=params.target_i_ip3,
      tol_i_ip3=params.tol_i_ip3,
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
      cut_up_segments=None,
      extend_segments=None,
      model_as_segment=None, # take the whole model as a segment
      verbose=None,
      out=sys.stdout):
    # Assumes model is just 1 chain of sequential residues
    #   obtained with split_model

    self.out=out
    self.extract_segments_from_pdb=extract_segments_from_pdb
    self.make_unique=make_unique
    self.cut_up_segments=cut_up_segments
    self.extend_segments=extend_segments
    self.model_as_segment=model_as_segment
    self.verbose=verbose
    self.params=params
    self.model=model
    self.name=params.name
    self.n_link_min=params.n_link_min

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
    diffs,norms,segment_dict=self.find_segments(params=params,sites=sites)

    # NOTE: separate find_segments method for each kind of segment
    # figure out how many residues really should be in these segments
    optimal_delta_length_dict,norm_dict=self.get_optimal_lengths(
       segment_dict=segment_dict,norms=norms)

    # create segment object (helix,strand) for each segment
    # If longer than the standard length, make a series of shorter segments

    keys=list(segment_dict.keys())
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
      # this segment starts at i and ends at segment_dict[i]
      overall_start_res=i
      overall_end_res=segment_dict[i]+self.last_residue_offset
      overall_length=overall_end_res-overall_start_res+1
      optimal_delta_length=optimal_delta_length_dict[overall_start_res]
      if (optimal_delta_length > 0 and not self.allow_insertions) or \
         (optimal_delta_length < 0 and not self.allow_deletions):
        optimal_delta_length=0
      # add to self.segments
      # cut into pieces of size self.standard_length
      if (not self.cut_up_segments) or overall_length<=self.standard_length:
        start_end_list.append([overall_start_res,overall_end_res])
      else:
        for start_res_use in range(
           overall_start_res,
           overall_start_res+overall_length-self.standard_length+1):
          start_end_list.append(
             [start_res_use,start_res_use+self.standard_length-1])
    self.segment_start_end_dict={}
    for start_res_use,end_res_use in start_end_list:
        self.segment_start_end_dict[start_res_use]=end_res_use
        ok=self.extract_segment(params=params,
         start_res=start_res_use,end_res=end_res_use,sites=sites,
         optimal_delta_length=optimal_delta_length)
        if not ok:
          del self.segment_start_end_dict[start_res_use]

  def show_summary(self,out=None):
    if not out: out=self.out
    for h in self.segments:
      print("Class: %12s  N: %d Start: %d End: %d " %(
          self.name,h.get_end_resno()-h.get_start_resno()+1,
           h.get_start_resno(),h.get_end_resno(),
         ) +" Rise:%5.2f A Dot:%5.2f" %(
          h.get_rise(),h.get_cosine()), file=out)

  def get_used_residues_list(self,end_buffer=1):
    # just return a list of used residues
    used_residues=[]
    if hasattr(self,'segment_dict'):
      for i in self.segment_start_end_dict.keys():
        for j in range(
             i+end_buffer,self.segment_start_end_dict[i]+1-end_buffer):
          # Note this segment_start_end_dict lists exact residues used,
          #  not starting points as in segment_dict
          if not j in used_residues:
            used_residues.append(j)
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
      if self.extract_segments_from_pdb or self.model_as_segment:
        start_resno=start_res_with_buffer+self.start_resno
        end_resno=end_res_with_buffer+self.start_resno
        # XXX NOTE: we assume that the model has been broken up so that there
        # is at most one residue with each residue number (broken up at
        # insertion codes).  If this changes this will not work properly.
        atom_selection="resseq %s:%s" %(
           resseq_encode(start_resno).replace(" ",""),
           resseq_encode(end_resno).replace(" ",""))
        hierarchy=apply_atom_selection(
           atom_selection,hierarchy=self.model.hierarchy)
        if hierarchy.overall_counts().n_residues==0:
          return False # did not find anything here and needed to
      else:
        hierarchy=None

      h=self.segment_class(params=params,
        sites=sites[start_res_with_buffer:end_res_with_buffer+1],
        start_resno=start_res_with_buffer+self.start_resno,
        optimal_delta_length=optimal_delta_length,hierarchy=hierarchy,
        is_n_terminus=is_n_terminus,
        is_c_terminus=is_c_terminus,
          )
      if h.is_ok() or self.model_as_segment:
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

  def find_segments(self,params=None,sites=None): # helix/strand
    # set up segment class for this kind of segment
    # for example, helix(sites=sites)
    h=self.segment_class(params=params,sites=sites)

    # get difference vectors i to i+2 (strands) or i to avg of i+3/i+4 (helix)
    diffs,norms=h.get_diffs_and_norms()  # difference vectors and lengths
    norms_3=h.get_diffs_and_norms_3()

    target,tol,dot_min,minimum_length,target_i_ip3,tol_i_ip3=h.get_targets()

    # now find sequence of N residues where diffs are parallel and values of
    # diffs are all about the same and about span*rise = 5.4 A for helix

    n=len(diffs)
    segment_dict={}
    used_residues=[]

    if len(sites)>0 and self.model_as_segment:  # take the whole thing
      segment_dict[0]=n-1
      return diffs,norms,segment_dict

    for i in range(n):
      if i in used_residues or abs(norms[i]-target)>tol: continue
      if i in self.previously_used_residues: continue
      # i is start of segment
      segment_dict[i]=i  # lists end of segment
      used_residues.append(i)
      for j in range(i+1,n):
        if abs(norms[j]-target)>tol: break
        if target_i_ip3 is not None and tol_i_ip3 is not None and \
           j<len(norms_3) and abs(norms_3[j]-target_i_ip3)>tol_i_ip3: break
        if j in self.previously_used_residues: break
        dot=col(diffs[j]).dot(col(diffs[j-1]))
        if dot < dot_min: break
        segment_dict[i]=j
        used_residues.append(j)

    # skip if only as long as a single span (requires 2 to be sure)
    for i in list(segment_dict.keys()):
      segment_length=segment_dict[i]+1+self.last_residue_offset-i
      if segment_length<minimum_length:
        del segment_dict[i]

    # prune out anything not really in this segment type based on direction

    segment_dict=self.remove_bad_residues(segment_dict=segment_dict,
      diffs=diffs,dot_min=dot_min,minimum_length=minimum_length)

    # merge any segments that can be combined
    segment_dict=self.merge_segments(params=params,
       segment_dict=segment_dict,sites=sites)

    # try to extend by 1 in each direction if it doesn't overlap anything
    if self.extend_segments:
      segment_dict=self.try_to_extend_segments(params=params,
        segment_dict=segment_dict,sites=sites)

    # trim ends of any segments that are connected by fewer than n_link_min=3
    #  residues and that go in opposite directions (connected by tight turns).

    segment_dict=self.trim_short_linkages(
     params=params,segment_dict=segment_dict,sites=sites)

    # again merge any segments that can be combined
    segment_dict=self.merge_segments(params=params,
       segment_dict=segment_dict,sites=sites)

    return diffs,norms,segment_dict


  def try_to_extend_segments(self,params=None,segment_dict=None,sites=None):
    # try to extend segments if they do not overlap
    found=True
    n_cycles=0
    while found and n_cycles <=len(segment_dict):
      found=False
      n_cycles+=1
      keys=list(segment_dict.keys())
      keys.sort()  # sorted on start
      n_space=2  # at least 2 apart in the end
      for i1,i2,i3 in zip([None]+keys[:-2],keys,keys[1:]+[None]):
        if found: break
        if i2 > 0 and ( i1 is None or
            i2 > segment_dict[i1]+self.last_residue_offset+1+n_space):
          # try adding before i2
          h=self.segment_class(params=params,
             sites=sites[i2-1:segment_dict[i2]+self.last_residue_offset+1],
             start_resno=i2-1+self.start_resno)
          if h.is_ok():
           found=True
           segment_dict[i2-1]=segment_dict[i2]
           del segment_dict[i2]
        if found: break
        if segment_dict[i2]+self.last_residue_offset+1+1 < sites.size()  and (
           i3 is None or
            segment_dict[i2]+self.last_residue_offset+1+n_space  < i3):
          # try adding after i2
          h=self.segment_class(params=params,
             sites=sites[i2:segment_dict[i2]+self.last_residue_offset+1+1],
             start_resno=i2+self.start_resno)
          if h.is_ok():
           found=True
           segment_dict[i2]=segment_dict[i2]+1

    return segment_dict


  def merge_segments(self,params=None,segment_dict=None,sites=None):

    # merge any segments that can be combined
    found=True
    n_cycles=0
    while found and n_cycles <=len(segment_dict):
      found=False
      n_cycles+=1
      keys=list(segment_dict.keys())
      keys.sort()  # sorted on start
      for i1,i2 in zip(keys[:-1],keys[1:]):
        if found: break
        if i2 <= segment_dict[i1]+self.last_residue_offset+1:
          # could merge. Test it
          h=self.segment_class(params=params,
             sites=sites[i1:segment_dict[i2]+self.last_residue_offset+1],
             start_resno=i1+self.start_resno)
          if h.is_ok():
           found=True
           segment_dict[i1]=segment_dict[i2]
           if i1 != i2:
             del segment_dict[i2]
    return segment_dict

  def remove_bad_residues(self,segment_dict=None,diffs=None,dot_min=None,
    minimum_length=None):
    still_changing=True
    n_cycles=0
    max_cycles=2
    while still_changing and n_cycles < max_cycles:
      still_changing=False
      new_segment_dict={}
      keys=list(segment_dict.keys())
      keys.sort()
      for i in keys:
        average_direction=get_average_direction(
          diffs=diffs,i=i,j=segment_dict[i])
        segment_start=None
        for j in range(i,segment_dict[i]+1):
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
    return segment_dict

  def trim_short_linkages(self,params=None,segment_dict=None,sites=None):
    found=True
    n_cycles=0
    while found and n_cycles <=len(segment_dict):
      found=False
      n_cycles+=1
      keys=list(segment_dict.keys())
      keys.sort()  # sorted on start
      for i1,i2 in zip(keys[:-1],keys[1:]):
        if found: break
        end_1=segment_dict[i1]+self.last_residue_offset
        start_2=i2
        delta=i2-(segment_dict[i1]+self.last_residue_offset)
        if delta >=self.n_link_min: continue
        h1=self.segment_class(params=params,
          sites=sites[i1:segment_dict[i1]+self.last_residue_offset+1],
             start_resno=i1+self.start_resno)
        h2=self.segment_class(params=params,
          sites=sites[i2:segment_dict[i2]+self.last_residue_offset+1],
             start_resno=i2+self.start_resno)
        if h1.segment_average_direction().dot(h2.segment_average_direction())>0:
          continue # (in generally the same direction...ignore)
        # in opposite direction
        residues_to_cut=(delta+1)//2
        found=True
        segment_dict[i1]=segment_dict[i1]-residues_to_cut
        ok=False
        if segment_dict[i1]>=i1: # check and keep
          h1=self.segment_class(params=params,
            sites=sites[i1:segment_dict[i1]+self.last_residue_offset+1],
             start_resno=i1+self.start_resno)
          ok=h1.is_ok()
        if not ok:
            del segment_dict[i1]

        ok=False
        new_i2=i2+residues_to_cut
        if segment_dict[i2]>=new_i2: # check and keep
          h2=self.segment_class(params=params,
            sites=sites[new_i2:segment_dict[i2]+self.last_residue_offset+1],
               start_resno=new_i2+self.start_resno)
          ok=h2.is_ok()
        if ok:
          segment_dict[new_i2]=segment_dict[i2]
        # remove original
        if new_i2 != i2:
          del segment_dict[i2]

    return segment_dict


  def get_optimal_lengths(self,segment_dict=None,norms=None):

    keys=list(segment_dict.keys())
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
      if mean_norm is None: # give up
        return optimal_delta_length_dict,norm_dict
      # guess number of residues it should be: standard rise
      segment_rise=mean_norm/self.span  # we used mean of i+3 and i+4 for helix
      segment_length=segment_dict[i]+1+self.last_residue_offset-i
      segment_distance=segment_rise*(segment_length-1)
      optimal_length=int(0.5+segment_distance/standard_rise)+1
      optimal_delta_length_dict[i]=optimal_length-segment_length
      norm_dict[i]=mean_norm
      if self.verbose:
        print("Chain start: %d  end: %d  N: %d  Rise: %7.2f  Optimal length: %d" %(
         i+self.start_resno,i+segment_length-1+self.start_resno,
         segment_length,segment_rise,optimal_length), file=self.out)

    return optimal_delta_length_dict,norm_dict

class find_helix(find_segment):

  def __init__(self,params=None,model=None,verbose=None,
    extract_segments_from_pdb=None,
    make_unique=None,
    previously_used_residues=None,
    model_as_segment=None,
    cut_up_segments=None,
    extend_segments=None,
    out=sys.stdout):
    self.previously_used_residues=previously_used_residues

    self.setup(params=params,model=model,segment_type='helix',
     extract_segments_from_pdb=extract_segments_from_pdb,
     make_unique=make_unique,
     model_as_segment=model_as_segment,
     cut_up_segments=cut_up_segments,
     extend_segments=extend_segments,
     verbose=verbose,out=out)

  def pdb_records(self,segment_list=None,last_id=0,helix_type='alpha',
     max_h_bond_length=None,
     force_secondary_structure_input=None,
     require_h_bonds=None,
     minimum_h_bonds=None,
     maximum_poor_h_bonds=None,
     allow_ca_only_model=None,out=sys.stdout): # helix

    records=[]
    number_of_good_h_bonds=0
    number_of_poor_h_bonds=0
    k=last_id
    for s in segment_list:
      if not s.hierarchy: continue
      f=StringIO()
      all_h_bonds,n_good,n_poor=self.list_h_bonds(
          segment=s,helix_type=helix_type,
          max_h_bond_length=max_h_bond_length,
          force_secondary_structure_input=force_secondary_structure_input,
          allow_ca_only_model=allow_ca_only_model,out=f)
      if require_h_bonds:
        if n_good<minimum_h_bonds or (
            maximum_poor_h_bonds is not None and n_poor>maximum_poor_h_bonds):
          continue

      print(f.getvalue(), end=' ', file=out)
      number_of_good_h_bonds+=n_good
      number_of_poor_h_bonds+=n_poor
      start=get_first_residue(s.hierarchy)
      end=get_last_residue(s.hierarchy)
      chain_id=get_chain_id(s.hierarchy)
      k=k+1
      record = secondary_structure.pdb_helix(
        serial=k,
        helix_id=k,
        start_resname=start.resname,
        start_chain_id=chain_id,
        start_resseq=start.resseq,
        start_icode=start.icode,
        end_resname=end.resname,
        end_chain_id=chain_id,
        end_resseq=end.resseq,
        end_icode=end.icode,
        helix_class=secondary_structure.pdb_helix.helix_class_to_int(
           helix_type), # 1=alpha 3=pi  5=3_10
        comment="",
        length=s.length())
      records.append(record)
    return records,number_of_good_h_bonds,number_of_poor_h_bonds

  def list_h_bonds(self,segment=None,helix_type='alpha',
     max_h_bond_length=None,
     force_secondary_structure_input=None,
     allow_ca_only_model=None,out=sys.stdout):

    helix_class=secondary_structure.pdb_helix.helix_class_to_int(
           helix_type) # 1=alpha 3=pi  5=3_10

    # residue that residue i is h-bonded to
    next_i_dict={
      1:4,   # alpha:  O of residue i H-bonds to N of residue i+4
      3:5,   # pi:     O of residue i H-bonds to N of residue i+5
      5:3,   # 3_10:   O of residue i H-bonds to N of residue i+3
     }
    next_i=next_i_dict[helix_class]

    all_h_bonds=[]
    number_of_good_h_bonds=0
    number_of_poor_h_bonds=0

    for i in range(segment.length()-next_i):

      cur_residue=get_indexed_residue(
        segment.hierarchy,index=i)
      cur_atom,cur_xyz=get_atom_from_residue(
        residue=cur_residue,
        atom_name=' O  ',allow_ca_only_model=allow_ca_only_model)

      next_residue=get_indexed_residue(
        segment.hierarchy,index=i+next_i)
      next_atom,next_xyz=get_atom_from_residue(
        residue=next_residue,
        atom_name=' N  ',allow_ca_only_model=allow_ca_only_model,
        skip_n_for_pro=True)

      if cur_xyz and next_xyz:  # both present at least
        dd=col(cur_xyz)-col(next_xyz)
        dist=dd.length()
        if dist <=max_h_bond_length:
          bad_one=""
          number_of_good_h_bonds+=1
        else:
          bad_one="**"
          number_of_poor_h_bonds+=1

        new_h_bond=h_bond(
           prev_atom=cur_atom,
           prev_resname=cur_residue.resname,
           prev_chain_id=get_chain_id(segment.hierarchy),
           prev_resseq=cur_residue.resseq,
           prev_icode=cur_residue.icode,
           cur_atom=next_atom,
           cur_resname=next_residue.resname,
           cur_chain_id=get_chain_id(segment.hierarchy),
           cur_resseq=next_residue.resseq,
           cur_icode=next_residue.icode,
           dist=dist,
           bad_one=bad_one,
           anything_is_ok=force_secondary_structure_input,
         )
        new_h_bond.show_summary(out=out,show_non_existent=False)
        all_h_bonds.append(new_h_bond)
      else:
        bad_one=None
        dist=None
    return all_h_bonds,number_of_good_h_bonds,number_of_poor_h_bonds

class find_beta_strand(find_segment):

  def __init__(self,params=None,model=None,verbose=None,
      extract_segments_from_pdb=None,
      make_unique=None,
      cut_up_segments=None,
      extend_segments=None,
      model_as_segment=None,
      previously_used_residues=None,
      out=sys.stdout):

    self.previously_used_residues=previously_used_residues

    self.setup(params=params,
      model=model,segment_type='strand',
      extract_segments_from_pdb=extract_segments_from_pdb,
      make_unique=make_unique,
      model_as_segment=model_as_segment,
      cut_up_segments=cut_up_segments,
      extend_segments=extend_segments,
      verbose=verbose,out=out)

  def get_pdb_strand(self,sheet_id=None,strand_id=1,segment=None,
     sense=0,start_index=None,end_index=None):

    if start_index is None:
      start=get_first_residue(segment.hierarchy)
    else:
      start=get_indexed_residue(segment.hierarchy,index=start_index)
    if end_index is None:
      end=get_last_residue(segment.hierarchy)
    else:
      end=get_indexed_residue(segment.hierarchy,index=end_index)
    if start is None or end is None:
      return None

    chain_id=get_chain_id(segment.hierarchy)
    pdb_strand = secondary_structure.pdb_strand(
        sheet_id=sheet_id,
        strand_id=strand_id,
        start_resname=start.resname,
        start_chain_id=chain_id,
        start_resseq=start.resseq,
        start_icode=start.icode,
        end_resname=end.resname,
        end_chain_id=chain_id,
        end_resseq=end.resseq,
        end_icode=end.icode,
        sense=sense)
    return pdb_strand

  def get_required_start_end(self,sheet=None,info_dict=None):
      start_dict={}
      end_dict={}
      for i in sheet:
        start_dict[i]=None
        end_dict[i]=None
      for i,j in zip(sheet[:-1],sheet[1:]):
        key="%d:%d" %(i,j)
        first_last_1_and_2=info_dict[key]
        first_ca_1,last_ca_1,first_ca_2,last_ca_2,is_parallel,i_index,j_index=\
           first_last_1_and_2

        if start_dict[i] is None or start_dict[i]>first_ca_1:
           start_dict[i]=first_ca_1
        if end_dict[i] is None or end_dict[i]<last_ca_1:
           end_dict[i]=last_ca_1

        if start_dict[j] is None or start_dict[j]>first_ca_2:
           start_dict[j]=first_ca_2
        if end_dict[j] is None or end_dict[j]<last_ca_2:
           end_dict[j]=last_ca_2

      return start_dict,end_dict

  def pdb_records(self,segment_list=None,   # sheet
     sheet_list=None,info_dict=None,allow_ca_only_model=None,
     force_secondary_structure_input=None,
     max_h_bond_length=None,
     require_h_bonds=None,
     minimum_h_bonds=None,
     maximum_poor_h_bonds=None,
     out=sys.stdout):

    # sheet_list is list of sheets. Each sheet is a list of strands (the index
    #  of the strand in segment_list). Info_dict has the relationship between
    #  pairs of strands, indexed with the key "%d:%d:" %(i,j) where i and j are
    #  the indices of the two strands. The dictionary returns
    #  [first_ca_1,last_ca_1,first_ca_2,last_ca_2,is_parallel,i_index,j_index]
    #    for the two strands

    records=[]
    number_of_good_h_bonds=0
    number_of_poor_h_bonds=0
    sheet_id=0
    for sheet in sheet_list:
      good_in_sheet=0
      poor_in_sheet=0
      f=StringIO()

      sheet_id+=1
      strand_id=1
      k=sheet[0]
      remainder=sheet[1:] # all the others
      s=segment_list[k]
      if not s.hierarchy or not s.hierarchy.overall_counts().n_residues:
        continue
      current_sheet = secondary_structure.pdb_sheet(
        sheet_id=sheet_id,
        n_strands=len(sheet),
        strands=[],
        registrations=[])
      ok=True

      # figure out what residues to include in each sheet. It is not a well-
      #   defined problem because a middle strand might H-bond to one but not
      #   both of its neighbors even if both neighbors are beta-strands
      start_dict,end_dict=self.get_required_start_end(sheet=sheet,
          info_dict=info_dict)

      first_strand = self.get_pdb_strand(sheet_id=sheet_id,strand_id=strand_id,
        segment=s,sense=0,start_index=start_dict[k],end_index=end_dict[k])
      if first_strand is None:
        print("Note: failed to identify strand %s in sheet %s index %s" %(
          k,strand_id,sheet_id), file=out)
        ok=False
        continue # found nothing (something went wrong in get_pdb_strand)

      current_sheet.add_strand(first_strand)
      current_sheet.add_registration(None)
      previous_s=s
      i=k
      for j in remainder: # previous strand is i, current is j
        if not ok: break
        s=segment_list[j]
        strand_id+=1
        if not s.hierarchy or not s.hierarchy.overall_counts().n_residues:
          continue

        key="%d:%d" %(i,j)
        first_last_1_and_2=info_dict[key]
        first_ca_1,last_ca_1,first_ca_2,last_ca_2,is_parallel,i_index,j_index=\
           first_last_1_and_2
        # sense is whether previous and current strands are parallel (1) or
        #   antiparallel (-1)

        if is_parallel:
          sense=1
        else:
          sense=-1
        next_strand = self.get_pdb_strand(sheet_id=sheet_id,strand_id=strand_id,
          segment=s,sense=sense,start_index=start_dict[j],end_index=end_dict[j])
        if next_strand is None:
          print("Note: failed to "+\
            "identify strand %s in sheet %s index %s" %(j,strand_id,sheet_id), file=out)
          ok=False
          continue# found nothing (something went wrong in get_pdb_strand)

        current_sheet.add_strand(next_strand)
        all_h_bonds,n_good,n_poor=self.list_h_bonds(
          segment=s,
          max_h_bond_length=max_h_bond_length,
          previous_segment=previous_s,first_last_1_and_2=first_last_1_and_2,
          force_secondary_structure_input=force_secondary_structure_input,
          allow_ca_only_model=allow_ca_only_model,out=f)
        good_in_sheet+=n_good
        poor_in_sheet+=n_poor

        register=self.get_pdb_strand_register(segment=s,
          previous_segment=previous_s,first_last_1_and_2=first_last_1_and_2,
          allow_ca_only_model=allow_ca_only_model,
          all_h_bonds=all_h_bonds)
        current_sheet.add_registration(register)
        previous_s=s
        i=j

      if not ok: continue # skip

      if require_h_bonds:
        if good_in_sheet<minimum_h_bonds or (
            maximum_poor_h_bonds is not None and n_poor>maximum_poor_h_bonds):
          sheet_id-=1
          continue
      print(f.getvalue(), end=' ', file=out)
      number_of_good_h_bonds+=good_in_sheet
      number_of_poor_h_bonds+=poor_in_sheet
      records.append(current_sheet)


    return records,number_of_good_h_bonds,number_of_poor_h_bonds


  def is_even(self,i):
    if 2*(i//2)==i: return True
    return False

  def get_pdb_strand_register(self,segment=None,previous_segment=None,
     first_last_1_and_2=None,allow_ca_only_model=None,
     all_h_bonds=None):

    for h_bond in all_h_bonds: # choose first that is ok
      if not h_bond.is_ok():
        continue

      from iotbx.pdb.secondary_structure import pdb_strand_register
      register=pdb_strand_register(
             cur_atom=h_bond.cur_atom,
             cur_resname=h_bond.cur_resname,
             cur_chain_id=h_bond.cur_chain_id,
             cur_resseq=h_bond.cur_resseq,
             cur_icode=h_bond.cur_icode,
             prev_atom=h_bond.prev_atom,
             prev_resname=h_bond.prev_resname,
             prev_chain_id=h_bond.prev_chain_id,
             prev_resseq=h_bond.prev_resseq,
             prev_icode=h_bond.prev_icode,)
      return register
      # just take the first good one

  def list_h_bonds(self,segment=None,previous_segment=None,
     max_h_bond_length=None,
     force_secondary_structure_input=None,
     first_last_1_and_2=None,allow_ca_only_model=None,out=sys.stdout):

    #  Looking down a strand in direction from N to C...
    #    the CA go up-down-up-down.
    #    The ones that are up have their O pointing to the right
    #    Those that are down have O pointing to the left
    #  So...if we orient strand n+1 from N to C...if
    #    strand n is to the right then choose an "up" residue of strand n+1 for
    #    the matching to strand n.  If strand n is to the left choose a "down"
    #    one.

    #  If CA residue i of strand n matches with residue i' of strand n+1:

    #  For antiparallel strands:
    #  O of residue i in strand n H-bonds to N of residue i' in strand n+1

    #  For parallel strands:
    #  O of residue i in strand n H-bonds to N of residue i'+1 in
    #   strand n+1.

    # Here strand n is previous_segment and n+1 is segment

    # Get entire list of H-bonded residues between these segments.
    # Residues in previous_segment go from first_ca_1 to last_ca_1.
    # We have already specified that i_index of previous_segment
    #   atom O H-bonds to j_index of segment atom N.
    # For antiparallel strands, other H-bond is N of i_index with O of j_index
    # For parallel strands, other H-bond is N of i_index with O of j_index-2

    # Increase i_index by 2 and decrease j_index
    #  by 2 and the same pattern occurs.

    # look at entire segments, not just the part we are including

    all_h_bonds=[]
    number_of_poor_h_bonds=0
    number_of_good_h_bonds=0

    first_ca_1,last_ca_1,first_ca_2,last_ca_2,is_parallel,i_index,j_index=\
           first_last_1_and_2

    for i in range(previous_segment.length()):
      if i_index is None or j_index is None: continue
      if not self.is_even(i-i_index): continue

      local_i_index=i_index+(i-i_index)
      if is_parallel:
        local_j_index=j_index+(i-i_index)
      else:
        local_j_index=j_index-(i-i_index)

      # make sure we are in range
      if local_i_index<0 or local_i_index>previous_segment.length()-1: continue

      local_prev_residue=get_indexed_residue(
        previous_segment.hierarchy,index=local_i_index)
      for o_to_n in [True,False]:
        if o_to_n:
          if local_j_index<0 or local_j_index>segment.length()-1: continue

          local_cur_residue=get_indexed_residue(
            segment.hierarchy,index=local_j_index)
          local_prev_atom,local_prev_xyz=get_atom_from_residue(
            residue=local_prev_residue,
            atom_name=' O  ',allow_ca_only_model=allow_ca_only_model)
          local_cur_atom,local_cur_xyz=get_atom_from_residue(
            residue=local_cur_residue,
            atom_name=' N  ',allow_ca_only_model=allow_ca_only_model,
            skip_n_for_pro=True)
        else:
          if is_parallel:
            if local_j_index-2<0 or local_j_index-2>segment.length()-1: continue
            local_cur_residue=get_indexed_residue(
              segment.hierarchy,index=local_j_index-2)
          else:
            if local_j_index<0 or local_j_index> segment.length()-1: continue
            local_cur_residue=get_indexed_residue(
              segment.hierarchy,index=local_j_index)
          local_prev_atom,local_prev_xyz=get_atom_from_residue(
            residue=local_prev_residue,
            atom_name=' N  ',allow_ca_only_model=allow_ca_only_model,
            skip_n_for_pro=True)
          local_cur_atom,local_cur_xyz=get_atom_from_residue(
            residue=local_cur_residue,
            atom_name=' O  ',allow_ca_only_model=allow_ca_only_model)

        # skip if there is no atom pair (e.g., if N and PRO )
        if not local_prev_atom or not local_cur_atom: continue

        if local_cur_xyz and local_prev_xyz:
          dd=col(local_cur_xyz)-col(local_prev_xyz)
          dist=dd.length()
          if dist <=max_h_bond_length:
            bad_one=""
          else:
            bad_one="**"
        else:
          bad_one=None
          dist=None

        # Save those that are outside range we are keeping only if good:
        if local_j_index<first_ca_2 or local_j_index>last_ca_2 or \
           local_i_index<first_ca_1 or local_i_index>last_ca_1:
          if bad_one!="": continue
          # and mark as not included
          bad_one="(Not included) "

        if bad_one=="":
          number_of_good_h_bonds+=1
        elif bad_one=="**":
          number_of_poor_h_bonds+=1

        new_h_bond=h_bond(
             prev_atom=local_prev_atom,
             prev_resname=local_prev_residue.resname,
             prev_chain_id=get_chain_id(previous_segment.hierarchy),
             prev_resseq=local_prev_residue.resseq,
             prev_icode=local_prev_residue.icode,
             cur_atom=local_cur_atom,
             cur_resname=local_cur_residue.resname,
             cur_chain_id=get_chain_id(segment.hierarchy),
             cur_resseq=local_cur_residue.resseq,
             cur_icode=local_cur_residue.icode,
             dist=dist,
             bad_one=bad_one,
             anything_is_ok=force_secondary_structure_input,
         )
        new_h_bond.show_summary(out=out,show_non_existent=False)
        all_h_bonds.append(new_h_bond)
    return all_h_bonds,number_of_good_h_bonds,number_of_poor_h_bonds

class h_bond:  # holder for a pair of atoms
  def __init__(self,
             prev_atom=None,
             prev_resname=None,
             prev_chain_id=None,
             prev_resseq=None,
             prev_icode=None,
             cur_atom=None,
             cur_resname=None,
             cur_chain_id=None,
             cur_resseq=None,
             cur_icode=None,
             dist=None,
             bad_one=None,
             anything_is_ok=False):
    adopt_init_args(self, locals())

  def is_ok(self):
    if self.anything_is_ok:
      return True
    if self.dist is None:  # was CA-only so no information
       return True
    elif self.dist and not self.bad_one: # was ok H-bond
       return True
    else:
      return False

  def show_summary(self,show_non_existent=False,out=sys.stdout):
    if self.dist is not None:
      print(" %4s%4s%4s%5s%s : %4s%4s%4s%5s%s :: %5.2f   %s" %(
             self.prev_atom,
             self.prev_resname,
             self.prev_chain_id,
             self.prev_resseq,
             self.prev_icode,
             self.cur_atom,
             self.cur_resname,
             self.cur_chain_id,
             self.cur_resseq,
             self.cur_icode,
             self.dist,
             self.bad_one), file=out)
    elif self.bad_one is not None:
      print(" %4s%4s%4s%5s%s : %4s%4s%4s%5s%s" %(
             self.prev_atom,
             self.prev_resname,
             self.prev_chain_id,
             self.prev_resseq,
             self.prev_icode,
             self.cur_atom,
             self.cur_resname,
             self.cur_chain_id,
             self.cur_resseq,
             self.cur_icode,), file=out)


class find_other_structure(find_segment):

  def __init__(self,previously_used_residues=None,
      params=None,model=None,
      extract_segments_from_pdb=None,
      make_unique=None,
      cut_up_segments=None,
      extend_segments=None,
      verbose=None,out=sys.stdout):

    self.previously_used_residues=previously_used_residues

    self.setup(params=params,model=model,segment_type='other',
      extract_segments_from_pdb=extract_segments_from_pdb,
      make_unique=make_unique,
      cut_up_segments=cut_up_segments,
      extend_segments=extend_segments,
      verbose=verbose,out=out)

  def pdb_records(self,last_id=0,out=sys.stdout):   #other (nothing)
    return []

  def get_optimal_lengths(self,segment_dict=None,norms=None):
    # always zero
    optimal_delta_length_dict={}
    norm_dict={}
    for key in segment_dict:
      optimal_delta_length_dict[key]=0
      norm_dict[key]=None
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

  def find_segments(self,params=None,sites=None): # other
    # find everything that is not alpha and not beta, put buffer_residues
    #  buffer on the end of each one.

    # set up segment class for this kind of segment SPECIFIC FOR OTHER
    # for example, helix(sites=sites)
    h=self.segment_class(params=params,sites=sites)

    n=len(sites)
    if n>0 and self.model_as_segment:  # take the whole thing
      segment_dict[0]=n-1-self.last_residue_offset
      return None,None,segment_dict

    # cross off used residues except for buffer of buffer_residues
    n_buf=params.buffer_residues
    used_residues=n*[False]
    used_residues=[]
    for i in range(n+1):
      if i in self.previously_used_residues:
        used_residues.append(True)
      else:
        used_residues.append(False)

    segment_dict={}
    segment_start=None
    segment_end=None
    for i,used in zip(range(n),used_residues):
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

class helix_strand_segments:
  def __init__(self):
    self.h_bond_text=""
    self.all_strands=[]
    self.all_alpha_helices=[]
    self.all_three_ten_helices=[]
    self.all_pi_helices=[]
    self.sheet_list=[]
    self.used_strands=[]
    self.pair_dict={}
    self.info_dict={}

    self.pdb_alpha_helix_list=[]
    self.pdb_three_ten_helix_list=[]
    self.pdb_pi_helix_list=[]
    self.pdb_sheet_list=[]

  def add_from_model(self,model):
      if model.find_beta:
        self.all_strands+=model.find_beta.segments
      if model.find_alpha:
        self.all_alpha_helices+=model.find_alpha.segments
      if model.find_three_ten:
        self.all_three_ten_helices+=model.find_three_ten.segments
      if model.find_pi:
        self.all_pi_helices+=model.find_pi.segments

  def find_sheets(self,out=sys.stdout,
     max_sheet_ca_ca_dist=6.,
     min_sheet_length=4,
     include_single_strands=None):
    if not self.all_strands: return
    print("\nFinding sheets from %d strands" %(len(
        self.all_strands)), file=out)
    self.get_strand_pairs(tol=max_sheet_ca_ca_dist,
        min_sheet_length=min_sheet_length)
    # pair_dict is list of all the strands that each strand matches with
    # self.info_dict is information on a particular pair
    #   of strands:
    #  self.info_dict["%d:%d" %(i,j)]=
    #     [first_ca_1,last_ca_1,first_ca_2,last_ca_2,is_parallel]

    # Create sheets from paired strands.
    # keep track of which ones we have assigned
    self.used_strands=[]

    # single_strands are those with no matching strands
    single_strands=self.get_strands_by_pairs(pairs=0)
    pair_strands=self.get_strands_by_pairs(pairs=1)
    triple_strands=self.get_strands_by_pairs(pairs=2)
    multiple_strands=self.get_strands_by_pairs(pairs=None)

    self.used_strands=[] # initialize again
    self.used_strands+=single_strands
    # we are going to ignore these

    self.sheet_list=[]

    if include_single_strands:# include singles
      for i in single_strands:
        self.sheet_list.append([i])

    # find all sheets with edges (beginning with a paired strand)
    self.sheet_list+=self.get_sheets_from_edges(
      pair_strands=pair_strands)
    self.sheet_list+=self.get_sheets_from_edges(
      pair_strands=triple_strands)

    # any pairs remaining? Create specialized "sheet" for each one
    existing_pairs_in_sheets=self.get_existing_pairs_in_sheets()
    missing_pairs=[]
    for i in range(len(self.all_strands)):
      for j in self.pair_dict.get(i,[]):
        if not [i,j] in existing_pairs_in_sheets and \
           not [i,j] in missing_pairs and not [j,i] in missing_pairs:
          missing_pairs.append([i,j])
          self.sheet_list.append([i,j])
          if not i in self.used_strands:
            self.used_strands.append(i)
          if not j in self.used_strands:
            self.used_strands.append(j)

    # Now we are ready to create sheets from sheet_list, self.pair_dict and
    #   self.info_dict

  def get_existing_pairs_in_sheets(self=None):
    existing_pairs=[]
    for sheet in self.sheet_list:
      for i,j in zip(sheet[:-1],sheet[1:]):
        existing_pairs.append([i,j])
        existing_pairs.append([j,i])
    return existing_pairs

  def get_sheets_from_edges(self,pair_strands=None):
    sheet_list=[]
    for i in pair_strands:
      if i in self.used_strands:continue
      strand_list=[i]
      current_strand=i
      while current_strand is not None:
        current_strand=self.get_available_strand(
          current_strand=current_strand,
          strand_list=strand_list)
        if current_strand is not None:
          strand_list.append(current_strand)
      if len(strand_list)>1: # require an actual sheet
        self.used_strands+=strand_list
        sheet_list.append(strand_list)
    return sheet_list


  def get_available_strand(self,current_strand=None,strand_list=None):
    for i in self.pair_dict.get(current_strand,[]):
      if not i in self.used_strands and not i in strand_list:
         return i
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
    for i in range(n):
      if i in used_strands: continue
      if pairs is None or len(self.pair_dict.get(i,[]))==pairs:
        return i
    return None

  def get_strand_pairs(self,tol=None,min_sheet_length=None):
    self.info_dict={}
    self.pair_dict={}
    for i in range(len(self.all_strands)):
      self.pair_dict[i]=[]
    for i in range(len(self.all_strands)):
      for j in range(i+1,len(self.all_strands)):
        self.ca1=None
        self.ca2=None
        if self.ca_pair_is_close(self.all_strands[i],self.all_strands[j],
            tol=tol):

          # figure out alignment and whether it really is ok
          first_last_1_and_2=self.align_strands(
            self.all_strands[i],self.all_strands[j],tol=tol,
            min_sheet_length=min_sheet_length)

          if first_last_1_and_2:
            # we have a match
            [first_ca_1,last_ca_1,first_ca_2,last_ca_2,is_parallel]=\
               first_last_1_and_2

            # figure out which "O" of strand i should H-bond with which "N" of j
            # it is either going to be first_ca_1 or first_ca_1+1

            i_index,j_index=self.get_ind_h_bond_sheet(
              first_last_1_and_2=first_last_1_and_2,i=i,j=j,switch_i_j=False)

            self.pair_dict[i].append(j)
            self.info_dict["%d:%d" %(i,j)]=\
               [first_ca_1,last_ca_1,first_ca_2,last_ca_2,is_parallel,
                i_index,j_index]

            # and make an entry for the other way around
            i_index,j_index=self.get_ind_h_bond_sheet(
              first_last_1_and_2=first_last_1_and_2,i=i,j=j,switch_i_j=True)
            self.pair_dict[j].append(i)
            self.info_dict["%d:%d" %(j,i)]=\
               [first_ca_2,last_ca_2,first_ca_1,last_ca_1,is_parallel,
                i_index,j_index]

  def get_ind_h_bond_sheet(self,
      first_last_1_and_2,i=None,j=None,switch_i_j=False,
      registration=None,force_secondary_structure_input=None,
      sense=None):

    if switch_i_j:
      xx=i
      i=j
      j=xx
      [first_ca_2,last_ca_2,first_ca_1,last_ca_1,is_parallel]=first_last_1_and_2
    else:
      [first_ca_1,last_ca_1,first_ca_2,last_ca_2,is_parallel]=first_last_1_and_2

    strand_i=self.all_strands[i]
    strand_j=self.all_strands[j]

    #-----------------force_secondary_structure_input------------
    if force_secondary_structure_input:
      if registration and strand_i.hierarchy and strand_j.hierarchy:
        # use the registration to id the paired atoms if provided
        if switch_i_j:
          prev_atom=registration.cur_atom
          i_bond=get_atom_index(hierarchy=strand_i.hierarchy,
           atom_name=registration.cur_atom,
           resname=registration.cur_resname,
           chain_id=registration.cur_chain_id,
           resseq=registration.cur_resseq,
           icode=registration.cur_icode)
          j_bond=get_atom_index(hierarchy=strand_j.hierarchy,
           atom_name=registration.prev_atom,
           resname=registration.prev_resname,
           chain_id=registration.prev_chain_id,
           resseq=registration.prev_resseq,
           icode=registration.prev_icode)
        else:
          prev_atom=registration.prev_atom
          strand_i=self.all_strands[i]
          strand_j=self.all_strands[j]
          i_bond=get_atom_index(hierarchy=strand_i.hierarchy,
           atom_name=registration.prev_atom,
           resname=registration.prev_resname,
           chain_id=registration.prev_chain_id,
           resseq=registration.prev_resseq,
           icode=registration.prev_icode)
          j_bond=get_atom_index(hierarchy=strand_j.hierarchy,
           atom_name=registration.cur_atom,
           resname=registration.cur_resname,
           chain_id=registration.cur_chain_id,
           resseq=registration.cur_resseq,
           icode=registration.cur_icode)

        assert sense in [-1,1]
        i_index=i_bond
        if sense==1: #  parallel strands:
          #  O of residue i in strand n H-bonds to N of i'+1 in strand n+1.
          #  N of residue i in strand n H-bonds to O of i'-1 in strand n+1.
          # O of i_index (prev) H-bonds to N of j_index (cur)
          # N of i_index H-bonds to O of j_index-2

          assert prev_atom.replace(" ","") in ["O","N"]
          if prev_atom.replace(" ","")=="N":
            # H-bond is to j_bond which is j_index-2
            j_index=j_bond+2
            if j_index>strand_j.sites.size()-1:
              j_index-=2
              i_index-=2
              if i_index>strand_i.sites.size()-1 or \
                  j_index>strand_j.sites.size()-1:
                return None,None # failed
          else:
            # H-bond is to j_bond which is j_index
            j_index=j_bond

        else: #  antiparallel strands:
          #  O of residue i in strand n H-bonds to N of residue i' in strand n+1
          #  N of residue i in strand n H-bonds to O of residue i' in strand n+1
          j_index=j_bond# everything is ok already


        return i_index,j_index

      else:
        return None,None
    #-----------------end force_secondary_structure_input------------



    i_index=first_ca_1
    if is_parallel:
      j_index=first_ca_2+1
    else:
      j_index=last_ca_2

    # View strand i from N to C with strand j to the right.  Every other
    #  residue in strand i has CA up/down/up/down.  Choose either
    #  residue first_ca_1 or first_ca_1+1, whichever is more up in this
    #  reference frame.
    #  "Up" here is
    # CA(i_index) -> CA(j_index)  X strand_i.segment_average_direction()

    inter_strand_vector=col(strand_j.get_sites()[j_index])-\
                        col(strand_i.get_sites()[i_index])
    if inter_strand_vector.is_zero():
      return None,None # give up (could not find a suitable H-bond)

    inter_strand_vector=inter_strand_vector.normalize()
    up_direction=inter_strand_vector.cross(
      strand_i.segment_average_direction())
    if up_direction.is_zero():
      return None,None # give up (could not find a suitable H-bond)

    up_direction=up_direction.normalize()

    n_dot=0.
    sum_dot=0.
    last_offset_index=len(strand_i.get_sites())-i_index-2
    for i in range(last_offset_index//2+1):
      offset=2*i
      delta=col(strand_i.get_sites()[i_index+1+offset])- \
            col(strand_i.get_sites()[i_index+offset])
      dot=up_direction.dot(delta)
      n_dot+=1
      sum_dot+=dot
    if n_dot:
      dot=sum_dot/n_dot
    else:
      return None,None # give up (could not find a suitable direction

    if dot > 0: # i_index is down. (dot is positive). Move 1 residue ahead

      i_index=i_index+1
      if is_parallel:
        j_index=j_index+1
      else:
        j_index=j_index-1

    if i_index+1>len(strand_i.get_sites()) or \
        j_index+1>len(strand_j.get_sites()) or \
        j_index< 0:
      return None,None # give up (could not find a suitable H-bond)
    return i_index,j_index


  def align_strands(self,s1,s2,tol=None,
     min_sheet_length=None):
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


    for offset in [0]:  # using other offsets did not help and sometimes worse
      if self.ca1+offset < 0 or self.ca1+offset > len(sites1)-1: continue
      dd_list,keep_1,keep_2=self.get_residue_pairs_in_sheet(sites1,sites2,
       center1=self.ca1+offset,center2=self.ca2,tol=tol)
      dd_list_reverse,keep_1_reverse,keep_2_reverse=\
         self.get_residue_pairs_in_sheet(sites1,sites2_reversed,
         center1=self.ca1+offset,center2=len(sites2)-self.ca2-1,tol=tol)

      if len(keep_1)<min_sheet_length and len(keep_1_reverse)<min_sheet_length:
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
      for offset in range(start_offset,end_offset+1):
        i1=center1+offset
        i2=center2+offset
        dd=(col(sites1[i1])-col(sites2[i2])).norm_sq()
        if dd <= tol**2:
          keep1_list.append(i1)
          keep2_list.append(i2)
          dd_list.append(dd**0.5)
        elif offset>0: # passed the middle, so end it
          break
        else: # have a bad one and have not gotten to middle yet. Start over
          keep_1_list=[]
          keep_2_list=[]
          dd_list=[]

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
      for i in range(jump//2,len(sites1),jump):
        for j in range(jump//2,len(sites2),jump):
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

  def set_up_pdb_records(self,allow_ca_only_model=None,
    max_h_bond_length=None,
    force_secondary_structure_input=None,
    require_h_bonds=None,
    minimum_h_bonds=None,
    maximum_poor_h_bonds=None,
    models=None,out=sys.stdout):

    # skip any secondary structure elements that have fewer than minimum_h_bonds
    #   if require_h_bonds is True

    number_of_good_h_bonds=0
    number_of_poor_h_bonds=0
    # save everything as pdb_alpha_helix or pdb_sheet objects
    if not models:
      return number_of_good_h_bonds,number_of_poor_h_bonds

    # determine if there are N and O atoms present.
    # If not, set allow_ca_only_model=True
    allow_ca_only_model=(not have_n_or_o(models))

    f=StringIO()
    print("\nList of H-bonds expected from helices and from strand pairings", file=f)
    print("Distances > %3.1f A indicated by **" %(max_h_bond_length), file=f)
    print("H-bonds not included in HELIX/SHEET records marked 'Not included'", file=f)
    print("\n      ATOM 1               ATOM 2           Dist (A)\n", file=f)
    def get_fa(find_what='find_alpha',models=None):
       for model in models:
        if getattr(model,find_what):
          return getattr(model,find_what)
    fa=get_fa(find_what='find_alpha',models=models)
    if fa and self.all_alpha_helices:
      self.pdb_alpha_helix_list,n_good,n_poor=fa.pdb_records(
         segment_list=self.all_alpha_helices,
         helix_type='alpha',
         max_h_bond_length=max_h_bond_length,
         force_secondary_structure_input=force_secondary_structure_input,
         require_h_bonds=require_h_bonds,
         minimum_h_bonds=minimum_h_bonds,
         maximum_poor_h_bonds=maximum_poor_h_bonds,
         allow_ca_only_model=allow_ca_only_model,out=f)
      number_of_good_h_bonds+=n_good
      number_of_poor_h_bonds+=n_poor

    fa=get_fa(find_what='find_three_ten',models=models)
    if fa and self.all_three_ten_helices:
      self.pdb_three_ten_helix_list,n_good,n_poor=fa.pdb_records(
         segment_list=self.all_three_ten_helices,
         helix_type='3_10',
         max_h_bond_length=max_h_bond_length,
         require_h_bonds=require_h_bonds,
         minimum_h_bonds=minimum_h_bonds,
         maximum_poor_h_bonds=maximum_poor_h_bonds,
         force_secondary_structure_input=force_secondary_structure_input,
         allow_ca_only_model=allow_ca_only_model,out=f)
      number_of_good_h_bonds+=n_good
      number_of_poor_h_bonds+=n_poor

    fa=get_fa(find_what='find_pi',models=models)
    if fa and self.all_pi_helices:
      self.pdb_pi_helix_list,n_good,n_poor=fa.pdb_records(
        segment_list=self.all_pi_helices,helix_type='pi',
        max_h_bond_length=max_h_bond_length,
        require_h_bonds=require_h_bonds,
        minimum_h_bonds=minimum_h_bonds,
        maximum_poor_h_bonds=maximum_poor_h_bonds,
        force_secondary_structure_input=force_secondary_structure_input,
        allow_ca_only_model=allow_ca_only_model,out=f)
      number_of_good_h_bonds+=n_good
      number_of_poor_h_bonds+=n_poor

    fa=get_fa(find_what='find_beta',models=models)
    if fa and self.sheet_list and self.all_strands:
      self.pdb_sheet_list,n_good,n_poor=fa.pdb_records(
        segment_list=self.all_strands,
       sheet_list=self.sheet_list,
       info_dict=self.info_dict,
       max_h_bond_length=max_h_bond_length,
       require_h_bonds=require_h_bonds,
       minimum_h_bonds=minimum_h_bonds,
       maximum_poor_h_bonds=maximum_poor_h_bonds,
       force_secondary_structure_input=force_secondary_structure_input,
       allow_ca_only_model=allow_ca_only_model,out=f)
      number_of_good_h_bonds+=n_good
      number_of_poor_h_bonds+=n_poor
    text=f.getvalue()
    print(text, file=out)
    self.h_bond_text=text
    return number_of_good_h_bonds,number_of_poor_h_bonds

  def save_pdb_records_as_string(self):
    all_pdb_records=StringIO()
    self.all_selection_records=[]
    for helix in self.pdb_alpha_helix_list:
      print(helix.as_pdb_str(), file=all_pdb_records)
      self.all_selection_records+=helix.as_atom_selections()

    for helix in self.pdb_three_ten_helix_list:
      print(helix.as_pdb_str(), file=all_pdb_records)
      self.all_selection_records+=helix.as_atom_selections()

    for helix in self.pdb_pi_helix_list:
      print(helix.as_pdb_str(), file=all_pdb_records)
      self.all_selection_records+=helix.as_atom_selections()

    for sheet in self.pdb_sheet_list:
      print(sheet.as_pdb_str(), file=all_pdb_records)
      self.all_selection_records+=sheet.as_atom_selections()

    self.all_pdb_records=all_pdb_records.getvalue()

  def get_pdb_alpha_helix_list(self):
    if hasattr(self,'pdb_alpha_helix_list'):
      return self.pdb_alpha_helix_list

  def get_pdb_three_ten_helix_list(self):
    if hasattr(self,'pdb_three_ten_helix_list'):
      return self.pdb_three_ten_helix_list

  def get_pdb_pi_helix_list(self):
    if hasattr(self,'pdb_pi_helix_list'):
      return self.pdb_pi_helix_list

  def get_pdb_sheet_list(self):
    if hasattr(self,'pdb_sheet_list'):
      return self.pdb_sheet_list

  def get_all_pdb_records(self):
    if hasattr(self,'all_pdb_records'):
      return self.all_pdb_records

  def get_annotation(self):
    self.save_pdb_records_as_string()
    records=self.all_pdb_records
    import iotbx.pdb.secondary_structure as ioss
    return ioss.annotation.from_records(records=flex.split_lines(records))

  def get_all_selection_records(self):
    if not hasattr(self,'all_selection_records'):
       return
    text='"'
    first=True
    for sel in self.all_selection_records:
      if not first:
        text+=" or "
      first=False
      text+=" ( "+sel.replace('"','')+") "
    text+='"'
    return text

class fss_result_object:
  # A holder for results of find_secondary_structure
  def __init__(self,
      id=None,
      chain_id=None,
      hierarchy=None,
      number_of_good_h_bonds=None,
      number_of_poor_h_bonds=None,
      h_bond_text=None,
      annotation=None,
      sequence=None,
      sites=None,
      start_resno=None,
      end_resno=None,
      max_rmsd=1):
    adopt_init_args(self, locals())

    if hierarchy:  #
      self.sites,self.sequence,self.start_resno,self.end_resno=\
         sites_and_seq_from_hierarchy( hierarchy)

    self.chain_id_list=[]
    if self.chain_id:
      self.chain_id_list.append(self.chain_id)

  def show_summary(self,out=sys.stdout):
    print("\nSummary of find_secondary_structure object %s" %(self.id), end=' ', file=out)
    print("for sequence: %s_%s::%s\n" %(self.start_resno,self.end_resno,
       self.sequence), file=out)
    print("Chain ID's where this applies: %s" %(
      " ".join(self.get_chain_id_list())), file=out)
    print("Good H-bonds: %s  Poor H-bonds: %s " %(
      self.number_of_good_h_bonds,self.number_of_poor_h_bonds), file=out)
    print("H-bond text:\n%s" %(self.h_bond_text), file=out)
    print("Annotation:\n%s" %(self.get_annotation()), file=out)

  def add_chain_id(self,chain_id=None):
    if chain_id is not None:
      self.chain_id_list.append(chain_id)


  def add_info(self,chain_id=None,
      number_of_good_h_bonds=None,
      number_of_poor_h_bonds=None,
      h_bond_text=None,
      annotation=None,):
    if chain_id is not None:
      self.chain_id_list.append(chain_id)
    if number_of_good_h_bonds is not None:
      self.number_of_good_h_bonds=number_of_good_h_bonds
    if number_of_poor_h_bonds is not None:
      self.number_of_poor_h_bonds=number_of_poor_h_bonds
    if h_bond_text is not None:
      self.h_bond_text=h_bond_text
    if annotation is not None:
      self.annotation=annotation

  def get_annotation(self):
    return self.annotation

  def get_chain_id_list(self):
    return self.chain_id_list

  def is_similar_fss_result(self,other):
    if self.sequence != other.sequence:
       return False
    if not sites_are_similar(self.sites,other.sites,max_rmsd=self.max_rmsd):
       return False
    if self.start_resno != other.start_resno:
       return False
    if self.end_resno != other.end_resno:
       return False
    return True

class conformation_group:
  # A group of fss_results with the same sequence but different
  #  conformations
  def __init__(self,
    fss_result=None,
    ):

    self.fss_results=[]
    self.last_id=0
    if fss_result:
       self.last_id+=1
       fss_result.id=self.last_id
       self.fss_results.append(fss_result)

  def __repr__(self):
    text="Conformation group with %s fss_results" %(
       len(self.get_fss_result_list()))
    if self.get_fss_result_list():
      text+="\nSequence: %s" %( self.get_fss_result_list()[0].sequence)
    return text

  def get_fss_result_list(self):
    return self.fss_results

  def add_fss_result(self,fss_result=None):
    if fss_result:
       self.last_id+=1
       fss_result.id=self.last_id
       self.fss_results.append(fss_result)

  def get_similar_fss_result(self,fss_result=None):
    for fss_r in self.fss_results:
      if fss_r.is_similar_fss_result(fss_result):
        return fss_r
    return None


class find_secondary_structure: # class to look for secondary structure

  def __init__(self,params=None,args=None,hierarchy=None,models=None,
      user_annotation_text=None,max_h_bond_length=None,
      force_secondary_structure_input=None,
      search_secondary_structure=None,
      combine_annotations=None,
      require_h_bonds=None,
      minimum_h_bonds=None,
      maximum_poor_h_bonds=None,
      helices_are_alpha=False,
      ss_by_chain=None,
      use_representative_chains=None,
      max_representative_chains=None,
      max_rmsd=None,
      verbose=None,out=sys.stdout):
    adopt_init_args(self, locals())

    if not args: args=[]

    if not params:  # get params from args if necessary
      params=self.get_params(args,out=out)

    if verbose is not None:
      params.control.verbose=verbose
    verbose=params.control.verbose


    if helices_are_alpha or params.find_ss_structure.helices_are_alpha:
      params.find_ss_structure.find_three_ten=False
      params.find_ss_structure.find_pi=False
    if require_h_bonds is not None:
      params.find_ss_structure.require_h_bonds=require_h_bonds
    if minimum_h_bonds is not None:
      params.find_ss_structure.minimum_h_bonds=minimum_h_bonds
    if maximum_poor_h_bonds is not None:
      params.find_ss_structure.maximum_poor_h_bonds=maximum_poor_h_bonds
    if max_h_bond_length is not None:
      params.find_ss_structure.max_h_bond_length=max_h_bond_length
    if combine_annotations is not None:
      params.find_ss_structure.combine_annotations=combine_annotations
    if force_secondary_structure_input is not None:
      params.input_files.force_secondary_structure_input=\
         force_secondary_structure_input
    force_secondary_structure_input=\
      params.input_files.force_secondary_structure_input
    if ss_by_chain is not None:
      params.find_ss_structure.ss_by_chain=ss_by_chain
    if max_rmsd is not None:
      params.find_ss_structure.max_rmsd=max_rmsd
    if use_representative_chains is not None:
      params.find_ss_structure.use_representative_chains=\
        use_representative_chains
    if max_representative_chains is not None:
      params.find_ss_structure.max_representative_chains=\
        max_representative_chains

    secondary_structure_input=params.input_files.secondary_structure_input

    if search_secondary_structure is not None:
      params.find_ss_structure.search_secondary_structure=\
         search_secondary_structure
    search_secondary_structure=\
       params.find_ss_structure.search_secondary_structure

    if (not params.find_ss_structure.search_secondary_structure) and (not
      params.input_files.secondary_structure_input) and (not
      user_annotation_text):
      raise Sorry(
       "Need either secondary_structure_input or search_secondary_structure=True")

    self.helix_strand_segments=helix_strand_segments()
    self.user_helix_strand_segments=helix_strand_segments()

    self.user_models=[]
    self.number_of_good_h_bonds=0
    self.number_of_poor_h_bonds=0
    self.user_number_of_good_h_bonds=0
    self.user_number_of_poor_h_bonds=0
    self.h_bond_text=""

    if verbose:
      local_out=out
    else:
      from libtbx.utils import null_out
      local_out=null_out()

    if hierarchy:
      hierarchy=hierarchy.deep_copy()
    elif (not models) or params.input_files.secondary_structure_input:
      # need to get a hierarchy
      if models:
        combined_model=merge_hierarchies_from_models(models=models)
        hierarchy=combined_model.hierarchy.deep_copy()
      else:  # read it in
        if not params.input_files.pdb_in or \
            not os.path.isfile(params.input_files.pdb_in):
          raise Sorry("Missing file: %s" %(params.input_files.pdb_in))
        hierarchy=get_pdb_hierarchy(text=open(params.input_files.pdb_in).read())
    if hierarchy:
      hierarchy.remove_alt_confs(always_keep_one_conformer=False)
      atom_selection="protein"
      try:
        hierarchy=apply_atom_selection(atom_selection,hierarchy=hierarchy)
      except Exception as e:
        hierarchy=None

    if hierarchy and params.find_ss_structure.auto_choose_ca_vs_ca_n_o:
      hierarchy=choose_ca_or_complete_backbone(hierarchy,params=params)

    if force_secondary_structure_input and not \
        (params.input_files.secondary_structure_input or user_annotation_text):
      raise Sorry(
         "Need secondary_structure_input for force_secondary_structure_input")

    # Get user ss information if any into composite_user_annotation
    if user_annotation_text or params.input_files.secondary_structure_input:
      composite_user_annotation=self.get_user_ss(
        params=params,hierarchy=hierarchy,
        user_annotation_text=user_annotation_text,out=out)
      if not params.input_files.secondary_structure_input:
        params.input_files.secondary_structure_input=True # so we can check
        secondary_structure_input=True # so we can check
    else:
      composite_user_annotation=None

    if models:
      self.models=models
    else:
      self.models=split_model(hierarchy=hierarchy)

    # Decide if we are going to run in parts and just extend those to all
    #   copies
    self.args=args
    self.params=params
    self.hierarchy=hierarchy
    if self.need_to_run_in_parts():
       self.run_in_parts()
       return  # done

    if force_secondary_structure_input or (not
       params.find_ss_structure.search_secondary_structure):
      working_annotation=composite_user_annotation
    else:
      for model in self.models:
        self.find_ss_in_model(params=params,model=model,out=out)
        self.helix_strand_segments.add_from_model(model)

      if self.helix_strand_segments and \
         params.find_ss_structure.set_up_helices_sheets:
        self.helix_strand_segments.find_sheets(
         include_single_strands=params.find_ss_structure.include_single_strands,
         max_sheet_ca_ca_dist=params.beta.max_sheet_ca_ca_dist,
         min_sheet_length=params.beta.min_sheet_length,
         out=out) # organize strands into sheets

      if self.helix_strand_segments:
        self.number_of_good_h_bonds,self.number_of_poor_h_bonds=\
          self.helix_strand_segments.set_up_pdb_records(models=self.models,
          max_h_bond_length=params.find_ss_structure.max_h_bond_length,
          force_secondary_structure_input=force_secondary_structure_input,
          allow_ca_only_model=params.beta.allow_ca_only_model,out=local_out)
        self.h_bond_text=self.helix_strand_segments.h_bond_text
        print("\nNumber of good H-bonds: %d  Number of poor H-bonds: %d" %(
          self.number_of_good_h_bonds,self.number_of_poor_h_bonds), file=local_out)

      # get annotation:
      print("\nNew working annotation:", file=out)
      working_annotation=self.helix_strand_segments.get_annotation()
      print(working_annotation.as_pdb_str(), file=out)

      if params.find_ss_structure.combine_annotations and \
          composite_user_annotation:
        print("\nMerging edited input annotation and working annotation", file=out)
        working_annotation=composite_user_annotation.combine_annotations(
          hierarchy=hierarchy, other=working_annotation)
        print("\nMerged annotation:\n", file=out)
        print(working_annotation.as_pdb_str(), file=out)


    #  Remove annotation that does not match model
    if params.find_ss_structure.remove_missing_atom_annotation:
      working_annotation=remove_bad_annotation(
        working_annotation,
        hierarchy=hierarchy,
        max_h_bond_length=params.find_ss_structure.max_h_bond_length,
        remove_overlaps=False, # XXX Required to prevent recursion
        out=out)

    # Now get final values of H-bonds etc with our final annotation

    if params.find_ss_structure.require_h_bonds:
      remove_text=""
      if params.find_ss_structure.minimum_h_bonds>0:
        remove_text+=\
         "\nRemoving any secondary structure with fewer than %d H-bonds"  %(
        params.find_ss_structure.minimum_h_bonds)
      if params.find_ss_structure.maximum_poor_h_bonds and \
         params.find_ss_structure.maximum_poor_h_bonds>0:
        remove_text+=\
         "\nRemoving any secondary structure with more than %d poor H-bonds"  %(
        params.find_ss_structure.maximum_poor_h_bonds)

    else:
      remove_text=""

    if self.helix_strand_segments and not secondary_structure_input:
      # Use analysis from working annotation (no user input)
      if remove_text: print(remove_text, file=out)
      print("\nGetting H-bonds from working annotation", file=out)
      self.number_of_good_h_bonds,self.number_of_poor_h_bonds=\
         self.helix_strand_segments.set_up_pdb_records(models=self.models,
         max_h_bond_length=params.find_ss_structure.max_h_bond_length,
         force_secondary_structure_input=force_secondary_structure_input,
         allow_ca_only_model=params.beta.allow_ca_only_model,
         require_h_bonds=params.find_ss_structure.require_h_bonds,
         minimum_h_bonds=params.find_ss_structure.minimum_h_bonds,
         maximum_poor_h_bonds=params.find_ss_structure.maximum_poor_h_bonds,
         out=out)
      self.h_bond_text=self.helix_strand_segments.h_bond_text
      working_annotation=self.helix_strand_segments.get_annotation()

    elif self.user_helix_strand_segments and secondary_structure_input and \
        not params.find_ss_structure.combine_annotations:
      # use analysis of user input (as is)
      if remove_text and not force_secondary_structure_input:
        print(remove_text, file=out)
        require_h_bonds=params.find_ss_structure.require_h_bonds
        minimum_h_bonds=params.find_ss_structure.minimum_h_bonds
        maximum_poor_h_bonds=params.find_ss_structure.maximum_poor_h_bonds,
      else:
        remove_text=None
        require_h_bonds=None
        minimum_h_bonds=None
        maximum_poor_h_bonds=None

      print("\nGetting H-bonds from user annotation", file=out)
      self.number_of_good_h_bonds,self.number_of_poor_h_bonds=\
         self.user_helix_strand_segments.set_up_pdb_records(
         models=self.user_models,
         max_h_bond_length=params.find_ss_structure.max_h_bond_length,
         force_secondary_structure_input=force_secondary_structure_input,
         require_h_bonds=require_h_bonds,
         minimum_h_bonds=minimum_h_bonds,
         maximum_poor_h_bonds=params.find_ss_structure.maximum_poor_h_bonds,
         allow_ca_only_model=params.beta.allow_ca_only_model,
         out=out)
      self.h_bond_text=self.user_helix_strand_segments.h_bond_text
      working_annotation=self.user_helix_strand_segments.get_annotation()

    else: # need to redo it from the beginning with our new annotation
      # user annotation combined with new annotation
      if remove_text: print(remove_text, file=out)
      print("\nRunning analysis with new annotation", file=out)
      if working_annotation and working_annotation.as_pdb_str():
        fss=find_secondary_structure(hierarchy=hierarchy,
          user_annotation_text=working_annotation.as_pdb_str(),
          force_secondary_structure_input=True,
          combine_annotations=False,
          max_h_bond_length=params.find_ss_structure.max_h_bond_length,
          require_h_bonds=params.find_ss_structure.require_h_bonds,
          minimum_h_bonds=params.find_ss_structure.minimum_h_bonds,
          maximum_poor_h_bonds=params.find_ss_structure.maximum_poor_h_bonds,
          ss_by_chain=params.find_ss_structure.ss_by_chain,
          use_representative_chains=\
            params.find_ss_structure.use_representative_chains,
          max_representative_chains=\
            params.find_ss_structure.max_representative_chains,
          max_rmsd=params.find_ss_structure.max_rmsd,
          out=local_out)
        print(fss.h_bond_text, file=out)
        self.number_of_good_h_bonds=fss.number_of_good_h_bonds
        self.number_of_poor_h_bonds=fss.number_of_poor_h_bonds
        working_annotation=fss.get_annotation()
      else:
        self.number_of_good_h_bonds=0
        self.number_of_poor_h_bonds=0

    print("\nNumber of good H-bonds: %d  Number of poor H-bonds: %d" %(
          self.number_of_good_h_bonds,self.number_of_poor_h_bonds), file=out)

    self.annotation=working_annotation

    self.show_summary(verbose=params.control.verbose,
      pdb_records_file=params.output_files.pdb_records_file,out=out)

  def need_to_run_in_parts(self,
     min_residues_for_parts=None,
     min_average_chain_length=None,):
    # run in parts if lots of ncs or big chains.
    # Don't if lots of little fragments or model objects are supplied.
    if not self.params.find_ss_structure.ss_by_chain:
      return # not going to do this at all
    if not self.hierarchy:
      return # not going to do this at all. Only from hierarchy

    oc=self.hierarchy.overall_counts()
    if min_residues_for_parts and oc.n_residues < min_residues_for_parts:
      return
    if min_average_chain_length and \
       oc.n_residues/max(1,oc.n_chains) < min_average_chain_length:
      return

    # If not a CA-only model, require N and O to be present on all residues
    #   to run with representative chains (otherwise there may be some N/O
    #   that are present in only some chains)
    if (not is_ca_only_hierarchy(self.hierarchy))  and (
          not ca_n_and_o_always_present(self.hierarchy)):
      return

    # Worth running on individual chains
    return True

  def run_in_parts(self):
    print("\nRunning on full chains (no "+\
        "intra-chain secondary structure)", file=self.out)

    # Just run through all the chains and get their ss.  If duplicate chains
    #   and use_representative_chains, copy results
    local_params=deepcopy(self.params)
    local_params.find_ss_structure.ss_by_chain=False
    if self.params.control.verbose:
      local_out=self.out
    else:
      from libtbx.utils import null_out
      local_out=null_out()
    result_dict={} # fss_conformation_groups keyed by sequence to find quickly
    unique_sequence_list=[]

    for model in self.hierarchy.models()[:1]:
      for chain in model.chains():
        chain_id=chain.id
        local_hierarchy=hierarchy_from_chain(chain)
        # get fss_result holder
        current_fss_result=fss_result_object(chain_id=chain_id,
           hierarchy=local_hierarchy,
           max_rmsd=self.params.find_ss_structure.max_rmsd)
        if self.params.find_ss_structure.use_representative_chains and \
          len(unique_sequence_list)< \
            self.params.find_ss_structure.max_representative_chains:
          # See if this sequence has been analyzed already:
          test_sequence_string="%s_%s::%s" %(
            current_fss_result.start_resno,current_fss_result.end_resno,
            current_fss_result.sequence)
          cg=result_dict.get(test_sequence_string,conformation_group())
          # cg is either empty or a conformation_group with current sequence
          existing_fss_result=cg.get_similar_fss_result(current_fss_result)
          # if present, existing_fss_result is same conformation as current
        else:
          cg=conformation_group()
          existing_fss_result=None
        if existing_fss_result:
          existing_fss_result.add_chain_id(chain_id=chain_id)
        else: # get the analysis of this chain
          fss=find_secondary_structure(
            params=local_params,hierarchy=local_hierarchy,
            user_annotation_text=self.user_annotation_text,
            max_h_bond_length=self.max_h_bond_length,
            force_secondary_structure_input=\
              self.force_secondary_structure_input,
            search_secondary_structure=self.search_secondary_structure,
            combine_annotations=self.combine_annotations,
            require_h_bonds=self.require_h_bonds,
            minimum_h_bonds=self.minimum_h_bonds,
            maximum_poor_h_bonds=self.maximum_poor_h_bonds,
            verbose=self.verbose,out=local_out)
          current_fss_result.add_info(  # add new information
             number_of_good_h_bonds=fss.number_of_good_h_bonds,
             number_of_poor_h_bonds=fss.number_of_poor_h_bonds,
             h_bond_text=fss.h_bond_text,
             annotation=fss.get_annotation())
          # add new fss_result to empty conformation_group and save
          cg.add_fss_result(fss_result=current_fss_result)
          sequence_string="%s_%s::%s" %(
            current_fss_result.start_resno,current_fss_result.end_resno,
            current_fss_result.sequence)
          result_dict[sequence_string]=cg
          if not sequence_string in unique_sequence_list:
            unique_sequence_list.append(sequence_string)


    # Go through all chains and save annotation and number of good/poor h bonds
    print("\nAnalysis using %s unique sequences:" %(
       len(unique_sequence_list)), file=self.out)
    print("Unique part of the analysis:", file=self.out)
    all_sheets=[]
    all_helices=[]
    number_of_good_h_bonds=0
    number_of_poor_h_bonds=0
    import iotbx.pdb.secondary_structure as ioss
    i=0
    for sequence_string in unique_sequence_list:
      cg=result_dict[sequence_string]
      i+=1
      print(80*"=", file=self.out)
      print("\nAnalysis of chains with sequence %s: %s\n" %(
        i,sequence_string), file=self.out)
      print(80*"=", file=self.out)
      for fss_result in cg.get_fss_result_list():
        fss_result.show_summary(out=self.out)
        chain_id_list=fss_result.get_chain_id_list()
        chain_id=chain_id_list[0]
        annotation=fss_result.get_annotation().deep_copy()
        n=len(chain_id_list)
        if len(chain_id_list)>1:
          chain_id_list=chain_id_list[1:]
          annotation.multiply_to_asu_2(chain_ids_dict={chain_id:chain_id_list})
        all_helices+=annotation.helices
        all_sheets+=annotation.sheets
        number_of_good_h_bonds+=n*fss_result.number_of_good_h_bonds
        number_of_poor_h_bonds+=n*fss_result.number_of_poor_h_bonds

    self.annotation=ioss.annotation(sheets=all_sheets,helices=all_helices)
    self.annotation.renumber_helices_and_sheets()
    self.number_of_good_h_bonds=number_of_good_h_bonds
    self.number_of_poor_h_bonds=number_of_poor_h_bonds
    print(80*"=", file=self.out)
    print("\nFinal annotation and selections", file=self.out)
    print(80*"=", file=self.out)
    self.show_summary(out=self.out)



  def show_summary(self,verbose=None,pdb_records_file=None,out=sys.stdout):

    for model in self.models:
      if verbose:
        print("\nModel %d  N: %d  Start: %d End: %d" %(
          model.info.get('chain_number',0),
          model.length(),model.first_residue(),model.last_residue()), file=out)
      if verbose:
        if model.find_alpha:
          model.find_alpha.show_summary(out=out)
        if model.find_three_ten:
          model.find_three_ten.show_summary(out=out)
        if model.find_pi:
          model.find_pi.show_summary(out=out)
        if model.find_beta:
          model.find_beta.show_summary(out=out)
        if model.find_other:
          model.find_other.show_summary(out=out)
    if self.annotation and self.annotation.as_pdb_str():
      print("\nFINAL PDB RECORDS:", file=out)
      print(self.annotation.as_pdb_str(), file=out)

      if self.params.control.verbose:
        print("\n\nFINAL HELIX selections:", file=out)
        print('"%s"' %(self.annotation.overall_helix_selection()), file=out)
        print("\n\nFINAL SHEET selections:", file=out)
        print('"%s"' %(self.annotation.overall_sheet_selection()), file=out)

      print("\n\nFINAL PDB selections:", file=out)
      print('"%s"' %(self.annotation.overall_selection()), file=out)


    if pdb_records_file and self.annotation:
      f=open(pdb_records_file,'w')
      print(self.annotation.as_pdb_str(), file=f)
      f.close()
      print("\nRecords written to %s\n" %(
         pdb_records_file), file=out)

  def get_results(self):
    return self.get_annotation()

  def get_annotation(self):
    if hasattr(self,'annotation'):
      return self.annotation

  def get_user_ss(self,params=None,hierarchy=None,
     user_annotation_text=None,out=sys.stdout):

    if not user_annotation_text:
      file_name=params.input_files.secondary_structure_input
      if file_name and not os.path.isfile(file_name):
       raise Sorry("The secondary_structure_input file '%s' is missing" %(
         str(file_name)))

      # Read ss structure from this file now
      print("\nReading secondary structure records from %s\n" %(file_name), file=out)
      user_annotation_text=open(file_name).read()

    import iotbx.pdb.secondary_structure as ioss
    user_annotation=ioss.annotation.from_records(
        records=flex.split_lines(user_annotation_text))

    print("\nUser helix/strand records as input:\n", file=out)
    print(user_annotation.as_pdb_str(), file=out)
    if params.input_files.force_secondary_structure_input:
      if params.find_ss_structure.combine_annotations:
        print("\nThis secondary structure annotation will be taken as is and"+\
         " then will be \ncombined with an edited version (updating H-bonding)", file=out)
      else:
        print("\nThis secondary structure annotation will be used as is.\n", file=out)
      remove_overlaps=False
    else:
      print("\nThis secondary structure annotation will be modified if necessary\n", file=out)
      remove_overlaps=True
    # Remove any parts of this annotation that do not exist in the hierarchy
    user_annotation=remove_bad_annotation(
        user_annotation,
        hierarchy=hierarchy,
        max_h_bond_length=params.find_ss_structure.max_h_bond_length,
        remove_overlaps=remove_overlaps,
        out=out)
    if not user_annotation:
        return None

    if params.control.verbose or \
        (not params.input_files.force_secondary_structure_input) or \
        params.find_ss_structure.combine_annotations:
      local_out=out
    else:
      from libtbx.utils import null_out
      local_out=null_out()

    # Set up our alpha_helix_list etc from this...(just copy)

    # Helix classes:   'alpha', 'pi', '3_10',
    for helix in user_annotation.helices:
      ph=apply_atom_selection(
       get_string_or_first_element_of_list(helix.as_atom_selections()),
       hierarchy=hierarchy)
      model=model_info(hierarchy=ph,info={'class':helix.helix_class})
      self.user_models.append(model)
      if helix.helix_class=='alpha':
        self.user_helix_strand_segments.pdb_alpha_helix_list.append(helix)
        model.find_alpha=find_helix(params=params.alpha,model_as_segment=True,
          model=model,verbose=params.control.verbose)
      elif helix.helix_class=='pi':
        self.user_helix_strand_segments.pdb_pi_helix_list.append(helix)
        model.find_pi=find_helix(params=params.three_ten,
          model_as_segment=True,
          model=model,verbose=params.control.verbose)
      elif helix.helix_class=='3_10':
        self.user_helix_strand_segments.pdb_three_ten_helix_list.append(helix)
        model.find_three_ten=find_helix(params=params.pi,
          model_as_segment=True,
          model=model,verbose=params.control.verbose)
      else:
        raise Sorry("Unknown helix type: '%s'" %(helix.helix_class))
      self.user_helix_strand_segments.add_from_model(model)

    self.user_helix_strand_segments.sheet_list=[]
    for sheet in user_annotation.sheets:
      strand_id_in_sheet=[]
      self.user_helix_strand_segments.sheet_list.append(strand_id_in_sheet)
      prev_strand_as_segment=None
      prev_strand_id=None
      if len(sheet.registrations)!=len(sheet.strands):
        raise Sorry("\nNot 1:1 registrations (%d) to strands (%d) " %(
          len(sheet.registrations),len(sheet.strands)))

      prev_strand=None
      prev_hierarchy=None
      for strand,registration in zip(sheet.strands,sheet.registrations):
        if prev_strand:
          is_parallel=(
             (prev_strand.sense==0 and strand.sense==1) or
             (strand.sense==prev_strand.sense)
          )
        else:
          is_parallel=None

        ph=apply_atom_selection(
         get_string_or_first_element_of_list(strand.as_atom_selections()),
         hierarchy=hierarchy)

        model=model_info(hierarchy=ph,info={'class':'strand'})
        self.user_models.append(model)
        model.find_beta=find_beta_strand(params=params.beta,
          model_as_segment=True,
          model=model,verbose=params.control.verbose)
        n_strands_prev=len(self.user_helix_strand_segments.all_strands)
        self.user_helix_strand_segments.add_from_model(model)
        n_strands_cur=len(self.user_helix_strand_segments.all_strands)
        assert n_strands_cur==n_strands_prev+1

        current_strand_as_segment=\
              self.user_helix_strand_segments.all_strands[-1]
        current_strand_id=len(self.user_helix_strand_segments.all_strands)-1
        strand_id_in_sheet.append(current_strand_id)

        if prev_strand is not None: # add entries to pair_dict
          first_ca_1=0
          last_ca_1=prev_strand_as_segment.length()-1
          first_ca_2=0
          last_ca_2=current_strand_as_segment.length()-1

          if not prev_strand_id in \
              self.user_helix_strand_segments.pair_dict:
            self.user_helix_strand_segments.pair_dict[prev_strand_id]=[]
          if not current_strand_id in \
              self.user_helix_strand_segments.pair_dict:
            self.user_helix_strand_segments.pair_dict[current_strand_id]=[]

          # identify residues i_index,j_index that are next to each other
          first_last_1_and_2=\
             [first_ca_1,last_ca_1,first_ca_2,last_ca_2,is_parallel]
          i_index,j_index=self.user_helix_strand_segments.get_ind_h_bond_sheet(
              first_last_1_and_2=first_last_1_and_2,
              i=prev_strand_id,j=current_strand_id,switch_i_j=False,
              registration=registration,
              sense=strand.sense,
              force_secondary_structure_input=\
                params.input_files.force_secondary_structure_input)
          self.user_helix_strand_segments.pair_dict[prev_strand_id].append(
             current_strand_id)
          key12="%d:%d" %(prev_strand_id,current_strand_id)
          self.user_helix_strand_segments.info_dict[key12]=\
            [first_ca_1,last_ca_1,first_ca_2,last_ca_2,is_parallel,i_index,j_index]

          # Now reversed
          i_index,j_index=self.user_helix_strand_segments.get_ind_h_bond_sheet(
              first_last_1_and_2=first_last_1_and_2,
              i=prev_strand_id,j=current_strand_id,switch_i_j=True,
              registration=registration,
              sense=strand.sense,
              force_secondary_structure_input=\
                params.input_files.force_secondary_structure_input)

          self.user_helix_strand_segments.pair_dict[current_strand_id].append(
             prev_strand_id)
          key21="%d:%d" %(current_strand_id,prev_strand_id)
          self.user_helix_strand_segments.info_dict[key21]=\
            [first_ca_2,last_ca_2,first_ca_1,last_ca_1,is_parallel,i_index,j_index]
        prev_strand=strand
        prev_strand_as_segment=self.user_helix_strand_segments.all_strands[-1]
        prev_strand_id=len(self.user_helix_strand_segments.all_strands)-1
        prev_hierarchy=ph

    self.user_number_of_good_h_bonds,self.user_number_of_poor_h_bonds=\
     self.user_helix_strand_segments.set_up_pdb_records(models=self.user_models,
         max_h_bond_length=params.find_ss_structure.max_h_bond_length,
         force_secondary_structure_input=\
                params.input_files.force_secondary_structure_input,
         allow_ca_only_model=params.beta.allow_ca_only_model,out=local_out)
    print("\nNumber of good H-bonds: %d  Number of poor H-bonds: %d" %(
          self.user_number_of_good_h_bonds,self.user_number_of_poor_h_bonds), file=out)
    edited_annotation=self.user_helix_strand_segments.get_annotation()


    if self.user_helix_strand_segments.get_all_pdb_records():
      if params.input_files.force_secondary_structure_input:
        print("\nWorking PDB RECORDS (equivalent to input records):", file=out)
      else:
        print("\nInput PDB RECORDS as modified:", file=out)
      print(edited_annotation.as_pdb_str(), file=out)

    if params.find_ss_structure.combine_annotations:
      print("\nMerging input and edited annotation", file=out)
      edited_annotation=self.user_helix_strand_segments.get_annotation()
      print("\nEdited annotation:", file=out)
      print(edited_annotation.as_pdb_str(), file=out)

      print("\nUser annotation:", file=out)
      print(user_annotation.as_pdb_str(), file=out)

      combined_annotation=edited_annotation.combine_annotations(
        hierarchy=hierarchy, other=user_annotation) # will take edited if equal
      if combined_annotation:
        print("\nMerged annotation (input and edited input annotation):\n", file=out)
        print(combined_annotation.as_pdb_str(), file=out)

      return combined_annotation
    else:
      return edited_annotation

  def find_ss_in_model(self,params=None,model=None,out=sys.stdout):

    previously_used_residues=[]
    if params.find_ss_structure.find_alpha:
      model.find_alpha=find_helix(params=params.alpha,
        model=model,verbose=params.control.verbose,
        make_unique=params.find_ss_structure.make_unique,
        extract_segments_from_pdb=params.extract_segments_from_pdb.extract,
        cut_up_segments=params.find_ss_structure.cut_up_segments,
        extend_segments=params.find_ss_structure.extend_segments,
        previously_used_residues=previously_used_residues,
        out=out)
      if params.find_ss_structure.exclude_alpha_in_beta:
        previously_used_residues+=model.find_alpha.get_used_residues_list()

    if params.find_ss_structure.find_three_ten:
      model.find_three_ten=find_helix(params=params.three_ten,
        model=model,verbose=params.control.verbose,
        make_unique=params.find_ss_structure.make_unique,
        extract_segments_from_pdb=params.extract_segments_from_pdb.extract,
        cut_up_segments=params.find_ss_structure.cut_up_segments,
        extend_segments=params.find_ss_structure.extend_segments,
        previously_used_residues=previously_used_residues,
        out=out)
      if params.find_ss_structure.exclude_alpha_in_beta:
        previously_used_residues+=model.find_three_ten.get_used_residues_list()

    if params.find_ss_structure.find_pi:
      model.find_pi=find_helix(params=params.pi,
        model=model,verbose=params.control.verbose,
        make_unique=params.find_ss_structure.make_unique,
        extract_segments_from_pdb=params.extract_segments_from_pdb.extract,
        cut_up_segments=params.find_ss_structure.cut_up_segments,
        extend_segments=params.find_ss_structure.extend_segments,
        previously_used_residues=previously_used_residues,
        out=out)
      if params.find_ss_structure.exclude_alpha_in_beta:
        previously_used_residues+=model.find_pi.get_used_residues_list()

    self.beta_list_by_model=[]
    if params.find_ss_structure.find_beta:

      model.find_beta=find_beta_strand(params=params.beta,
        model=model,verbose=params.control.verbose,
        make_unique=params.find_ss_structure.make_unique,
        extract_segments_from_pdb=params.extract_segments_from_pdb.extract,
        cut_up_segments=params.find_ss_structure.cut_up_segments,
        extend_segments=params.find_ss_structure.extend_segments,
        previously_used_residues=previously_used_residues,
        out=out)
      if params.find_ss_structure.exclude_alpha_in_beta:
        previously_used_residues+=model.find_beta.get_used_residues_list()

    if params.find_ss_structure.make_unique: # do this before finding other
      # get unique residues in alpha, beta if desired
      self.make_unique(model)

    if params.find_ss_structure.find_other:
      model.find_other=find_other_structure(params=params.other,
        model=model,verbose=params.control.verbose,
        make_unique=params.find_ss_structure.make_unique,
        extract_segments_from_pdb=params.extract_segments_from_pdb.extract,
        cut_up_segments=params.find_ss_structure.cut_up_segments,
        extend_segments=params.find_ss_structure.extend_segments,
        previously_used_residues=previously_used_residues,
        out=out)

  def make_unique(self,model):
    # iteratively work down residues in this model starting with helix, then
    #  strand, then anything left
    # trim in from end if overlapping with previous, delete to nearest end
    # if a used residue in the middle

    is_used_list=model.length()*[None]
    first_res=model.first_residue()

    for find_group in ['find_alpha','find_three_ten','find_pi',
         'find_beta','find_other']:
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
        for i in range(first_pos+new_start,first_pos+new_end+1):
          is_used_list[i]=True

        #save it
        new_segment_list.append(s) # keep it if we get this far
      fg.segments=new_segment_list


  def get_start_end(self,already_used=None):
    # goes from first available to end of available (not necessarily optimal)
    new_start=None
    new_end=None
    for i in range(len(already_used)):
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
    print("\nFind secondary structure in hierarchy", file=out)
    master_phil.format(python_object=params).show(out=out)
    return params


if __name__=="__main__":
  args=sys.argv[1:]
  fss=find_secondary_structure(args=args,out=sys.stdout)

  """
  # How to get cctbx helix/sheet objects:
  alpha_helices=fss.get_pdb_alpha_helix_list()
  sheets=fss.get_pdb_sheet_list()
  print "\nHelix Summary"
  for helix in alpha_helices:
    print helix.as_pdb_str()
  print "\nSheet Summary"
  for sheet in sheets:
    print sheet.as_pdb_str()
  """
