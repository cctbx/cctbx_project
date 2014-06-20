from __future__ import division
from phenix.command_line.simple_ncs_from_pdb import simple_ncs_from_pdb
from phenix.command_line.simple_ncs_from_pdb import ncs_master_params
from iotbx import crystal_symmetry_from_any
import iotbx.pdb.hierarchy
from scitbx.array_family import flex
from libtbx.utils import null_out
from scitbx.math import superpose
from libtbx.utils import Sorry
from libtbx.phil import parse
from mmtbx.ncs import ncs
from scitbx import matrix
from iotbx import pdb
import string
import math
import os


# parameters for manual specification of NCS - ASU mapping
master_phil = parse("""
  ncs_group_selection {
    ncs_group
      .multiple = True
      {
      master_ncs_selection = ''
        .type = str
        .help = 'Residue selection string for the complete master NCS copy'
      selection_copy = ''
        .type = str
        .help = 'Residue selection string for each NCS copy location in ASU'
        .multiple = True
    }
  }
  """)


class multimer(object):
  '''
  Reconstruction of either the biological assembly or the crystallographic asymmetric unit

  Reconstruction of the biological assembly multimer by applying
  BIOMT transformation, from the pdb file, to all chains.

  Reconstruction of the crystallographic asymmetric unit by applying
  MTRIX transformations, from the pdb file, to all chains.

  self.assembled_multimer is a pdb.hierarchy object with the multimer information

  The method write generates a PDB file containing the multimer

  since chain names string length is limited to two, this process does not maintain any
  referance to the original chain names.

  Several useful attributes:
  --------------------------
  self.write() : To write the new object to a pdb file
  self.assembled_multimer : The assembled or reconstructed object is in
  self.number_of_transforms : Number of transformations that where applied

  @author: Youval Dar (LBL )
  '''
  def __init__(self,
               reconstruction_type = 'cau',
               file_name=None,
               pdb_str=None,
               error_handle=True,
               eps=1e-3,
               round_coordinates=True):
    '''
    Arguments:
    reconstruction_type -- 'ba' or 'cau'
                           'ba': biological assembly
                           'cau': crystallographic asymmetric unit
    file_name -- the name of the pdb file we want to process.
                 a string such as 'pdb_file_name.pdb'
    pdb_str -- a string containing pdb information, when not using a pdb file
    error_handle -- True: will stop execution on improper rotation matrices
                    False: will continue execution but will replace the values
                          in the rotation matrix with [0,0,0,0,0,0,0,0,0]
    eps -- Rounding accuracy for avoiding numerical issue when when testing
           proper rotation
    round_coordinates -- round coordinates of new NCS copies, for sites_cart
                         constancy

    @author: Youval Dar (2013)
    '''
    assert file_name or pdb_str
    # Read and process the pdb file
    if file_name:
      self.pdb_input_file_name = file_name
      pdb_obj = pdb.hierarchy.input(file_name=file_name)
      pdb_inp = pdb.input(file_name=file_name)
    else:
      self.pdb_input_file_name = pdb_str
      pdb_obj = pdb.hierarchy.input(pdb_string=pdb_str)
      pdb_inp = pdb.input(lines=pdb_str)
    pdb_obj_new = pdb_obj.hierarchy.deep_copy()
    self.assembled_multimer = pdb_obj_new
    if reconstruction_type == 'ba':
      transform_info = pdb_inp.process_BIOMT_records(
        error_handle=error_handle,
        eps=eps)
      self.transform_type = 'biological_assembly'
    elif reconstruction_type == 'cau':
      transform_info = pdb_inp.process_mtrix_records(
        error_handle=error_handle,
        eps=eps)
      self.transform_type = 'crystall_asymmetric_unit'
    else:
      raise Sorry('Sorry, wrong reconstruction type is given \n' + \
                  'Reconstruction type can be: \n' + \
                  "'ba': biological assembly \n" + \
                  "'cau': crystallographic asymmetric unit \n")
    if len(pdb_obj_new.models()) > 1:
      raise Sorry('Sorry, this feature currently supports on single models ' +
                  'hierarchies')

    self.transforms_obj = ncs_group_object()
    # Read the relevant transformation matrices
    self.transforms_obj.build_ncs_obj_from_pdb_ncs(
      transform_info=transform_info,
      pdb_hierarchy_inp=pdb_obj)

    # Calculate ASU (if there are any transforms to apply)
    if self.transforms_obj.transform_to_be_used:
      self.ncs_unique_chains_ids = self.transforms_obj.ncs_copies_chains_names
      self.number_of_transforms = len(self.transforms_obj.transform_to_be_used)
      new_sites = self.transforms_obj.apply_transforms(
        ncs_coordinates = pdb_obj.hierarchy.atoms().extract_xyz(),
        round_coordinates = round_coordinates)
      # apply the transformation
      model = pdb_obj_new.models()[0]
      for tr in self.transforms_obj.transform_chain_assignment:
        key = tr.split('_')[0]
        ncs_selection = self.transforms_obj.asu_to_ncs_map[key]
        new_part = pdb_obj.hierarchy.select(ncs_selection).deep_copy()
        new_chain = iotbx.pdb.hierarchy.ext.chain()
        new_chain.id = self.transforms_obj.ncs_copies_chains_names[tr]
        for res in new_part.residue_groups():
          new_chain.append_residue_group(res.detached_copy())
        model.append_chain(new_chain)
      self.assembled_multimer.atoms().set_xyz(new_sites)


  def get_ncs_hierarchy(self):
    """
    Retrieve the Original PDB hierarchy from the ASU
    """
    ncs_selection = self.transforms_obj.ncs_atom_selection
    return self.assembled_multimer.select(ncs_selection)

  def sites_cart(self):
    """ () -> flex.vec3
    Returns the reconstructed hierarchy sites cart (atom coordinates)
    """
    return self.assembled_multimer.atoms().extract_xyz()


  def write(self,pdb_output_file_name='',crystal_symmetry=None):
    ''' (string) -> text file
    Writes the modified protein, with the added chains, obtained by the BIOMT/MTRIX
    reconstruction, to a text file in a pdb format.
    self.assembled_multimer is the modified pdb object with the added chains

    Argumets:
    pdb_output_file_name -- string. 'name.pdb'
    if no pdn_output_file_name is given pdb_output_file_name=file_name

    >>> v = multimer('name.pdb','ba')
    >>> v.write('new_name.pdb')
    Write a file 'new_name.pdb' to the current directory
    >>> v.write(v.pdb_input_file_name)
    Write a file 'copy_name.pdb' to the current directory
    '''
    input_file_name = os.path.basename(self.pdb_input_file_name)
    if pdb_output_file_name == '':
      pdb_output_file_name = input_file_name
    # Aviod writing over the original file
    if pdb_output_file_name == input_file_name:
      # if file name of output is the same as the input, add 'copy_' in front of the name
      self.pdb_output_file_name = self.transform_type + '_' + input_file_name
    else:
      self.pdb_output_file_name = pdb_output_file_name
    # we need to add crystal symmetry to the new file since it is
    # sometimes needed when calulating the R-work factor (r_factor_calc.py)
    if not crystal_symmetry:
      crystal_symmetry = crystal_symmetry_from_any.extract_from(
        self.pdb_input_file_name)
    # using the pdb hierarchy pdb file writing method
    self.assembled_multimer.write_pdb_file(file_name=self.pdb_output_file_name,
      crystal_symmetry=crystal_symmetry)


class ncs_group_object(object):

  def __init__(self):
    """
    process MTRIX, BOIMT and PHIL parameters and produce an object
    with information for NCS refinement

    Argument:
    transform_info:  an object produced from the PDB MTRIX or BIOMT records
    ncs_refinement_params: an object produced by the PHIL parameters
    """
    # summery of NCS groups selection and operators
    # [[[ncs group 1 master selection],[ncs group 2 master selection],..]]
    # ,[[[r1_1,  t1_1],[r2_1,t2_1]],[[r1_2, t1_2],[r2_2,t2_2],[r2_3,t2_3]],..]
    self.group_selection_and_transforms = [[],[]]

    # Todo:  check if all attributes are being used
    # When ncs is known, reproduce asu even if the PDB file already contains
    # the complete asu
    self.total_asu_length = None

    # maps each chain in the ncs to its copies position in the asu, and the
    # asu back to the ncs.
    # the keys are chain-id_transform-serial-number, ie A_3
    # values are flex iselection of the corresponding atoms in the asu
    self.ncs_to_asu_map = {}

    # iselection of all part in ASU that are not related via NCS operators
    self.non_ncs_region_selection = None

    # keys are items in ncs_chain_selection, values are lists of selection str
    self.ncs_to_asu_selection = {}
    # values are flex iselection of the corresponding atoms in the ncs
    self.asu_to_ncs_map = {}

    self.ncs_copies_chains_names = {}
    # map all dictionaries key to chain ID or ncs_selection
    self.map_keys_to_ncs_selection = {}
    # dictionary of transform names, same keys as ncs_to_asu_map
    self.chain_transform_assignment = {}
    self.number_of_ncs_groups = 1
    self.ncs_group_map = {}
    # map transform name (s1,s2,...) to transform object
    self.ncs_transform = {}

    # list of which transform is applied to which chain
    # (the keys of ncs_to_asu_map, asu_to_ncs_map and ncs_group)
    self.transform_chain_assignment = []

    # map transformation to atoms in ncs. keys are s+transform serial number
    # ie s1,s2.... values are list of ncs_to_asu_map keys
    self.transform_to_ncs = {}

    # flex.bool ncs_atom_selection
    # (the atom selection includes any atoms not included in any
    # ncs operation. The ncs_selection_str is not)
    self.ncs_atom_selection = None
    self.ncs_selection_str = ''

    # selection of all chains in NCS
    self.all_pdb_selection = ''
    self.ncs_chain_selection = []

    # unique identifiers
    self.model_unique_chains_ids = tuple()
    self.selection_ids = set()
    # transform application order
    self.model_order_chain_ids = []
    self.transform_to_be_used = set()

    # del: consider deleting self.selection_names_index
    # Use to produce new unique names for atom selections
    self.selection_names_index = [65,65]
    # order of transforms  - used when concatenating or separating them
    self.transform_order = []

  def build_ncs_obj_from_pdb_ncs(self,
                                 pdb_hierarchy_inp,
                                 transform_info=None,
                                 rotations = None,
                                 translations = None):
    """
    Build transforms objects and NCS <-> ASU mapping using PDB file containing
    a single NCS copy and MTRIX  or BIOMT records

    Arguments:
    pdb_hierarchy_inp : iotbx.pdb.hierarchy.input
    transform_info : an object containing MTRIX or BIOMT transformation info
    rotations : matrix.sqr 3x3 object
    translations : matrix.col 3x1 object
    """
    self.collect_basic_info_from_pdb(pdb_hierarchy_inp=pdb_hierarchy_inp)
    self.ncs_selection_str = self.all_pdb_selection
    assert bool(transform_info or (rotations and translations))
    if rotations:
      # add rotations,translations to ncs_refinement_groups
      self.add_transforms_to_ncs_refinement_groups(
        rotations=rotations,
        translations=translations)
    else:
      # use only MTRIX/BIOMT records from PDB
      assert_test = all_ncs_copies_not_present(transform_info)
      msg = 'MTRIX record error : Mixed present and not-present NCS copies'
      assert assert_test,msg
      self.process_pdb(transform_info=transform_info)

    self.transform_chain_assignment = get_transform_order(self.transform_to_ncs)

    self.ncs_copies_chains_names = self.make_chains_names(
      transform_assignment=self.transform_chain_assignment,
      unique_chain_names = self.model_unique_chains_ids)

    # build self.ncs_to_asu_selection
    for k in self.transform_chain_assignment:
      selection_str = 'chain ' + self.ncs_copies_chains_names[k]
      key =  k.split('_')[0]
      if self.ncs_to_asu_selection.has_key(key):
        self.ncs_to_asu_selection[key].append(selection_str)
      else:
        self.ncs_to_asu_selection[key] = [selection_str]

    self.finalize_pre_process(pdb_hierarchy_inp=pdb_hierarchy_inp)

  def build_ncs_obj_from_phil(self,
                              ncs_selection_params = None,
                              ncs_phil_groups = None,
                              pdb_hierarchy_inp = None):
    """
    Build transforms objects and NCS <-> ASU mapping using phil selection
    strings and complete ASU

    Arguments:
    ncs_selection_params : Phil parameters
    pdb_hierarchy_inp : iotbx.pdb.hierarchy.input

    Phil structure
    ncs_group_selection {
      ncs_group (multiple)
      {
        master_ncs_selection = ''
        selection_copy = ''   (multiple)
      }
    }
    """
    # process params
    if ncs_selection_params:
      if isinstance(ncs_selection_params,str):
        ncs_selection_params = parse(ncs_selection_params)
      phil_param =  master_phil.fetch(
        source=ncs_selection_params,track_unused_definitions=True)
      working_phil = phil_param[0].extract()
      assert  phil_param[1] == [],'Check phil parameters...'
      ncs_phil_groups = working_phil.ncs_group_selection.ncs_group
    else:
      assert not (ncs_phil_groups is None)
    assert self.ncs_selection_str == ''
    unique_selections = set()
    transform_sn = 0
    ncs_group_id = 0
    # populate ncs selection and ncs to copies location
    for group in ncs_phil_groups:
      gns = group.master_ncs_selection
      self.ncs_chain_selection.append(gns)
      unique_selections = uniqueness_test(unique_selections,gns)
      master_ncs_ph = get_pdb_selection(ph=pdb_hierarchy_inp,selection_str=gns)
      ncs_group_id += 1
      transform_sn += 1
      self.add_identity_transform(
        ncs_selection=gns,
        ncs_group_id=ncs_group_id,
        transform_sn=transform_sn)
      asu_locations = []
      for asu_select in group.selection_copy:
        unique_selections = uniqueness_test(unique_selections,asu_select)
        asu_locations.append(asu_select)
        transform_sn += 1
        key = 's' + format_num_as_str(transform_sn)
        ncs_copy_ph = get_pdb_selection(
          ph=pdb_hierarchy_inp,selection_str=asu_select)
        r, t = get_rot_trans(master_ncs_ph,ncs_copy_ph)
        tr = transform(
          rotation = r,
          translation = t,
          serial_num = transform_sn,
          coordinates_present = True,
          ncs_group_id = ncs_group_id)
        # Update ncs_group dictionary and transform_to_ncs list
        self.build_transform_dict(
          transform_id = key,
          transform = tr,
          selection_id = gns)
        self.ncs_group_map = update_ncs_group_map(
          ncs_group_map=self.ncs_group_map,
          ncs_group_id = ncs_group_id,
          selection_ids = gns,
          transform_id = key)
        assert not self.ncs_transform.has_key(key)
        self.ncs_transform[key] = tr
        self.selection_ids.add(gns)
      self.ncs_to_asu_selection[gns] = asu_locations
      self.number_of_ncs_groups = ncs_group_id

    self.ncs_selection_str = '('+ self.ncs_chain_selection[0] +')'
    for i in range(1,len(self.ncs_chain_selection)):
      self.ncs_selection_str += ' or (' + self.ncs_chain_selection[i] + ')'

    self.transform_chain_assignment = get_transform_order(self.transform_to_ncs)
    self.finalize_pre_process(pdb_hierarchy_inp=pdb_hierarchy_inp)

  def build_ncs_obj_from_pdb_asu(self,pdb_hierarchy_inp):
    """
    Build transforms objects and NCS <-> ASU mapping using PDB file which
    contains all NCS copies and MTRIX.
    Note that the MTRIX record are ignored, they are being produced in the
    process of identifying the master NCS

    Arguments:
    pdb_hierarchy_inp : iotbx.pdb.hierarchy.input
    """
    # identify the NCS copies of each group
    # Todo: evaluate if there is a simple NCS groups. replace the
    # simple_find_ncs_operators
    ncs_phil_groups = []
    n_atoms_per_chain = []
    n_atoms_per_chain_dict = {}
    # evaluate possible master NCS possibilities
    for chain in pdb_hierarchy_inp.hierarchy.chains():
      n_atoms = chain.atoms_size()
      n_atoms_per_chain.append(n_atoms)
      if n_atoms_per_chain_dict.has_key(n_atoms):
        n_atoms_per_chain_dict[n_atoms].append('chain ' + chain.id)
      else:
        n_atoms_per_chain_dict[n_atoms] = ['chain ' + chain.id]

    # Todo: make more general
    # Basic case: Assume first chain is the master

    if len(n_atoms_per_chain_dict.keys()) == 1:
      # Note: even when there all chain ae of the same length, master NCS
      # might be made of several chains
      new_ncs_group = ncs_groups_container()
      key = n_atoms_per_chain_dict.keys()[0]
      chains = n_atoms_per_chain_dict[key]
      new_ncs_group.master_ncs_selection = chains[0]
      new_ncs_group.selection_copy = chains[1:]
      ncs_phil_groups.append(new_ncs_group)

    # Todo: evaluate possible groups base on number of atoms
    # evaluate possible groups by testing rotations and build a selection
    # list that can use the phil processing method
    # build_ncs_obj_from_phil

    #
    if ncs_phil_groups != []:
      self.build_ncs_obj_from_phil(
        pdb_hierarchy_inp=pdb_hierarchy_inp,
        ncs_phil_groups=ncs_phil_groups)
    else:
      params = ncs_master_params.extract()
      params.simple_ncs_from_pdb.min_length = 1
      ncs_from_pdb=simple_ncs_from_pdb(
        pdb_inp=pdb_hierarchy_inp.input,
        hierarchy=pdb_hierarchy_inp.hierarchy,
        quiet=True,
        log=null_out(),
        params=params)
      spec_ncs_groups = ncs_from_pdb.ncs_object.ncs_groups()
      self.build_ncs_obj_from_spec_file(
        pdb_hierarchy_inp=pdb_hierarchy_inp,
        spec_ncs_groups=spec_ncs_groups)


  def build_ncs_obj_from_spec_file(self,file_name=None,file_path='',
                                   file_str='',source_info="",quiet=True,
                                   pdb_hierarchy_inp=None,
                                   spec_ncs_groups=None):
    """
    read .spec files and build transforms object and NCS <-> ASU mapping

    Arguments:
    spec file information using a file, file_str (str) of a source_info (str)
    quite: (bool) quite output when True
    pdb_hierarchy_inp: pdb hierarchy input
    spec_ncs_groups: ncs_groups object
    """
    if spec_ncs_groups is None:
      fn = os.path.join(file_path,file_name)
      lines = file_str.splitlines()
      assert lines != [] or os.path.isfile(fn)
      spec_object = ncs.ncs()
      spec_object.read_ncs(
        file_name=fn,lines=lines,source_info=source_info,quiet=quiet)
      spec_ncs_groups = spec_object.ncs_groups()

    transform_sn = 0
    ncs_group_id = 0
    for gr in spec_ncs_groups:
      ncs_group_id += 1
      # create selection
      ncs_group_selection_list = get_ncs_group_selection(gr.chain_residue_id())
      gs = ncs_group_selection_list[0]
      self.ncs_chain_selection.append(gs)
      asu_locations = []
      for i,ncs_copy_select in enumerate(ncs_group_selection_list):
        r = gr.rota_matrices()[i]
        t = gr.translations_orth()[i]
        if (not r.is_r3_identity_matrix()) and (not t.is_r3_identity_matrix()):
          asu_locations.append(ncs_copy_select)
        transform_sn += 1
        key = 's' + format_num_as_str(transform_sn)
        tr = transform(
          rotation = r,
          translation = t,
          serial_num = transform_sn,
          coordinates_present = True,
          ncs_group_id = ncs_group_id)
        # Update ncs_group dictionary and transform_to_ncs list
        self.build_transform_dict(
          transform_id = key,
          transform = tr,
          selection_id = gs)
        self.ncs_group_map = update_ncs_group_map(
          ncs_group_map=self.ncs_group_map,
          ncs_group_id = ncs_group_id,
          selection_ids = gs,
          transform_id = key)
        assert not self.ncs_transform.has_key(key)
        self.ncs_transform[key] = tr
        self.selection_ids.add(gs)
      self.ncs_to_asu_selection[gs] = asu_locations
      self.number_of_ncs_groups = ncs_group_id

    self.ncs_selection_str = '('+ self.ncs_chain_selection[0] +')'
    for i in range(1,len(self.ncs_chain_selection)):
      self.ncs_selection_str += ' or (' + self.ncs_chain_selection[i] + ')'

    self.transform_chain_assignment = get_transform_order(self.transform_to_ncs)
    self.finalize_pre_process(pdb_hierarchy_inp=pdb_hierarchy_inp)

  def build_ncs_obj_from_mmcif(self):
    """
    Build transforms objects and NCS <-> ASU mapping using mmcif file
    """
    # TODO: build function

  def add_transforms_to_ncs_refinement_groups(self,rotations,translations):
    """
    Add rotation matrices and translations vectors
    to ncs_refinement_groups
    """
    assert len(rotations) == len(translations)
    assert not self.ncs_transform, 'ncs_transform should be empty'
    sn = {1}
    self.add_identity_transform(ncs_selection=self.ncs_selection_str)
    n = 1
    assert len(rotations) == len(translations)
    for (r,t) in zip(rotations,translations):
      # check if transforms are the identity transform
      identity = r.is_r3_identity_matrix() and t.is_col_zero()
      if not identity:
        n += 1
        sn.add(n)
        key = 's' + format_num_as_str(n)
        tr = transform(
          rotation = r,
          translation = t,
          serial_num = n,
          coordinates_present = False,
          ncs_group_id = 1)
        assert not self.ncs_transform.has_key(key)
        self.ncs_transform[key] = tr
        for select in self.ncs_chain_selection:
          self.build_transform_dict(
            transform_id = key,
            transform = tr,
            selection_id = select)
          self.selection_ids.add(select)
        self.ncs_group_map = update_ncs_group_map(
          ncs_group_map=self.ncs_group_map,
          ncs_group_id = 1,
          selection_ids = self.ncs_chain_selection,
          transform_id = key)

  def collect_basic_info_from_pdb(self,pdb_hierarchy_inp):
    """  Build chain selection string and collect chains IDs from pdb """
    if pdb_hierarchy_inp:
      model  = pdb_hierarchy_inp.hierarchy.models()[0]
      chain_ids = {x.id for x in model.chains()}
      # Collect order if chains IDs and unique IDs
      self.model_unique_chains_ids = tuple(sorted(chain_ids))
      self.model_order_chain_ids = [x.id for x in model.chains()]
      s = ' or chain '.join(self.model_unique_chains_ids)
      self.all_pdb_selection = 'chain ' + s
      assert self.ncs_chain_selection == []
      self.ncs_chain_selection =\
        ['chain ' + s for s in self.model_unique_chains_ids]
      self.ncs_chain_selection.sort()

      # Note: probably not needed - used for asu entry
      # modify the code to only modify existing hierarchy
      self.ncs_copies_chains_names = [x.id for x in model.chains()]

  def compute_ncs_asu_coordinates_map(self,pdb_hierarchy_inp):
    """ Calculates coordinates maps from ncs to asu and from asu to ncs """
    temp = pdb_hierarchy_inp.hierarchy.atom_selection_cache()
    # check if pdb_hierarchy_inp contain only the master NCS copy
    pdb_length = len(pdb_hierarchy_inp.hierarchy.atoms())
    self.ncs_atom_selection = temp.selection(self.ncs_selection_str)
    ncs_length = self.ncs_atom_selection.count(True)
    # keep track on the asu copy number
    copy_count = {}

    if pdb_length > ncs_length:
      self.total_asu_length = pdb_length
      selection_ref = flex.bool([False]*pdb_length)
      for k in self.transform_chain_assignment:
        key =  k.split('_')[0]
        ncs_selection = temp.selection(key)
        if not self.asu_to_ncs_map.has_key(key):
          copy_count[key] = 0
          selection_ref = update_selection_ref(selection_ref,ncs_selection)
          self.asu_to_ncs_map[key] = ncs_selection.iselection()
        else:
          copy_count[key] += 1
        asu_copy_ref = self.ncs_to_asu_selection[key][copy_count[key]]
        asu_selection = temp.selection(asu_copy_ref)
        selection_ref = update_selection_ref(selection_ref,asu_selection)
        self.ncs_to_asu_map[k] = asu_selection.iselection()
      self.non_ncs_region_selection = (~selection_ref).iselection()
      # add the non ncs regions to the master ncs copy
      self.ncs_atom_selection = self.ncs_atom_selection | (~selection_ref)
    elif pdb_length == ncs_length:
      # this case is when the pdb hierarchy contain only the master NCS copy
      self.total_asu_length = self.get_asu_length(temp)
      ns = [True]*pdb_length + [False]*(self.total_asu_length - pdb_length)
      self.ncs_atom_selection = flex.bool(ns)
      sorted_keys = sorted(self.transform_to_ncs)
      i = 0
      for k in sorted_keys:
        v = self.transform_to_ncs[k]
        i += 1
        for transform_key in v:
          key =  transform_key.split('_')[0]
          ncs_basic_selection = temp.selection(key)
          ncs_selection = flex.bool(self.total_asu_length,temp.iselection(key))
          # asu_key = self.ncs_copies_chains_names[transform_key]
          if not self.asu_to_ncs_map.has_key(key):
            self.asu_to_ncs_map[key] = ncs_selection.iselection()
          # make the selection at the proper location at the ASU
          selection_list = list(ncs_basic_selection)
          temp_selection = flex.bool([False]*i*ncs_length + selection_list)
          asu_selection =flex.bool(self.total_asu_length,
                                   temp_selection.iselection())
          self.ncs_to_asu_map[transform_key] = asu_selection.iselection()

    for select,tr in self.ncs_group_map.itervalues():
      if len(select) == 1:
        s = select.pop()
      else:
        s = ' or '.join(select)
      ncs_group_selection = temp.selection(s)
      selection = ncs_group_selection.iselection()
      nt = self.ncs_transform
      transforms = [[nt[x].r, nt[x].t] for x in tr]
      self.group_selection_and_transforms[0].append(selection)
      self.group_selection_and_transforms[1].append(transforms)


  def add_identity_transform(self,ncs_selection,ncs_group_id=1,transform_sn=1):
    """    Add identity transform
    Argument:

    ncs_selection: (str) selection string for the NCS master copy
    ncs_group_id: (int) the NCS group ID
    transform_sn: (int) Over all transform serial number
    """
    transform_obj = transform(
      rotation = matrix.sqr([1,0,0,0,1,0,0,0,1]),
      translation = matrix.col([0,0,0]),
      serial_num = transform_sn,
      coordinates_present = True,
      ncs_group_id = ncs_group_id)
    id_str = 's' + format_num_as_str(transform_sn)
    self.ncs_transform[id_str] = transform_obj
    self.build_transform_dict(
      transform_id = id_str,
      transform = transform_obj,
      selection_id = ncs_selection)
    self.selection_ids.add(ncs_selection)
    self.ncs_group_map = update_ncs_group_map(
      ncs_group_map=self.ncs_group_map,
      ncs_group_id = ncs_group_id,
      selection_ids = ncs_selection,
      transform_id = id_str)

  def process_pdb(self,transform_info):
    """
    Process PDB Hierarchy object

    Remarks:
    The information on a chain in a PDB file does not have to be continuous.
    Every time the chain name changes in the pdb file, a new chain is added
    to the model, even if the chain ID already exist. so there model.
    chains() might contain several chains that have the same chain ID
    """

    transform_info_available = bool(transform_info) and bool(transform_info.r)
    if transform_info_available:
      ti = transform_info
      for (r,t,n,cp) in zip(ti.r,ti.t,ti.serial_number,ti.coordinates_present):
        n = int(n)
        key = 's' + format_num_as_str(n)
        tr = transform(
          rotation = r,
          translation = t,
          serial_num = n,
          coordinates_present = cp,
          ncs_group_id = 1)
        for select in self.ncs_chain_selection:
          self.build_transform_dict(
            transform_id = key,
            transform = tr,
            selection_id = select)
          self.selection_ids.add(select)
        self.ncs_group_map = update_ncs_group_map(
          ncs_group_map=self.ncs_group_map,
          ncs_group_id = 1,
          selection_ids = self.ncs_chain_selection,
          transform_id = key)
        # if ncs selection was not provided in phil parameter
        assert not self.ncs_transform.has_key(key)
        self.ncs_transform[key] = tr

  def build_transform_dict(self,
                           transform_id,
                           transform,
                           selection_id):
    """
    Apply all non-identity transforms

    Arguments:
    transform_id : (str) s001,s002...
    transform : transform object, containing information on transformation
    selection_id : (str) NCS selection string

    Build transform_to_ncs dictionary, which provides the location of the
    particular chains or selections in the NCS
    """
    if not (transform.r.is_r3_identity_matrix() and transform.t.is_col_zero()):
      self.transform_to_be_used.add(transform.serial_num)
      key = selection_id + '_' + format_num_as_str(transform.serial_num)
      self.chain_transform_assignment[key] = transform_id
      if self.transform_to_ncs.has_key(transform_id):
        self.transform_to_ncs[transform_id].append(key)
      else:
        self.transform_to_ncs[transform_id] = [key]

  def get_asu_length(self,atom_selection_cache):
    """" Collect the length of all ncs copies """
    asu_total_length = self.ncs_atom_selection.count(True)
    for k in self.transform_chain_assignment:
      key =  k.split('_')[0]
      ncs_selection = atom_selection_cache.selection(key)
      asu_total_length += ncs_selection.count(True)
    return asu_total_length

  def build_MTRIX_object(self):
    """
    Build a MTRIX object from ncs_group_object
    Used for testing
    """
    assert  self.number_of_ncs_groups == 1
    result = iotbx.pdb._._mtrix_and_biomt_records_container()
    tr_dict = self.ncs_transform
    tr_sorted = sorted(tr_dict,key=lambda k:tr_dict[k].serial_num)
    for key in tr_sorted:
      transform = self.ncs_transform[key]
      result.add(
        r=transform.r,
        t=transform.t,
        coordinates_present=transform.coordinates_present,
        serial_number=transform.serial_num)
    return result

  def produce_selection_name(self):
    """
    Create a new unique name each time called, to name the NCS selections,
    since they do not have to be the chain names
    """
    top_i = ord('Z')
    i,j = self.selection_names_index
    new_selection_name = 'S' + chr(i) + chr(j)
    i += (j == ord('Z')) * 1
    j = 65 + (j - 65 + 1) % 26
    assert i <= top_i
    self.selection_names_index = [i,j]
    return new_selection_name

  def make_chains_names(self,
                        transform_assignment,
                        unique_chain_names):
    """
    Create a dictionary names for the new NCS copies
    keys: (str) chain_name + '_s' + serial_num
    values: (str) (one or two chr long)

    Chain names might repeat themselves several times in a pdb file
    We want copies of chains with the same name to still have the
    same name after similar BIOMT/MTRIX transformation

    Arguments:
    transform_assignment : (list) transformation assignment order
    unique_chain_names : (tuple) a set of unique chain names

    Returns:
    new_names : a dictionary. {'A_1': 'G', 'A_2': 'H',....} map a chain
    name and a transform number to a new chain name

    >>> self.make_chains_names((1,2),['chain A_002','chain B_002'],('A','B'))
    {'A_1': 'C', 'A_2': 'D', 'B_1': 'E', 'B_2': 'F'}
    """
    if not transform_assignment or not unique_chain_names: return {}
    # create list of character from which to assemble the list of names
    # total_chains_number = len(i_transforms)*len(unique_chain_names)
    total_chains_number = len(transform_assignment)
    # start naming chains with a single letter
    chr_list1 = list(set(string.ascii_uppercase) - set(unique_chain_names))
    chr_list2 = list(set(string.ascii_lowercase) - set(unique_chain_names))
    chr_list1.sort()
    chr_list2.sort()
    new_names_list = chr_list1 + chr_list2
    # check if we need more chain names
    if len(new_names_list) < total_chains_number:
      n_names =  total_chains_number - len(new_names_list)
      # the number of character needed to produce new names
      chr_number = int(math.sqrt(n_names)) + 1
      # build character list
      chr_list = list(string.ascii_uppercase) + \
                 list(string.ascii_lowercase) + \
                 list(string.digits)
      # take only as many characters as needed
      chr_list = chr_list[:chr_number]
      extra_names = set([ x+y for x in chr_list for y in chr_list])
      # make sure not using existing names
      extra_names = list(extra_names - unique_chain_names)
      extra_names.sort()
      new_names_list.extend(extra_names)
    assert len(new_names_list) >= total_chains_number
    dictionary_values = new_names_list[:total_chains_number]
    # create the dictionary
    zippedlists = zip(transform_assignment,dictionary_values)
    new_names_dictionary = {x:y for (x,y) in zippedlists}
    return new_names_dictionary

  def apply_transforms(self,
                      ncs_coordinates,
                      round_coordinates = True):
    """
    Apply transformation to the information in the transforms_obj to the
    ncs_coordinates (flex.vec3), and round the results if
    round_coordinates is True

    Returns:
    complete asymmetric or the biological unit
    """
    asu_xyz = flex.vec3_double([(0,0,0)]*self.total_asu_length)
    s = flex.size_t_range(len(ncs_coordinates))
    asu_xyz.set_selected(s,ncs_coordinates)

    for trans in self.transform_chain_assignment:
      asu_selection = self.ncs_to_asu_map[trans]
      ncs_selection = self.asu_to_ncs_map[trans.split('_')[0]]
      ncs_xyz = ncs_coordinates.select(ncs_selection)
      tr_key = 's' + trans.split('_')[1]
      assert self.ncs_transform.has_key(tr_key)
      tr = self.ncs_transform[tr_key]
      new_sites = tr.r.elems* ncs_xyz+tr.t
      asu_xyz.set_selected(asu_selection,new_sites)
    if round_coordinates:
      return flex.vec3_double(asu_xyz).round(3)
    else:
      return flex.vec3_double(asu_xyz)

  def identify_ncs(self,transform_info):
    """
    Find the ncs group and create total NCS selection string

    Arguments:
    transform_info : an object containing MTRIX
    """
    # Todo: add find NCS from PDB ASU
    # Todo: update self.ncs_chain_selection to include only the NCS chains selection
    params = ncs_master_params.extract()
    params.simple_ncs_from_pdb.min_length = 1
    ncs_from_pdb=simple_ncs_from_pdb(
    pdb_inp=pdb_obj.input,
    hierarchy=pdb_obj.hierarchy,
    quiet=True,
    log=null_out(),
    params=params)

    # update ncs_chain_selection to include only master NCS chains
    self.ncs_chain_selection =['chain ' + s for s in ncs_chain_ids]
    self.ncs_chain_selection.sort()
    # set str selection
    self.ncs_selection_str = '('+ self.ncs_chain_selection[0] +')'
    for i in range(1,len(self.ncs_chain_selection)):
      self.ncs_selection_str += ' or (' + self.ncs_chain_selection[i] + ')'


  def finalize_pre_process(self,pdb_hierarchy_inp=None):
    """
    Steps that are common to most method of transform info
    """
    if pdb_hierarchy_inp:
      # get the mapping between NCS and ASU
      self.map_keys_to_ncs_selection=\
        {k:k.split('_')[0] for k in self.chain_transform_assignment}
      self.compute_ncs_asu_coordinates_map(pdb_hierarchy_inp=pdb_hierarchy_inp)

    self.transform_order = sorted(self.transform_to_ncs)

class transform(object):

  def __init__(self,
               rotation,
               translation,
               serial_num,
               coordinates_present,
               ncs_group_id):
    """
    Basic transformation properties

    Argument:
    rotation : Rotation matrix object
    translation: Translation matrix object
    serial_num : (int) Transform serial number
    coordinates_present: equals 1 when coordinates are presents in PDB file
    ncs_group_id : (int) ncs groups in, all the selections that the same
                   transforms are being applied to
    """
    self.r = rotation
    self.t = translation
    self.serial_num = serial_num
    self.coordinates_present = bool(coordinates_present)
    self.ncs_group_id = ncs_group_id

def format_num_as_str(n):
  """  return a 3 digit string of n  """
  if n > 999 or n < 0:
    raise IOError('Input out of the range 0 - 999.')
  else:
    s1 = n//100
    s2 = (n%100)//10
    s3 = n%10
    return str(s1) + str(s2) + str(s3)

def all_ncs_copies_not_present(transform_info):
  """
  Check if a single NCS is in transform_info, all transforms, but the
  identity should no be present
  :param transform_info: (transformation object)
  :return: (bool)
  """
  test = False
  ti = transform_info
  for (r,t,n,cp) in zip(ti.r,ti.t,ti.serial_number,ti.coordinates_present):
    if not (r.is_r3_identity_matrix() and t.is_col_zero()):
      test = test or cp
  return not test

def all_ncs_copies_present(transform_info):
  """
  Check if all transforms coordinates are present,
  if the complete ASU is present
  :param transform_info: (transformation object)
  :return: (bool)
  """
  test = True
  for cp in transform_info.coordinates_present:
    test = test and cp
  return test

def order_transforms(transform_to_ncs):
  """ set transformation application order  """
  return sorted(transform_to_ncs,key=lambda key: int(key[1:]))

def uniqueness_test(unique_selection_set,new_item):
  """
  Insert new item to set, if not there, raise an error if already in set

  :param unique_selection_set: (set)
  :param new_item: (str)
  :return: updated set
  """
  if new_item in unique_selection_set:
    raise IOError,'Phil selection strings are not unique !!!'
  else:
    unique_selection_set.add(new_item)
    return unique_selection_set

def get_rot_trans(master_ncs_ph,ncs_copy_ph):
  """
  Get rotation and translation using superpose

  :param ph: pdb hierarchy input object
  :return r: rotation matrix
  :return t: translation vector
  """
  t1 = not (master_ncs_ph is None)
  t2 = not (ncs_copy_ph is None)
  if t1 and t2:
    test = len(master_ncs_ph.atoms()) == len(ncs_copy_ph.atoms())
    assert test,'Master and copy lengths are different. Check phil parameters..'
    lsq_fit_obj = superpose.least_squares_fit(
      reference_sites = ncs_copy_ph.atoms().extract_xyz(),
      other_sites     = master_ncs_ph.atoms().extract_xyz())
    r = lsq_fit_obj.r
    t = lsq_fit_obj.t
  else:
    r = matrix.sqr([0]*9)
    t = matrix.col([0,0,0])
  return r,t

def get_pdb_selection(ph,selection_str):
  """
  :param ph: pdb hierarchy input object
  :param selection_str:
  :return: Portion of the hierarchy according to the selection
  """
  if not (ph is None):
      temp = ph.hierarchy.atom_selection_cache()
      return ph.hierarchy.select(temp.selection(selection_str))
  else:
    return None

def update_selection_ref(selection_ref,new_selection):
  """
  Test for overlapping selection and then updates the selection_ref
  with the new_selection

  Both received and return arguments are flex.bool
  """
  test = (selection_ref & new_selection).count(True) == 0
  assert test,'Overlapping atom selection. Check phil parameters..'
  return selection_ref | new_selection

def get_ncs_group_selection(chain_residue_id):
  """
  :param chain_residue_id: [[chain id's],[[[residues range]],[[...]]]
  :return: selection lists, with selection string for each ncs copy in the group
  """
  chains = chain_residue_id[0]
  res_ranges = chain_residue_id[1]
  assert len(chains) == len(res_ranges)
  ncs_selection = []
  for c,rr in zip(chains, res_ranges):
    c = c.strip()
    assert c.find(' ') < 0,'Multiple chains in a single spec ncs group'
    ch_selection = 'chain ' + c
    res_range = ['resseq {0}:{1}'.format(s,e) for s,e in rr]
    res_range = '(' + ' or '.join(res_range) + ')'
    ncs_selection.append(ch_selection + ' and ' + res_range)
  return ncs_selection

def get_transform_order(transform_to_ncs):
  """ order transforms mainly for proper chain naming """
  transform_order = order_transforms(transform_to_ncs)
  transform_chain_assignment = []
  for tr_id in transform_order:
    for tr_selection in transform_to_ncs[tr_id]:
      transform_chain_assignment.append(tr_selection)
  return transform_chain_assignment

def update_ncs_group_map(
        ncs_group_map, ncs_group_id, selection_ids, transform_id):
  """  Update ncs_group_map  """
  if isinstance(selection_ids, str): selection_ids = [selection_ids]
  if ncs_group_map.has_key(ncs_group_id):
    ncs_group_map[ncs_group_id][0].update(set(selection_ids))
    ncs_group_map[ncs_group_id][1].add(transform_id)
  else:
    ncs_group_map[ncs_group_id] = [set(selection_ids),{transform_id}]
  return ncs_group_map


class ncs_groups_container(object):

  def __init__(self):
    self.master_ncs_selection = ''
    self.selection_copy = []
