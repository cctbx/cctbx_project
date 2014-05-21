from __future__ import division
from iotbx import crystal_symmetry_from_any
import iotbx.pdb.hierarchy
from scitbx.array_family import flex
from libtbx.utils import Sorry
from libtbx.phil import parse
from scitbx import matrix
from iotbx import pdb
import string
import math
import os

master_phil = parse("""
  ncs_refinement {
    dont_apply_when_coordinates_present = True
      .type = bool
      .help='''
      Can apply transformations even if coordinates are already present.
      Need to provide NCS selection'''
    ncs_selection = None
      .type = str
      .help = '''
      When applying transforms to coordinates that are present or when
      applying selected transformations to selected chains, you need to
      provide selection string for the NCS portion of the pdb hierarchy'''
    ncs_group
      .multiple = True
      {
        transform
        .multiple = True
        {
          rotation = None
            .type = floats(size=9, value_min=-1, value_max=1)
          translation = None
            .type = floats(size=3)
          coordinates_present = None
            .type = bool
        }
        apply_to_selection = None
          .help = '''Selection syntax similar to what is used in CNS or PyMOL.
          multiple selection lines will be concatenates.'''
          .type = str
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
               reconstruction_type,
               file_name=None,
               pdb_str=None,
               error_handle=True,
               eps=1e-3,
               round_coordinates=True,
               params=None ):
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
    params -- information about how to apply transformations, with the
              following phil structure:
                ncs_refinement {
                  dont_apply_when_coordinates_present = (bool)
                  ncs_selection = (str)
                  ncs_group {
                    transform {
                      rotation = tuple(float), size 9
                      translation = tuple(float), size 3
                      coordinates_present = (Bool / None)
                    }
                    apply_to_selection = (bool)
                  }
                }

    @author: Youval Dar (2013)
    '''
    assert file_name or pdb_str
    # Read and process the pdb file
    if file_name:
      self.pdb_input_file_name = file_name
      pdb_inp = pdb.input(file_name=file_name)
      pdb_obj = pdb.hierarchy.input(file_name=file_name)
    else:
      self.pdb_input_file_name = pdb_str
      pdb_inp = pdb.input(source_info=None, lines=pdb_str)
      pdb_obj = pdb.hierarchy.input(pdb_string=pdb_str)
    pdb_obj_new = pdb_obj.hierarchy.deep_copy()
    self.assembled_multimer = pdb_obj_new

    # Read the relevant transformation matrices
    if reconstruction_type == 'ba':
      # Read BIOMT info
      transform_info = pdb_inp.process_BIOMT_records(
        error_handle=error_handle,eps=eps)
      self.transform_type = 'biological_assembly'
    elif reconstruction_type == 'cau':
      # Read MTRIX info
      transform_info = pdb_inp.process_mtrix_records(
        error_handle=error_handle,eps=eps)
      self.transform_type = 'crystall_asymmetric_unit'
    else:
      raise Sorry('Worg reconstruction type is given \n' + \
                  'Reconstruction type can be: \n' + \
                  "'ba': biological assembly \n" + \
                  "'cau': crystallographic asymmetric unit \n")
    if len(pdb_obj_new.models()) > 1:
      raise Sorry('Sorry, this feature currently supports on single models ' +
                  'hierarchies')

    self.transforms_obj = ncs_group_object()
    self.transforms_obj.populate_ncs_group_object(
      ncs_refinement_params = params,
      transform_info = transform_info,
      pdb_hierarchy_inp = pdb_obj)
    # Calculate ASU (if there are any transforms to apply)
    if self.transforms_obj.transform_to_be_used:
      self.ncs_unique_chains_ids = self.transforms_obj.model_unique_chains_ids
      self.number_of_transforms = len(self.transforms_obj.transform_to_be_used)
      new_sites = self.transforms_obj.apply_transforms(
        ncs_coordinates = pdb_obj.hierarchy.atoms().extract_xyz(),
        round_coordinates = round_coordinates)
      # apply the transformation
      model = pdb_obj_new.models()[0]
      for tr in self.transforms_obj.transform_chain_assignment:
        ncs_selection = \
          self.transforms_obj.asu_to_ncs_map[tr.split('_')[0]]
        new_part = pdb_obj.hierarchy.select(ncs_selection).deep_copy()
        new_chain = iotbx.pdb.hierarchy.ext.chain()
        new_chain.id = self.transforms_obj.ncs_copies_chains_names[tr]
        for res in new_part.residue_groups():
          new_chain.append_residue_group(res.detached_copy())
        model.append_chain(new_chain)
      self.assembled_multimer.atoms().set_xyz(new_sites)


  def get_source_hierarchy(self):
    """
    Retrieve the Original PDB hierarchy from the ASU
    """
    ncs_length = self.transforms_obj.total_length_extended_ncs
    asu_length = self.transforms_obj.total_asu_length
    source_atom_selection = flex.bool([True]*ncs_length)
    temp = flex.bool([False]*(asu_length - ncs_length))
    source_atom_selection.extend(temp)
    return self.assembled_multimer.select(source_atom_selection)

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
    self.ncs_refinement_groups = None
    self.apply_to_all_chains = True
    # When ncs is known, reproduce asu even if the PDB file already contains
    # the complete asu
    self.dont_apply_when_coordinates_present = True

    # maps each chain in the ncs to its copies position in the asu, and the
    # asu back to the ncs.
    # the keys are chain-id_transform-serial-number, ie A_3
    # values are list [start pos, end pos]
    self.ncs_to_asu_map = {}
    # values are flex iselection of the corresponding atoms in the ncs
    self.asu_to_ncs_map = {}
    self.ncs_copies_chains_names = {}
    # map all dictionaries key to chain ID or ncs_selection
    self.map_keys_to_selection = {}
    # dictionary of transform names, same keys as ncs_to_asu_map
    # TODO: Consider changing to list
    self.ncs_group = {}
    self.ncs_group_pdb = {}
    # map transform name (s1,s2,...) to transform object
    self.ncs_transform = {}
    self.ncs_transform_pdb = {}

    # dictionary kes: selection names. values: number_of_copies_in_asu
    self.number_of_ncs_copies = {}

    # list of which transform is applied to which chain
    # (the keys of ncs_to_asu_map, asu_to_ncs_map and ncs_group)
    self.transform_chain_assignment = []

    # map transformation to atoms in ncs. keys are s+transform serial number
    # ie s1,s2.... values are list of ncs_to_asu_map keys
    self.transform_to_ncs = {}
    # sorted list of transformation
    self.transform_order = []

    # length of all chains in ncs
    self.total_ncs_length = None
    # ncs + all coordinates that no transform is applied to
    self.total_length_extended_ncs = None
    # length of all chains in asu
    self.total_asu_length = None
    # flex.bool ncs_atom_selection
    self.ncs_atom_selection = None
    self.ncs_selection_str = None
    self.model_unique_chains_ids = tuple()
    self.model_order_chain_ids = []
    self.transform_to_be_used = set()

    # Use to produce new unique names for atom selections
    self.selection_names_index = [65,65]



  def populate_ncs_group_object(self,
                                ncs_refinement_params = None,
                                transform_info = None,
                                pdb_hierarchy_inp = None):
    """
    Use the iotbx.pdb.hierarchy.input, transforms info and custom
    instructions to build an ncs_group_object
    """
    # process params
    if ncs_refinement_params:
      if isinstance(ncs_refinement_params,str):
        ncs_refinement_params = parse(ncs_refinement_params)
      working_phil = master_phil.fetch(source=ncs_refinement_params).extract()
      self.ncs_refinement_groups = working_phil.ncs_refinement.ncs_group
      self.ncs_selection_str = working_phil.ncs_refinement.ncs_selection
      self.dont_apply_when_coordinates_present = \
        working_phil.ncs_refinement.dont_apply_when_coordinates_present
      self.process_phil_param()

    if pdb_hierarchy_inp:
      self.process_pdb(
        transform_info=transform_info,pdb_hierarchy_inp=pdb_hierarchy_inp)
      if not self.ncs_selection_str and not self.ncs_atom_selection:
        self.ncs_group = self.ncs_group_pdb
        s = ' or chain '.join(self.model_unique_chains_ids)
        self.ncs_selection_str = 'chain ' + s
        temp = pdb_hierarchy_inp.hierarchy.atom_selection_cache()
        self.ncs_atom_selection = temp.selection(self.ncs_selection_str)
        self.map_keys_to_selection = {k:k.split('_')[0] for k in self.ncs_group}
    # transformation application order
    self.transform_order =sorted(
      self.transform_to_ncs,key=lambda key: int(key[1:]))

    # TODO:replace collective transform on chain ID with consideration to order
    for tr_id in self.transform_order:
      for ch_id in self.model_unique_chains_ids:
        self.transform_chain_assignment.append(ch_id + '_' + tr_id[1:])

    self.ncs_copies_chains_names = self.make_chains_names(
      i_transforms = self.transform_to_be_used,
      unique_chain_names = self.model_unique_chains_ids)
    if pdb_hierarchy_inp:
      self.compute_ncs_asu_coordinates_map(pdb_hierarchy_inp=pdb_hierarchy_inp)
    if not self.ncs_refinement_groups:
      self.ncs_transform = self.ncs_transform_pdb
    for k,v in self.ncs_group.iteritems():
      key = k.split('_')[0]
      if self.number_of_ncs_copies.has_key(key):
        self.number_of_ncs_copies[key] += 1
      else:
         self.number_of_ncs_copies[key] = 1

  def compute_ncs_asu_coordinates_map(self,pdb_hierarchy_inp):
    """
    Calculates coordinates maps from ncs to asu and from asu to ncs
    """
    # in the case of regular pdb processing
    if self.apply_to_all_chains:
      if self.dont_apply_when_coordinates_present:
        self.total_ncs_length = len(pdb_hierarchy_inp.hierarchy.atoms())

    if not self.total_length_extended_ncs:
      self.total_length_extended_ncs = self.total_ncs_length

    i = self.total_length_extended_ncs
    temp = pdb_hierarchy_inp.hierarchy.atom_selection_cache()
    # TODO:  Change asu_to_ncs_map to selection type
    for k in self.transform_chain_assignment:
      key =  k.split('_')[0]
      temp_selection = temp.selection('chain ' + key)
      self.asu_to_ncs_map[key] = temp_selection.iselection()
      d = len(self.asu_to_ncs_map[key])
      self.ncs_to_asu_map[k] = [i,i + d]
      i += d

    self.total_asu_length = i
    # update ncs_atom_selection to the correct asu length
    ncs_length = self.total_length_extended_ncs
    asu_length = self.total_asu_length
    temp = flex.bool([False]*(asu_length - ncs_length))
    self.ncs_atom_selection.extend(temp)

  def process_phil_param(self):
    """ """
    # When using phil parameters, the selection string is used as a chain ID
    self.apply_to_all_chains = not self.ncs_refinement_groups
    if not self.ncs_selection_str: self.ncs_selection_str = ''
    selection_ids = set()
    serial_num = 2
    for ncsg in self.ncs_refinement_groups:
      assert ncsg.apply_to_selection
      # create ncs selection
      selection_key = ' '.join(ncsg.apply_to_selection)
      if not selection_key in selection_ids:
        self.model_order_chain_ids.append(selection_key)
      selection_ids.add(selection_key)

      # build transform objects
      self.ncs_selection_str += ' or (' + selection_key + ')'
      for tr in ncsg.transform:
        r = tr.rotation
        t = tr.translation
        cp = tr.coordinates_present
        cp = (cp is None) or cp
        assert len(r) == 9
        assert len(t) == 3
        tranform_obj = transform(
          rotation = matrix.sqr(r),
          translation = matrix.col(t),
          serial_num = serial_num,
          coordinates_present = cp)
        key = 's' + format_num_as_str(serial_num)
        serial_num += 1
        self.ncs_transform[key] = tranform_obj

    # apply all transforms that are not present to all chains
    for k,v in self.ncs_transform.iteritems():
      if not v.coordinates_present:
        self.transform_to_be_used.add(v.serial_num)
        for chain_id in self.model_unique_chains_ids:
          key = chain_id + '_' + format_num_as_str(v.serial_num)
          transform_key = 's' + format_num_as_str(v.serial_num)
          self.ncs_group[key] = transform_key
          if self.transform_to_ncs.has_key(transform_key):
            self.transform_to_ncs[transform_key].append(key)
          else:
            self.transform_to_ncs[transform_key] = [key]



  def process_pdb(self,transform_info,pdb_hierarchy_inp):
    """
    Process PDB Hierarchy object

    Remarks:
    The information on a chain in a PDB file does not have to be continuous.
    Every time the chain name changes in the pdb file, a new chain is added
    to the model, even if the chain ID already exist. so there model.
    chains() might contain several chains that have the same chain ID
    """
    model  = pdb_hierarchy_inp.hierarchy.models()[0]
    chain_ids = {x.id for x in model.chains()}
    # Collect order if chains IDs and unique IDs
    self.model_unique_chains_ids = tuple(sorted(chain_ids))
    self.model_order_chain_ids = [x.id for x in model.chains()]

    # number of chains in hierarchy (more than the actual chains in the model)
    # collect transforms
    temp_sn = 1
    sn = set()
    ti = transform_info
    for (r,t,n,cp) in zip(ti.r,ti.t,ti.serial_number,ti.coordinates_present):
      n = int(n)
      t = transform(
        rotation = r,
        translation = t,
        serial_num = n,
        coordinates_present = cp)
      if n:
        sn.add(n)
        key = 's' + format_num_as_str(n)
      else:
        key = 't' + format_num_as_str(temp_sn)
        temp_sn += 1
      assert not self.ncs_transform_pdb.has_key(key)
      self.ncs_transform_pdb[key] = t
    # replace any temp keys (keys that start with t), if there are any
    temp_sn = 1
    max_sn = max(sn) + 1
    while self.ncs_transform_pdb.has_key('t' + format_num_as_str(temp_sn)):
      key = 's' + format_num_as_str(max_sn)
      temp_key = 't' + format_num_as_str(temp_sn)
      self.ncs_transform_pdb[key] = self.ncs_transform_pdb[temp_key]
      self.ncs_transform_pdb[key].serial_num = max_sn
      self.ncs_transform_pdb.pop(temp_key)
      max_sn += 1
      temp_sn += 1
    # apply all transforms that are not present to all chains
    for k,v in self.ncs_transform_pdb.iteritems():
      if not v.coordinates_present:
        self.transform_to_be_used.add(v.serial_num)
        for chain_id in self.model_unique_chains_ids:
          key = chain_id + '_' + format_num_as_str(v.serial_num)
          transform_key = 's' + format_num_as_str(v.serial_num)
          self.ncs_group_pdb[key] = transform_key
          if self.transform_to_ncs.has_key(transform_key):
            self.transform_to_ncs[transform_key].append(key)
          else:
            self.transform_to_ncs[transform_key] = [key]


  def build_MTRIX_object(self):
    """
    Build a MTRIX object from ncs_group_object
    """
    result = iotbx.pdb._._mtrix_and_biomt_records_container()
    tr_dict = self.ncs_transform_pdb
    tr_sorted = sorted(tr_dict,key=lambda k:tr_dict[k].serial_num)
    for key in tr_sorted:
      transform = self.ncs_transform_pdb[key]
      result.add(
        r=transform.rotation,
        t=transform.translation,
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


  def make_chains_names(self,i_transforms, unique_chain_names):
    """
    Create a dictionary names for the new NCS copies
    keys: (str) chain_name + '_s' + serial_num
    values: (str) (one or two chr long)

    Chain names might repeat themselves several times in a pdb file
    We want copies of chains with the same name to still have the
    same name after similar BIOMT/MTRIX transformation

    Arguments:
    i_transform : (set) set of BIOMT/MTRIX serial number
    unique_chain_names : (tuple) a set of unique chain names or selection name

    Returns:
    new_names : a dictionary. {'A_1': 'G', 'A_2': 'H',....} map a chain
    name and a transform number to a new chain name

    >>> self.make_chains_names((1,2),('A','B'))
    {'A_1': 'C', 'A_2': 'D', 'B_1': 'E', 'B_2': 'F'}
    """
    if not i_transforms or not unique_chain_names: return {}
    # create list of character from which to assemble the list of names
    total_chains_number = len(i_transforms)*len(unique_chain_names)
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
    dictinary_key = list(set(
      [x + '_' + format_num_as_str(y) for x in unique_chain_names for y in
                                   i_transforms]))
    dictinary_key.sort()
    # create the dictionary
    zippedlists = zip(dictinary_key,dictionary_values)
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
    asu_xyz = flex.vec3_double(ncs_coordinates)

    for trans in self.transform_chain_assignment:
      s_asu,e_asu = self.ncs_to_asu_map[trans]
      ncs_selection = self.asu_to_ncs_map[trans.split('_')[0]]
      ncs_xyz = ncs_coordinates.select(ncs_selection)
      tr_key = 's' + trans.split('_')[1]
      assert self.ncs_transform.has_key(tr_key)
      tr = self.ncs_transform[tr_key]
      new_sites = tr.rotation.elems* ncs_xyz+tr.translation
      asu_xyz.extend(new_sites)
      # make sure that the new ncs copies are added in the expected order
      assert len(new_sites) == (e_asu - s_asu)
      assert len(asu_xyz) == e_asu
    if round_coordinates:
      return flex.vec3_double(asu_xyz).round(3)
    else:
      return flex.vec3_double(asu_xyz)


class transform(object):

  def __init__(self,
               rotation,
               translation,
               serial_num,
               coordinates_present):
    """
    Basic transformation properties

    Argument:
    rotation : Rotation matrix object
    translation: Translation matrix object
    serial_num : (int) Transform serial number
    coordinates_present: equals 1 when coordinates are presents in PDB file
    """
    self.rotation = rotation
    self.translation = translation
    self.serial_num = serial_num
    self.coordinates_present = False
    if coordinates_present and coordinates_present == 1:
      self.coordinates_present = True


def format_num_as_str(n):
  """  return a 3 digit string of n  """
  if n > 999 or n < 0:
    raise IOError('Input out of the range 0 - 999.')
  else:
    s1 = n//100
    s2 = (n%100)//10
    s3 = n%10
    return str(s1) + str(s2) + str(s3)
