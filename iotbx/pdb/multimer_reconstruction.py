from __future__ import division
from iotbx import crystal_symmetry_from_any
from libtbx.utils import Sorry
from scitbx import matrix
from iotbx import pdb
import string
import math
import os

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
  def __init__(self,pdb_input_file_name,reconstruction_type,error_handle=True,eps=1e-4):
    ''' (str) -> NoType
    Arguments:
    pdb_input_file_name -- the name of the pdb file we want to process. a string such as 'pdb_file_name.pdb'
    reconstruction_type -- 'ba' or 'cau'
                           'ba': biological assembly
                           'cau': crystallographic asymmetric unit
    error_handle -- True: will stop execution on improper retation matrices
                    False: will continue execution but will replace the values in the
                           rotation matrix with [0,0,0,0,0,0,0,0,0]
    eps -- Rounding accuracy for avoiding numerical issue when when testing proper rotation

    @author: Youval Dar (2013)
    '''
    # Read and process the pdb file
    self.pdb_input_file_name = pdb_input_file_name
    pdb_inp = pdb.input(file_name=pdb_input_file_name)
    pdb_obj = pdb.hierarchy.input(file_name=pdb_input_file_name)
    pdb_obj_new = pdb_obj.hierarchy.deep_copy()
    self.assembled_multimer = pdb_obj_new

    # Read the relevant transformation matrices
    if reconstruction_type == 'ba':
      # Read BIOMT info
      TRASFORM_info = pdb_inp.process_BIOMT_records(
        error_handle=error_handle,eps=eps)
      self.transform_type = 'biological_assembly'
    elif reconstruction_type == 'cau':
      # Read MTRIX info
      TRASFORM_info = pdb_inp.process_mtrix_records(
        error_handle=error_handle,eps=eps)
      self.transform_type = 'crystall_asymmetric_unit'
    else:
      raise Sorry('Worg reconstruction type is given \n' + \
                  'Reconstruction type can be: \n' + \
                  "'ba': biological assembly \n" + \
                  "'cau': crystallographic asymmetric unit \n")
    # Collect only transformations which coordinates are not present
    i_transforms = []
    for i,coordinates_present in enumerate(TRASFORM_info.coordinates_present):
      if not coordinates_present: i_transforms.append(i)
    self.number_of_transforms = len(i_transforms)
    if (self.number_of_transforms>0):
      # number of chains in hierachy (more than the actual chains in the model)
      chains_number = pdb_obj.hierarchy.overall_counts().n_chains
      # apply the transformation
      for model in pdb_obj_new.models():
        # The information on a chain in a PDB file does not have to be continuous.
        # Every time the chain name changes in the pdb file, a new chain is added to the model,
        # even if the chain ID already exist.
        # so there model.chains() might containe several chains that have the same chain ID
        # collect the unique chain names
        unique_chain_names = {x.id for x in model.chains()}
        nChains = len(model.chains())
        # get a dictionary for new chains naming
        new_chains_names = self._chains_names(i_transforms,\
          nChains,unique_chain_names)
        for chain in model.chains():
          # iterating over the TRASFORM transforms that are not present
          for i_transform in i_transforms:
            new_chain = chain.detached_copy()
            new_chain.id = new_chains_names[new_chain.id + str(i_transform)]
            new_sites = TRASFORM_info.r[i_transform].elems*\
              new_chain.atoms().extract_xyz()+TRASFORM_info.t[i_transform]
            new_chain.atoms().set_xyz(new_sites)
            # add a new chain to current model
            model.append_chain(new_chain)

  def _chains_names(self, i_transforms,nChains, unique_chain_names):
    ''' (int, int, set) -> dictionary

    Create a dictionary
    keys: a string made of chain_name + str(TRASFORM_transform_number)
    values: two letters and digits string combination

    The total number of new chains can be large (order of hundereds)
    Chain names might repeat themselves several times in a pdb file
    We want copies of chains with the same name to stil have the same name after
    similar BIOMT/MTRIX transformation

    Arguments:
    TRASFORM_transform_numbers -- an integer. the number of BIOMT/MTRIX transformation operations in the pdb object
    nChains -- an integer. the number of chains as interpreted by pdb.hierarchy.input
    unique_chain_names -- a set. a set of unique chain names

    Returns:
    new_names -- a dictionary. {'A1': 'aa', 'A2': 'gq',....} map a chain name and
    a TRASFORM transform number to a new chain name

    >>> self._chains_names(3,4,{'A','B'})
    {'A1': 'aa', 'A3': 'ac', 'A2': 'ab', 'B1': 'ba', 'B2': 'bb', 'B3': 'bc'}
    '''
    # create list of character from which to assemble the list of names
    total_chains_number = len(i_transforms)*len(unique_chain_names)
    chr_number = int(math.sqrt(total_chains_number)) + 1        # the number of charater needed to produce new names
    chr_list = list(string.ascii_letters) + list(string.digits) # build character list
    chr_list = chr_list[:chr_number]                            # take only as many characters as needed
    dictionary_values = set([ x+y for x in chr_list for y in chr_list])
    dictinary_key = set([x+str(y) for x in unique_chain_names for y in i_transforms])
    # create the dictionary
    new_names_dictionary ={x:y for (x,y) in zip(dictinary_key,dictionary_values)}
    return new_names_dictionary

  def get_xyz(self):
    '''(multimer) -> list of arrays

    Returns:
    xyz -- a list of x,y,z coordinates for each atoms in the resconstructed multimer
    '''
    xyz = []
    for model in self.assembled_multimer.models():
      for chain in model.chains():
        xyz.extend(list(chain.atoms().extract_xyz()))
    return xyz

  def write(self,pdb_output_file_name=''):
    ''' (string) -> text file
    Writes the modified protein, with the added chains, obtained by the BIOMT/MTRIX
    reconstruction, to a text file in a pdb format.
    self.assembled_multimer is the modified pdb object with the added chains

    Argumets:
    pdb_output_file_name -- string. 'name.pdb'
    if no pdn_output_file_name is given pdb_output_file_name=pdb_input_file_name

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
    crystal_symmetry = crystal_symmetry_from_any.extract_from(
      self.pdb_input_file_name)
    # using the pdb hierarchy pdb file writing method
    self.assembled_multimer.write_pdb_file(file_name=self.pdb_output_file_name,
      crystal_symmetry=crystal_symmetry)
