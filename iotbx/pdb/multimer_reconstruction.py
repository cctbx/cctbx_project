from __future__ import division
from iotbx import pdb
from scitbx import matrix
from scitbx.array_family import flex
from libtbx.utils import Sorry
import string
import math

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

  @author: Youval Dar (LBL)
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
    self.pdb_input_file_name = pdb_input_file_name              # store input file name
    pdb_inp = pdb.input(file_name=pdb_input_file_name)          # read the pdb file data
    pdb_obj = pdb.hierarchy.input(file_name=pdb_input_file_name)# read the pdb file hierarchy data
    pdb_obj_new = pdb_obj.hierarchy.deep_copy()                 # create a copy to be modified
    self.assembled_multimer = pdb_obj_new

    # Read the relevant transformation matrices
    if reconstruction_type == 'ba':
      # Read BIOMT info
      TRASFORM_info = pdb_inp.process_BIOMT_records(error_handle=error_handle,eps=eps)
      self.transform_type = 'biological_assembly'
    elif reconstruction_type == 'cau':
      # Read MTRIX info
      TRASFORM_info = pdb_inp.process_mtrix_records()
      self.transform_type = 'crystall_asymmetric_unit'
    else:
      raise Sorry('Worg reconstruction type is given \n' + \
                  'Reconstruction type can be: \n' + \
                  "'ba': biological assembly \n" + \
                  "'cau': crystallographic asymmetric unit \n")

    # convert TRASFORM object to a list of matrices
    TRASFORM = self._convert_lists_to_matrices(TRASFORM_info)
    TRASFORM_transform_number = len(TRASFORM)
    self.number_of_transforms = TRASFORM_transform_number
    if TRASFORM_transform_number > 0:
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
        new_chains_names = self._chains_names(TRASFORM_transform_number,nChains, unique_chain_names)
        for chain in model.chains():
          # iterating over the TRASFORM transform strating from 1 (not 0)
          # since TRASFORM_info[0] is a unit transformation
          for i_transform in range(TRASFORM_transform_number):
            new_chain = chain.detached_copy()
            new_chain.id = new_chains_names[new_chain.id + str(i_transform+1)]
            sites = new_chain.atoms().extract_xyz()
            # calculating new sites
            new_sites = TRASFORM[i_transform].r.elems*sites + tuple(TRASFORM[i_transform].t)
            new_chain.atoms().set_xyz(new_sites)
            # add a new chain to current model
            model.append_chain(new_chain)

  def _convert_lists_to_matrices(self, TRASFORM_info):
    ''' (object contianing lists) -> list of scitbx matrices and vectors

    Convert the rotation and translation information from lists to rotation and
    translation components of matrix.rt object.

    Argumnets:
    TRASFORM_info -- a object of lists containing the BIOMT/MTRIX transformation information obtained
    from a pdb file.
    The transformation information in TRASFORM_info looks like
    TRASFORM_info[0].values = [[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 0.0]]

    Returns:
    TRASFORM_matrices_list -- a list of matrix.rt objects

    >>> TRASFORM_matrices_list.r[0]
    matrix.rec(elems=(1,0,0,0,1,0,0,0,1), n=(3,3))

    >>> TRASFORM_matrices_list.r[0]
    matrix.rec(elems=(0,0,0), n=(3,1))

    >>> len(TRASFORM_matrices_list) == len(TRASFORM_info)
    True
    '''
    # TRASFORM_info.coordinates_present=True is transformation already included in pdb file
    TRASFORM_matrices_list = [matrix.rt(x.values) for x in TRASFORM_info if not x.coordinates_present]
    return TRASFORM_matrices_list

  def _chains_names(self, TRASFORM_transform_numbers,nChains, unique_chain_names):
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
    total_chains_number = TRASFORM_transform_numbers*len(unique_chain_names)
    chr_number = int(math.sqrt(total_chains_number)) + 1        # the number of charater needed to produce new names
    chr_list = list(string.ascii_letters) + list(string.digits) # build character list
    chr_list = chr_list[:chr_number]                            # take only as many characters as needed
    dictionary_values = set([ x+y for x in chr_list for y in chr_list])
    dictinary_key = set([x+str(y) for x in unique_chain_names for y in range(1,TRASFORM_transform_numbers+1)])
    # create the dictionary
    new_names_dictionary  = {x:y for (x,y) in zip(dictinary_key,dictionary_values)}
    return new_names_dictionary

  def __eq__(self, other):
    ''' (multimer, multimer) -> bool

    return true iff all the coordinates of all atoms are the same
    '''
    x = self.get_xyz()
    y = other.get_xyz()
    return x == y


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


  def __str__(self):
    ''' (multimer) -> string

    Retruns a string with information about this multimer

    '''
    return 'the multimer attribute self.assembled_multimer is a pdb.hierarchy object that contains: \n' + \
           'a model in assembled_multimer.models()  \n' + \
           '   a chain in model.chains()  \n' + \
           '      the sites which are the coordinates of the atoms \n' + \
           '      where  sites = chain.atoms().extract_xyz() ' + \
           '      The atoms xyz coordinates are tuples (x,y,z) \n' + \
           ' \n' + \
           'self.get_xyz() will return a list of x,y,z coordinates for each atoms in the resconstructed multimer \n'


  def write(self,*pdb_output_file_name):
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
    if len(pdb_output_file_name) == 0:
      pdb_output_file_name = self.pdb_input_file_name
    # Aviod writing over the original file
    if self.pdb_input_file_name == pdb_output_file_name:
      # if file name of output is the same as the input, add 'copy_' in front of the name
      self.pdb_output_file_name = self.transform_type + '_' +pdb_output_file_name
    else:
      self.pdb_output_file_name = pdb_output_file_name

    # using the pdb hierarchy pdb file writing method
    self.assembled_multimer.write_pdb_file(file_name=self.pdb_output_file_name)
