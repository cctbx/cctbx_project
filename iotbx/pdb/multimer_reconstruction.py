from iotbx import pdb
from scitbx import matrix 
import string
import math


class multimer(object):
  '''
  Reconstruction of the biological assembly multimer by applying
  BIOMT transformation, from the pdb file, to all chains.
  
  self.assembled_multimer is a pdb.hierarchy object with the multimer information
  
  The method write gatheres a PDB file containing the multimer  
  
  since chain names string length is limited to two, this process does not maintain any 
  referance to the original chain names.

  '''  
  def __init__(self,pdb_input_file_name):
    ''' (str) -> NoType
    Arguments:
    pdb_input_file_name -- the name of the pdb file we want to process. a string such as 'pdb_file_name.pdb'
    '''          
    # Read and process the pdb file
    self.pdb_input_file_name = pdb_input_file_name		# store input file name
    pdb_inp = pdb.input(file_name=pdb_input_file_name)		# read the pdb file data
    pdb_obj = pdb.hierarchy.input(file_name=pdb_input_file_name)# read the pdb file hierarchy data
    pdb_obj_new = pdb_obj.hierarchy.deep_copy()			# create a copy to be modified
    self.assembled_multimer = pdb_obj_new
    
    # take care of potenitial input issues
    
    
    # Read BIOMT info
    BIOMT_info = pdb_inp.process_BIOMT_records()
    # convert BIOMT object to a list of matrices
    BIOMT = self._convert_lists_to_matrices(BIOMT_info)
    BIOMT_transform_number = len(BIOMT_info)
    # number of chains
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
      new_chains_names = self._chains_names(BIOMT_transform_number,nChains, unique_chain_names)
      for iChain in range(nChains):	#instead of: for chain in model.chains():
        # iterating over the BIOMT transform strating from 1 (not 0)
        # since BIOMT_info[0] is a unit transformation
        for i_transform in range(1,BIOMT_transform_number):
          new_chain = model.chains()[iChain].detached_copy()
          new_chain.id = new_chains_names[new_chain.id + str(i_transform)]
          for res_group in new_chain.residue_groups():
            for ag in res_group.atom_groups():
              for atom in ag.atoms():
                # apply transform               
                xyz = matrix.col(atom.xyz)	# get current x,y,z coordinate
                xyz = BIOMT[i_transform][0]*xyz # apply rotation
                xyz +=  BIOMT[i_transform][1]	# apply translation
                atom.xyz = xyz.elems            # replace the new chain x,y,z coordinates
          # add a new chain to current model
          model.append_chain(new_chain)
    
  def _convert_lists_to_matrices(self, BIOMT_info):  
    ''' (object contianing lists) -> list of scitbx matrices and vectors
    
    Convert the rotation and translation information from lists to matrices and vectors
    
    Argumnets:
    BIOMT_info -- a object of lists containing the BIOMT transformation information obtained from a pdb file
    The transformation information in BIOMT_info looks like
    BIOMT_info[0].values = [[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 0.0]]
    
    Returns:
    BIOMT_matrices_list -- a list of matrices and vectors
    [[rotation matrix 1,translation vector 1], [rotation matrix 2,translation vector 2],...]
    
    >>> BIOMT_matrices_list[0][0]
    matrix([[ 1.,  0.,  0.],
        [ 0.,  1.,  0.],
        [ 0.,  0.,  1.]])
    >>> len(BIOMT_matrices_list) == len(BIOMT_info)
    True
    '''
    
    BIOMT_matrices_list = [[matrix.sqr(x.values[0]), matrix.col(x.values[1])] for x in BIOMT_info]
    return BIOMT_matrices_list
      
  def _chains_names(self, BIOMT_transform_numbers,nChains, unique_chain_names):
    ''' (int, int, set) -> dictionary
    
    Create a dictionary 
    keys: a string made of chain_name + str(BIOMT_transform_number)
    values: two letters and digits string combination
    
    The total number of new chains can be large (order of hundereds)
    Chain names might repeat themselves several times in a pdb file
    We want copies of chains with the same name to stil have the same name after similar BIOMT transformation
    
    Arguments:
    BIOMT_transform_numbers -- an integer. the number of BIOMT transformation operations in the pdb object
    nChains -- an integer. the number of chains as interpreted by pdb.hierarchy.input
    unique_chain_names -- a set. a set of unique chain names
    
    Returns:
    new_names -- a dictionary. {'A1': 'aa', 'A2': 'gq',....} map a chain name and 
    a BIOMT transform number to a new chain name
    
    >>> self._chains_names(3,4,{'A','B'})
    {'A1': 'aa', 'A3': 'ac', 'A2': 'ab', 'B1': 'ba', 'B2': 'bb', 'B3': 'bc'}
    '''
    # create list of character from which to assemble the list of names
    total_chains_number = BIOMT_transform_numbers*len(unique_chain_names)
    chr_number = int(math.sqrt(total_chains_number)) + 1	# the number of charater needed to produce new names
    chr_list = list(string.ascii_letters) + list(string.digits)	# build character list
    chr_list = chr_list[:chr_number]				# take only as many characters as needed
    dictionary_values = set([ x+y for x in chr_list for y in chr_list])
    dictinary_key = set([x+str(y) for x in unique_chain_names for y in range(1,BIOMT_transform_numbers+1)])
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
        for res_group in chain.residue_groups():
          for ag in res_group.atom_groups():
            for atom in ag.atoms():
              xyz.append(list(atom.xyz))	# get x,y,z coordinates
    return xyz

        
  def __str__(self):
    ''' (multimer) -> string
    
    Retruns a string with information about this multimer
    
    '''
    return 'the multimer attribute self.assembled_multimer is a pdb.hierarchy object that contains: \n' + \
           'a model in assembled_multimer.models()  \n' + \
           'a chain in model.chains()  \n' + \
           'a residue_group in chain.residue_groups()  \n' + \
           'an atom_group in residue_group.atom_groups() \n' + \
           'an atom in atom_group.atoms() \n' + \
           'atom.xyz is a tuple with the atom coordinates: (x,y,z)'
    
    
  def write(self,pdb_output_file_name):
    ''' (string) -> text file
    Writes the modified protein, with the added chains, obtained by the BIOMT
    reconstruction, to a text file in a pdb format. 
    self.assembled_multimer is the modified pdb object with the added chains
    
    Argumets:
    pdb_output_file_name -- string. 'name.pdb'
    
    >>> v = multimer('name.pdb')
    >>> v.write('new_name.pdb')
    should write a file 'new_name.pdb' to the current directory
    >>> v.write(v.pdb_input_file_name)
    should write a file 'copy_name.pdb' to the current directory
    ''' 
    # Aviod writing over the original file
    if self.pdb_input_file_name == pdb_output_file_name:
      # if file name of output is the same as the input, add 'copy_' in front of the name
      pdb_output_file_name = 'copy_' + pdb_output_file_name
        
    # using the pdb hierarchy pdb file writing method
    self.assembled_multimer.write_pdb_file(file_name = pdb_output_file_name)
    
    

