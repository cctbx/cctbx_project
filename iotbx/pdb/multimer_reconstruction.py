from __future__ import division
from iotbx.ncs.ncs_preprocess import insure_identity_is_in_transform_info
from  iotbx.ncs.ncs_preprocess import ncs_group_object
from  iotbx.ncs.ncs_preprocess import ncs_only
from iotbx import crystal_symmetry_from_any
from libtbx.utils import Sorry
# import iotbx.pdb.hierarchy
from iotbx import pdb
import os


class multimer(object):
  '''
  Reconstruction of either the biological assembly or the crystallographic
  asymmetric unit

  Reconstruction of the biological assembly multimer by applying
  BIOMT transformation, from the pdb file, to all chains.

  Reconstruction of the crystallographic asymmetric unit by applying
  MTRIX transformations, from the pdb file, to all chains.

  self.assembled_multimer is a pdb.hierarchy object with the multimer
  information

  The method write generates a PDB file containing the multimer

  since chain names string length is limited to two, this process does
  not maintain any referance to the original chain names.

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
    # Todo : add ability to reconstruct ASU for multiple groups
    assert file_name or pdb_str
    # Read and process the pdb file
    if file_name:
      self.pdb_input_file_name = file_name
      pdb_obj = pdb.hierarchy.input(file_name=file_name)
      pdb_inp = pdb.input(file_name=file_name)
    else:
      self.pdb_input_file_name = 'complete_reconstructed_unit.pdb'
      pdb_obj = pdb.hierarchy.input(pdb_string=pdb_str)
      pdb_inp = pdb.input(lines=pdb_str, source_info=None)
    if reconstruction_type == 'ba':
      transform_info = pdb_inp.process_BIOMT_records(
        error_handle=error_handle,
        eps=eps)
      self.transform_type = 'biological_assembly'
    elif reconstruction_type == 'cau':
      transform_info = pdb_inp.process_mtrix_records(
        error_handle=error_handle,
        eps=eps)
      if transform_info.as_pdb_string() == '' or (not ncs_only(transform_info)):
        transform_info = None
      else:
        transform_info = insure_identity_is_in_transform_info(transform_info)
      self.transform_type = 'crystall_asymmetric_unit'
    else:
      raise Sorry('Sorry, wrong reconstruction type is given \n' + \
                  'Reconstruction type can be: \n' + \
                  "'ba': biological assembly \n" + \
                  "'cau': crystallographic asymmetric unit \n")
    if len(pdb_obj.hierarchy.models()) > 1:
      raise Sorry('Sorry, this feature currently supports on single models ' +
                  'hierarchies')

    self.transforms_obj = ncs_group_object()
    # Read the relevant transformation matrices
    self.transforms_obj.build_ncs_obj_from_pdb_ncs(
      transform_info=transform_info,
      pdb_hierarchy_inp=pdb_obj)

    # Calculate ASU (if there are any transforms to apply)
    self.number_of_transforms = len(self.transforms_obj.transform_to_be_used)
    self.assembled_multimer = self.transforms_obj.build_asu_hierarchy(
      pdb_hierarchy=pdb_obj.hierarchy,
      round_coordinates=round_coordinates)
    annot = pdb_inp.extract_secondary_structure()
    self.new_annotation = None
    if annot is not None:
      annot.multiply_to_asu(
          ncs_copies_chain_names=self.transforms_obj.ncs_copies_chains_names,
          n_copies=self.number_of_transforms)
      self.new_annotation = annot

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
    Writes the modified protein, with the added chains, obtained by the
    BIOMT/MTRIX reconstruction, to a text file in a pdb format.
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
    # Avoid writing over the original file
    if pdb_output_file_name == input_file_name:
      # if file name of output is the same as the input, add 'copy_' in front of the name
      self.pdb_output_file_name = self.transform_type + '_' + input_file_name
    else:
      self.pdb_output_file_name = pdb_output_file_name
    # we need to add crystal symmetry to the new file since it is
    # sometimes needed when calculating the R-work factor (r_factor_calc.py)
    if not crystal_symmetry and os.path.isfile(self.pdb_input_file_name):
      crystal_symmetry = crystal_symmetry_from_any.extract_from(
        self.pdb_input_file_name)
    # using the function to write whole pdb file
    pdb.write_whole_pdb_file(
        file_name=self.pdb_output_file_name,
        pdb_hierarchy=self.assembled_multimer,
        crystal_symmetry=crystal_symmetry,
        ss_annotation=self.new_annotation)

  def get_ncs_restraints_group_list(self,raise_sorry=True):
    get_nrgl = self.transforms_obj.get_ncs_restraints_group_list
    return get_nrgl(raise_sorry=raise_sorry)

  def total_asu_length(self):
    return self.transforms_obj.total_asu_length
