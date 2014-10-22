from __future__ import division
import iotbx.ncs.ncs_preprocess
import sys


def input(pdb_hierarchy_inp=None,
          pdb_inp=None,
          transform_info=None,
          rotations = None,
          translations = None,
          ncs_selection_params = None,
          ncs_phil_groups = None,
          file_name=None,
          file_path='',
          spec_file_str='',
          spec_source_info='',
          cif_string = '',
          quiet=True,
          spec_ncs_groups=None,
          pdb_string=None,
          use_cctbx_find_ncs_tools=True,
          use_simple_ncs_from_pdb=False,
          use_minimal_master_ncs=True,
          rms_eps=0.02,
          write_messages=False,
          log=sys.stdout,
          process_similar_chains=True,
          min_fraction_domain=0.2,
          min_contig_length=10,
          check_atom_order=False):
    """
    Select method to build ncs_group_object

    order of implementation:
    1) rotations,translations
    2) transform_info
    3) ncs_selection_params
    4) ncs_phil_groups
    5) spec file
    6) mmcif file
    7) iotbx.pdb.hierarchy.input object

    Args:
    -----
      pdb_hierarchy_inp: iotbx.pdb.hierarchy.input
      transform_info: object containing MTRIX or BIOMT transformation info
      rotations: matrix.sqr 3x3 object
      translations: matrix.col 3x1 object
      ncs_selection_params: Phil parameters
        Phil structure
           ncs_group (multiple)
           {
             master_selection = ''
             copy_selection = ''   (multiple)
           }
      ncs_phil_groups: a list of ncs_groups_container object, containing
        master NCS selection and a list of NCS copies selection
      file_name: (str) .ncs_spec or .mmcif  or .pdb file name
      file_path: (str)
      spec_file_str: (str) spec format data
      spec_source_info:
      quiet: (bool) When True -> quiet output when processing files
      spec_ncs_groups: ncs_groups object as produced by simple_ncs_from_pdb
      cif_string: (str) string of cif type data
      use_cctbx_find_ncs_tools: (bool) Enable using of chain base NCS search
      use_simple_ncs_from_pdb: (bool) Enable using use_simple_ncs_from_pdb
      use_minimal_master_ncs: (bool) use maximal or minimal common chains
        in master ncs groups
      rms_eps (float): limit of rms difference between chains to be considered
        as copies
      write_messages (bool): When True, write messages to log
        nearly the same length (but not exactly the same) and are NCS related.
        Raise error if NCS relations are not found
      process_similar_chains (bool): When True, process chains that are close
       in length without raising errors
      min_fraction_domain (float): Threshold for similarity between chains
        similarity define as:
        (number of matching res) / (number of res in longer chain)
      min_contig_length (int): minimum length of matching chain segments
      check_atom_order (bool): check atom order in matching residues.
        When False, matching residues with different number of atoms will be
        excluded from matching set
    """
    ncs_group_obj = iotbx.ncs.ncs_preprocess.ncs_group_object()
    ncs_group_obj.preprocess_ncs_obj(
      pdb_hierarchy_inp=pdb_hierarchy_inp,
      pdb_inp=pdb_inp,
      transform_info=transform_info,
      rotations=rotations,
      translations=translations,
      ncs_selection_params=ncs_selection_params,
      ncs_phil_groups=ncs_phil_groups,
      file_name=file_name,
      file_path=file_path,
      spec_file_str=spec_file_str,
      spec_source_info=spec_source_info,
      cif_string=cif_string,
      quiet=quiet,
      spec_ncs_groups=spec_ncs_groups,
      pdb_string=pdb_string,
      use_cctbx_find_ncs_tools=use_cctbx_find_ncs_tools,
      use_simple_ncs_from_pdb=use_simple_ncs_from_pdb,
      use_minimal_master_ncs=use_minimal_master_ncs,
      rms_eps=rms_eps,
      write_messages=write_messages,
      log=log,
      process_similar_chains=process_similar_chains,
      min_fraction_domain=min_fraction_domain,
      min_contig_length=min_contig_length,
      check_atom_order=check_atom_order)
    return ncs_group_obj

