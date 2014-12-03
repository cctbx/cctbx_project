from __future__ import division
import sys
import ncs_preprocess


def input(pdb_hierarchy_inp=None,
          pdb_inp=None,
          hierarchy=None,
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
          use_minimal_master_ncs=True,
          max_rmsd=2.0,
          write_messages=False,
          log=None,
          process_similar_chains=True,
          min_percent=0.85,
          chain_similarity_limit=0.95,
          min_contig_length=10,
          check_atom_order=False,
          allow_different_size_res=True,
          exclude_misaligned_residues=False,
          max_dist_diff=4.0):
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
      max_rmsd (float): limit of rms difference between chains to be considered
        as copies
      write_messages (bool): When True, write messages to log
        nearly the same length (but not exactly the same) and are NCS related.
        Raise error if NCS relations are not found
      process_similar_chains (bool): When True, process chains that are close
       in length without raising errors
      min_percent (float): Threshold for similarity between chains
        similarity define as:
        (number of matching res) / (number of res in longer chain)
      chain_similarity_limit (float): min similarity between matching chains
      min_contig_length (int): minimum length of matching chain segments
      check_atom_order (bool): check atom order in matching residues.
        When False, matching residues with different number of atoms will be
        excluded from matching set
      allow_different_size_res (bool): keep matching residue with different
        number of atoms
      exclude_misaligned_residues (bool): check and exclude individual residues
        alignment quality
      max_dist_diff (float): max allow distance difference between pairs of matching
        atoms of two residues
    """
    if not log: log = sys.stdout
    ncs_group_obj = ncs_preprocess.ncs_group_object()
    ncs_group_obj.preprocess_ncs_obj(
      pdb_hierarchy_inp=pdb_hierarchy_inp,
      pdb_inp=pdb_inp,
      hierarchy=hierarchy,
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
      use_minimal_master_ncs=use_minimal_master_ncs,
      max_rmsd=max_rmsd,
      write_messages=write_messages,
      log=log,
      process_similar_chains=process_similar_chains,
      min_percent=min_percent,
      min_contig_length=min_contig_length,
      check_atom_order=check_atom_order,
      allow_different_size_res=allow_different_size_res,
      exclude_misaligned_residues=exclude_misaligned_residues,
      max_dist_diff=max_dist_diff,
      chain_similarity_limit=chain_similarity_limit)
    return ncs_group_obj

