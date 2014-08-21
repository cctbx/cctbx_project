from __future__ import division
import iotbx.ncs.ncs_preprocess


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
          use_simple_ncs_from_pdb=True,
          use_minimal_master_ncs=True):
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
             master_ncs_selection = ''
             selection_copy = ''   (multiple)
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
      use_minimal_master_ncs=use_minimal_master_ncs)
    return ncs_group_obj

