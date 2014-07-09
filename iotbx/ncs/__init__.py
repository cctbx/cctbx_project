from __future__ import division
import iotbx.ncs.ncs_preprocess
from scitbx.array_family import flex


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
          use_simple_ncs_from_pdb=False):
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

    :param pdb_hierarchy_inp: iotbx.pdb.hierarchy.input
    :param transform_info: object containing MTRIX or BIOMT transformation info
    :param rotations: matrix.sqr 3x3 object
    :param translations: matrix.col 3x1 object
    :param ncs_selection_params: Phil parameters
           Phil structure
              ncs_group (multiple)
              {
                master_ncs_selection = ''
                selection_copy = ''   (multiple)
              }
    :param ncs_phil_groups: a list of ncs_groups_container object, containing
           master NCS selection and a list of NCS copies selection
    :param file_name: (str) .ncs_spec or .mmcif  or .pdb file name
    :param file_path: (str)
    :param spec_file_str: (str) spec format data
    :param spec_source_info:
    :param quiet: (bool) When True -> quiet output when processing files
    :param spec_ncs_groups: ncs_groups object as produced by simple_ncs_from_pdb
    :param cif_string: (str) string of cif type data
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
      use_simple_ncs_from_pdb=use_simple_ncs_from_pdb)
    return ncs_group_obj

