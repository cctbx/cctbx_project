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
          pdb_string=None):
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
              ncs_group_selection {
                ncs_group (multiple)
                {
                  master_ncs_selection = ''
                  selection_copy = ''   (multiple)
                }
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
      pdb_string=pdb_string)
    return ncs_group_obj

def apply_transforms(ncs_coordinates,
                     ncs_restraints_group_list,
                     total_asu_length,
                     round_coordinates = True):
  """
  Apply transformation to ncs_coordinates,
  and round the results if round_coordinates is True

  Argument:
  ncs_coordinates: (flex.vec3) master ncs coordinates
  ncs_restraints_group_list: list of ncs_restraint_group objects
  total_asu_length: (int) Complete ASU length

  Returns:
  complete asymmetric or the biological unit
  """
  s = flex.size_t_range(len(ncs_coordinates))
  asu_xyz = flex.vec3_double([(0,0,0)]*total_asu_length)
  asu_xyz.set_selected(s,ncs_coordinates)

  for nrg in ncs_restraints_group_list:
    master_ncs_selection = nrg.master_ncs_iselection
    asu_selection = nrg.ncs_copy_iselection
    ncs_xyz = ncs_coordinates.select(master_ncs_selection)
    new_sites = nrg.r.elems* ncs_xyz+nrg.t
    asu_xyz.set_selected(asu_selection,new_sites)
  if round_coordinates:
    return flex.vec3_double(asu_xyz).round(3)
  else:
    return flex.vec3_double(asu_xyz)
