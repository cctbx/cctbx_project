from __future__ import absolute_import, division, print_function

import iotbx.phil
import mmtbx.geometry_restraints.torsion_restraints.torsion_ncs

global_ncs_params = iotbx.phil.parse("""\
  ncs
    .caption = Parameter controls for NCS refinement (restrained, constrained)
    .short_caption = Global NCS
    .style = scrolled menu_item auto_align box
  {
    type = *torsion cartesian constraints
      .type = choice(multi=False)
      .short_caption = NCS type
      .caption = torsion-angle global constraints
      .style = bold noauto
    constraints
      .caption = Parameter controls for strict NCS refinement
      .expert_level=2
    {
      refine_operators = True
        .type = bool
        .expert_level=2
      apply_to_coordinates = True
        .type = bool
        .expert_level=2
      apply_to_adp = True
        .type = bool
        .expert_level=2
    }
    coordinate_sigma=0.05
      .type = float
      .expert_level=1
    restrain_b_factors = False
      .type = bool
      .short_caption = Restrain NCS-related B-factors
      .help = If enabled, b-factors will be restrained for NCS-related atoms. \
        Otherwise, atomic b-factors will be refined independently, and \
        b_factor_weight will be set to 0.0
      .style = bold
    b_factor_weight=10
      .type=float
      .short_caption = Weight on NCS-related B-factor restraints
    excessive_distance_limit = 1.5
      .type = float
      .style = bold
    special_position_warnings_only = False
      .type = bool
      .expert_level=2
    torsion
      .style = noauto menu_item auto_align box
      .short_caption = Torsion-angle NCS
    {
      include scope mmtbx.geometry_restraints.torsion_restraints.torsion_ncs.torsion_ncs_params
    }
  }""",process_includes=True)
