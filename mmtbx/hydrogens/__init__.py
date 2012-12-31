from __future__ import division
import cctbx.array_family.flex # import dependency

hydrogens_master_params_str = """
    refine = individual *riding
      .type = choice
      .help = Choice for refinement: riding model or full (H is refined as \
              other atoms, useful at very high resolutions only)
      .short_caption = Hydrogen refinement model
      .expert_level=1
    optimize_scattering_contribution = True
      .type = bool
    contribute_to_f_calc = True
      .type = bool
      .help = Add H contribution to Xray (Fcalc) calculations
      .short_caption=Include hydrogens in Fcalc
      .expert_level=1
    high_resolution_limit_to_include_scattering_from_h = 1.6
      .type = float
      .short_caption = High-resolution limit to include scattering from H
      .expert_level=2
    real_space_optimize_x_h_orientation = True
      .type = bool
      .short_caption = Optimize X-H orientation in real-space
      .expert_level = 1
    xh_bond_distance_deviation_limit = 0.0
      .type = float
      .help = Idealize XH bond distances if deviation from ideal is greater \
              than xh_bond_distance_deviation_limit
      .short_caption = X-H bond distance deviation limit
      .expert_level=2
"""
