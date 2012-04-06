import iotbx.phil

tw_customizations_str = """\
weight_selection_criteria
  .multiple=True
{
  d_min = 0.0
    .type=float
  d_max = 1.5
    .type=float
  bonds_rmsd = 0.025
    .type=float
  angles_rmsd = 3.0
    .type=float
  r_free_minus_r_work = 4
    .type=float
  r_free_range_width = 0
    .type=float
  mean_diff_b_iso_bonded_fraction = 0.1
    .type=float
  min_diff_b_iso_bonded = 10
    .type=float
}
weight_selection_criteria
  .multiple=True
{
  d_min = 1.5
    .type=float
  d_max = 2.0
    .type=float
  bonds_rmsd = 0.02
    .type=float
  angles_rmsd = 2.5
    .type=float
  r_free_minus_r_work = 5
    .type=float
  r_free_range_width = 0
    .type=float
  mean_diff_b_iso_bonded_fraction = 0.1
    .type=float
  min_diff_b_iso_bonded = 10
    .type=float
}
weight_selection_criteria
  .multiple=True
{
  d_min = 2.0
    .type=float
  d_max = 2.5
    .type=float
  bonds_rmsd = 0.015
    .type=float
  angles_rmsd = 2.0
    .type=float
  r_free_minus_r_work = 6
    .type=float
  r_free_range_width = 1
    .type=float
  mean_diff_b_iso_bonded_fraction = 0.1
    .type=float
  min_diff_b_iso_bonded = 10
    .type=float
}
weight_selection_criteria
  .multiple=True
{
  d_min = 2.5
    .type=float
  d_max = 3.5
    .type=float
  bonds_rmsd = 0.015
    .type=float
  angles_rmsd = 2.0
    .type=float
  r_free_minus_r_work = 6
    .type=float
  r_free_range_width = 1.5
    .type=float
  mean_diff_b_iso_bonded_fraction = 0.1
    .type=float
  min_diff_b_iso_bonded = 10
    .type=float
}
weight_selection_criteria
  .multiple=True
{
  d_min = 3.5
    .type=float
  d_max = 1.e+6
    .type=float
  bonds_rmsd = 0.015
    .type=float
  angles_rmsd = 2.0
    .type=float
  r_free_minus_r_work = 7
    .type=float
  r_free_range_width = 1.5
    .type=float
  mean_diff_b_iso_bonded_fraction = 0.2
    .type=float
  min_diff_b_iso_bonded = 20
    .type=float
}
"""
tw_customizations_params = iotbx.phil.parse(tw_customizations_str)
