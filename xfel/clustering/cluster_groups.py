""" Utilitites for dealing with lists of clusters. """
from __future__ import absolute_import, division, print_function
__author__ = 'zeldin'

def unit_cell_info(sub_clusters):
  """
  Print unit cell information for a list of clusters.

  :param sub_clusters: a list of cluster objects
  :return: a string containing median unit cells, standard deviations and
   point group composition of each cluster.
  """
  from libtbx.utils import plural_s
  # 3. print out some information that is useful.
  out_str = "\n\n{:<16} {:<8} {:<13} {:<13} {:<13} {:<12} {:<12} {:<12}{:<8}\n".format(
    "Cluster_id",
    "N_xtals",
    "Med_a", "Med_b", "Med_c",
    "Med_alpha", "Med_beta", "Med_gamma","Delta(deg)")
  singletons = []
  for cluster in sub_clusters:
    if len(cluster.members) != 1:
      # New approach, takes niggli setting of the cluster median and converts
      # back to reference setting for cluster report. Fixes cctbx#97.
      from cctbx import crystal
      from cctbx.uctbx import unit_cell
      from cctbx.sgtbx.lattice_symmetry import metric_subgroups

      input_symmetry = crystal.symmetry(
          unit_cell=unit_cell(cluster.medians[0:6]),
          space_group_symbol="P 1")
      groups = metric_subgroups(input_symmetry, 3.00,
        enforce_max_delta_for_generated_two_folds=True)
      group = groups.result_groups[0]
      # suppress stdout output for now
      from six.moves import StringIO
      SS = StringIO()
      import sys
      sys.stdout = SS
      group['best_subsym'].space_group_info().show_summary()
      sys.stdout=sys.__stdout__
      print("                       Unit cell:", group['best_subsym'].unit_cell())
      uc_params_conv = group['best_subsym'].unit_cell().parameters()

      sorted_pg_comp = sorted(cluster.pg_composition.items(),
                              key=lambda x: -1 * x[1])
      pg_strings = ["{} in {}".format(pg[1], pg[0])
                    for pg in sorted_pg_comp]
      point_group_string = ", ".join(pg_strings) + "."
      out_str += point_group_string
      out_str += ("\n{:<16} {:<8} {:<6.2f}({:<5.2f}) {:<6.2f}({:<5.2f})"
                  " {:<6.2f}({:<5.2f}) {:<6.2f}({:<4.2f}) {:<6.2f}"
                  "({:<4.2f}) {:<6.2f}({:<4.2f})").format(
        cluster.cname,
        len(cluster.members),
        cluster.medians[0], cluster.stdevs[0],
        cluster.medians[1], cluster.stdevs[1],
        cluster.medians[2], cluster.stdevs[2],
        cluster.medians[3], cluster.stdevs[3],
        cluster.medians[4], cluster.stdevs[4],
        cluster.medians[5], cluster.stdevs[5])
      out_str += ("\n{:>24}  {:<6.2f}{:<7} {:<6.2f}{:<7}"
                  " {:<6.2f}{:<7} {:<6.2f}{:<6} {:<6.2f}"
                  "{:<6} {:<6.2f}{:<6}  {:<6.2}").format(
        SS.getvalue().strip()[13:],
        uc_params_conv[0], "",
        uc_params_conv[1], "",
        uc_params_conv[2], "",
        uc_params_conv[3], "",
        uc_params_conv[4], "",
        uc_params_conv[5], "",
        group["max_angular_difference"]) + "\n\n"

    else:
      singletons.append("".join([("{:<14} {:<11.2f} {:<11.2f} {:<11.2f}"
                                  "{:<12.1f} {:<12.1f} {:<12.1f}").format(
        cluster.pg_composition.keys()[0],
        cluster.members[0].uc[0], cluster.members[0].uc[1],
        cluster.members[0].uc[2], cluster.members[0].uc[3],
        cluster.members[0].uc[4], cluster.members[0].uc[5]),
                                 '\n']))
  out_str += "\nStandard deviations are in brackets."
  explanation = """\nEach cluster:
Input lattice count, with integration Bravais setting space group.
Cluster median with Niggli cell parameters (std dev in brackets).
Highest possible metric symmetry and unit cell using LePage (J Appl Cryst 1982, 15:255) method, maximum delta 3deg."""
  out_str += explanation
  singleton_str = "\n%i singleton%s:" %plural_s(len(singletons))
  singleton_str += "\n\n{:<14} {:<11} {:<11} {:<11}{:<12} {:<12} {:<12}\n".format(
    "Point group",
    "a", "b", "c",      "alpha", "beta", "gamma")
  singleton_str += "".join(singletons)
  n_clusters = len(sub_clusters) - len(singletons)
  out_str = "\n%i cluster%s:" %plural_s(n_clusters) + out_str
  return singleton_str + out_str
