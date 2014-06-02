__author__ = 'zeldin'
""" Utilitites for dealing with lists of clusters. """

def unit_cell_info(sub_clusters):
      # 3. print out some information that is useful.
      out_str = "{} clusters.".format(len(sub_clusters))
      out_str += "\n{:^5} {:^14} {:<11} {:<11} {:<11} {:<12} {:<12} {:<12}".format(
        "C_id",
        "Num in cluster",
        "Med_a", "Med_b", "Med_c",
        "Med_alpha", "Med_beta", "Med_gamma")
      singletons = []
      for cluster in sub_clusters:
        if len(cluster.members) != 1:

          sorted_pg_comp = sorted(cluster.pg_composition.items(),
                                  key=lambda x: -1 * x[1])
          pg_strings = ["{} in {}".format(pg[1], pg[0])
                        for pg in sorted_pg_comp]
          point_group_string = ", ".join(pg_strings) + "."
          out_str += ("\n{:^5} {:^14} {:<5.1f}({:<4.1f}) {:<5.1f}({:<4.1f})"
                      " {:<5.1f}({:<4.1f}) {:<6.2f}({:<4.2f}) {:<6.2f}"
                      "({:<4.2f}) {:<6.2f}({:<4.2f})").format(
            cluster.cname,
            len(cluster.members),
            cluster.medians[0], cluster.stdevs[0],
            cluster.medians[1], cluster.stdevs[1],
            cluster.medians[2], cluster.stdevs[2],
            cluster.medians[3], cluster.stdevs[3],
            cluster.medians[4], cluster.stdevs[4],
            cluster.medians[5], cluster.stdevs[5])
          out_str += "\n" + point_group_string
        else:
          singletons.append("".join([("{:<14} {:<11.1f} {:<11.1f} {:<11.1f}"
                                      "{:<12.1f} {:<12.1f} {:<12.1f}").format(
            cluster.pg_composition.keys()[0],
            cluster.members[0].uc[0], cluster.members[0].uc[1],
            cluster.members[0].uc[2], cluster.members[0].uc[3],
            cluster.members[0].uc[4], cluster.members[0].uc[5]),
                                     '\n']))
      out_str += "\nStandard deviations are in brackets."
      out_str += "\n" + str(len(singletons)) + " singletons:"
      out_str += "\n{:^14} {:<11} {:<11} {:<11} {:<12} {:<12} {:<12}".format(
        "Point group",
        "a", "b", "c",
        "alpha", "beta", "gamma")
      out_str += "".join(singletons)
      return out_str
