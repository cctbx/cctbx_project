"""
Determination of sets of space groups that cannot be distinguished by
inspection of systematic absences.

From first principles; inelegant theoretically, but compact and practical.

See also: International Tables for Crystallography, Volume A, section 3.
"""

def run(args):
  assert args in [[], ["python"], ["c++"]]
  #
  sgno_list_by_index_list_by_cs = {}
  from cctbx import sgtbx
  from cctbx import miller
  for symbols in sgtbx.space_group_symbol_iterator():
    psgi = sgtbx.space_group_info(symbols.universal_hermann_mauguin()) \
      .primitive_setting()
    p_indices = miller.index_generator(
      space_group_type=psgi.type(),
      anomalous_flag=False,
      max_index=[4]*3).to_array()
        # 4 is the smallest value leading to correct results; any larger
        # value will work, too, but will make this procedure slower
    p1_indices = miller.expand_to_p1_iselection(
      space_group=psgi.group(),
      anomalous_flag=False,
      indices=p_indices,
      build_iselection=False).indices
    from cctbx.array_family import flex
    sort_perm = flex.sort_permutation(
      data=miller.index_span(p1_indices).pack(p1_indices))
    p1_indices = p1_indices.select(sort_perm)
    index_list = tuple(p1_indices)
    sgno = psgi.type().number()
    sgno_list_by_index_list = sgno_list_by_index_list_by_cs \
      .setdefault(symbols.crystal_system(), {})
    sgno_list_by_index_list.setdefault(index_list, []).append(sgno)
  from scitbx.graph import tardy_tree
  cluster_manager = tardy_tree.cluster_manager(n_vertices=231)
  for cs,sgno_list_by_index_list in sgno_list_by_index_list_by_cs.items():
    for sgno_list in sgno_list_by_index_list.values():
      i = sgno_list[0]
      for j in sgno_list[1:]:
        cluster_manager.connect_vertices(i=i, j=j, optimize=True)
  cluster_manager.tidy()
  #
  # everything below is just to format the results
  #
  if (args == []):
    for cluster in cluster_manager.clusters:
      if (len(cluster) == 1): break
      print cluster
  else:
    note = ("""\
Output of: cctbx/examples/find_sys_abs_equiv_space_groups.py %s
If you have to edit this table, please send email to: cctbx@cci.lbl.gov
""" % args[0]).splitlines()
    #
    if (args == ["python"]):
      print "space_group_numbers = ["
      for line in note:
        print "  #", line
      ci = cluster_manager.cluster_indices
      cl = cluster_manager.clusters
      for sgno in xrange(231):
        cluster = list(cl[ci[sgno]])
        cluster.remove(sgno)
        if (len(cluster) == 0): s = "None"
        else:                   s = str(tuple(cluster))
        if (sgno == 230): comma = ""
        else:             comma = ","
        print "  %s%s" % (s, comma)
      print "]"
    else:
      print """\
#ifndef CCTBX_SGTBX_SYS_ABS_EQUIV_H
#define CCTBX_SGTBX_SYS_ABS_EQUIV_H

namespace cctbx { namespace sgtbx { namespace sys_abs_equiv {
"""
      data = []
      ci = cluster_manager.cluster_indices
      cl = cluster_manager.clusters
      for line in note:
        print "  //", line
      for sgno in xrange(231):
        cluster = list(cl[ci[sgno]])
        cluster.remove(sgno)
        if (len(cluster) == 0):
          data.append("0")
        else:
          cid = "data_%03d" % sgno
          data.append(cid)
          print "  static const unsigned %s[] = {%d, %s};" % (
            cid, len(cluster), ", ".join([str(i) for i in cluster]))
      print ""
      print "  static const unsigned* space_group_numbers[] = {"
      print "   ", ",\n    ".join(data)
      print """\
    };

}}}

#endif // GUARD"""

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
