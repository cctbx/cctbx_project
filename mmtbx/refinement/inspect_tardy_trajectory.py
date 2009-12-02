from libtbx import easy_pickle
import sys, os
op = os.path

def run(args):
  assert len(args) == 1, "trajectory_directory"
  traj_dir = args[0]
  labels = easy_pickle.load(file_name=op.join(traj_dir, "labels"))
  fn = op.join(traj_dir, "xray_structure")
  if (op.isfile(fn)):
    xray_structure = easy_pickle.load(file_name=fn)
  sites_file_names_by_index = {}
  for fn in os.listdir(traj_dir):
    if (not fn.startswith("sites_")): continue
    try: i = int(fn[6:])
    except ValueError: continue
    assert i not in sites_file_names_by_index
    sites_file_names_by_index[i] = fn
  assert 0 in sites_file_names_by_index
  sites_0 = easy_pickle.load(
    file_name=op.join(traj_dir, sites_file_names_by_index[0]))
  max_i = max(sites_file_names_by_index.keys())
  sites_final = easy_pickle.load(
    file_name=op.join(traj_dir, sites_file_names_by_index[i]))
  print "rms difference start vs. final: %.3f" % \
    sites_final.rms_difference(sites_0)
  #
  import scitbx.math.superpose
  ls = scitbx.math.superpose.least_squares_fit(
    reference_sites=sites_0, other_sites=sites_final)
  print "ls r:", ls.r
  fm = scitbx.math.r3_rotation_axis_and_angle_from_matrix(r=ls.r)
  print "ls r axis, angle:", fm.axis, "%.2f" % fm.angle(deg=True)
  print "ls t:", ls.t.transpose()
  sites_fit = ls.other_sites_best_fit()
  print "rms difference start vs. best fit: %.3f" % \
    sites_fit.rms_difference(sites_0)
  #
  from cStringIO import StringIO
  out = StringIO()
  fit_dists_sq = (sites_fit - sites_0).dot()
  from scitbx.array_family import flex
  perm = flex.sort_permutation(fit_dists_sq, reverse=True)
  for lbl,fit_d,final_d in zip(
        flex.std_string(labels).select(perm),
        flex.sqrt(fit_dists_sq.select(perm)),
        flex.sqrt((sites_final - sites_0).dot().select(perm))):
    print >> out, lbl, "%.2f %.2f" % (fit_d, final_d)
  if (0):
    sys.stdout.write(out.getvalue())
  #
  if (xray_structure is not None):
    from cctbx import crystal
    site_symmetry_table = xray_structure.site_symmetry_table()
    for i_seq in site_symmetry_table.special_position_indices():
      print "special position:", labels[i_seq]
      try:
        crystal.correct_special_position(
          crystal_symmetry = xray_structure,
          special_op       = site_symmetry_table.get(i_seq).special_op(),
          site_cart        = sites_final[i_seq])
      except AssertionError, e:
        from libtbx.str_utils import prefix_each_line
        print prefix_each_line(prefix="  ", lines_as_one_string=str(e))
        print

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
