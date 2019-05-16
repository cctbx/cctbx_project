from __future__ import absolute_import, division, print_function
from six.moves import zip
def run(args):
  from rstbx.simage import create
  from scitbx.array_family import flex
  work_params = create.process_args(args=args)
  i_calc, image_info = create.compute(
    work_params=work_params,
    store_miller_index_i_seqs=True,
    store_spots=True,
    store_signals=True)
  if (work_params.wavelength_2 is None):
    for i_seq,x,s in zip(
          image_info.miller_index_i_seqs,
          image_info.spots,
          image_info.signals):
      h = i_calc.p1_anom.indices()[i_seq]
      print("%3d %3d %3d" % h, "(%8.2f, %8.2f)" % tuple(x[:2]), "%.6g" % s)
    print()
  else:
    i_calc_2, image_info_2 = create.compute(
      work_params=work_params,
      use_wavelength_2=True,
      store_miller_index_i_seqs=True,
      store_spots=True,
      store_signals=True)
    assert i_calc.p1_anom.indices().all_eq(i_calc_2.p1_anom.indices())
    lookup_dict = {}
    for ii,i_seq in enumerate(image_info_2.miller_index_i_seqs):
      lookup_dict[i_seq] = ii
    d_spacings = i_calc.p1_anom.d_spacings().data()
    d_array = flex.double()
    streak_array = flex.double()
    i_mean_array = flex.double()
    for i_seq,x,s in zip(
          image_info.miller_index_i_seqs,
          image_info.spots,
          image_info.signals):
      ii = lookup_dict.get(i_seq)
      if (ii is not None):
        x = tuple(x[:2])
        x2 = tuple(image_info_2.spots[ii][:2])
        s2 = image_info_2.signals[ii]
        dx = tuple([b-a for a,b in zip(x,x2)])
        streak = (dx[0]**2 + dx[1]**2)**0.5
        h = i_calc.p1_anom.indices()[i_seq]
        print("%3d %3d %3d" % h, " %8.5f" % d_spacings[i_seq], \
          "(%7.1f, %7.1f)" % x, "%.1f" % s)
        print("                     ", \
          "(%7.1f, %7.1f)" % x2, "%.1f" % s2, " %.1f" % (s2-s))
        print("                      (%7.1f, %7.1f)" % dx, " %.1f" % streak)
        d_array.append(d_spacings[i_seq])
        streak_array.append(streak)
        i_mean_array.append((s+s2)*0.5)
    print()
    print("Writing file: d_vs_streak.xy")
    f = open("d_vs_streak.xy", "w")
    for x,y in zip(d_array, streak_array):
      print(x, y, file=f)
    del f
    print("Writing file: i_vs_streak.xy")
    perm = flex.sort_permutation(i_mean_array)
    f = open("i_vs_streak.xy", "w")
    for x,y in zip(i_mean_array.select(perm), streak_array.select(perm)):
      print(x, y, file=f)
    del f
  _ = i_calc.p1_anom
  d = _.d_spacings().data()
  i = _.data()
  print("Writing file: d_vs_i_calc.xy")
  f = open("d_vs_i_calc.xy", "w")
  for x,y in zip(d, i):
    print(x, y, file=f)
  del f
  print("Writing file: i_calc_vs_d.xy")
  perm = flex.sort_permutation(i)
  f = open("i_calc_vs_d.xy", "w")
  for x,y in zip(i.select(perm), d.select(perm)):
    print(x, y, file=f)
  del f
  print()
  _ = i_calc.p1_anom
  for cutoff in [0.1, 0.05, 0.02, 0.01]:
    print("Completeness with i_calc >= %.2f:" % cutoff)
    visible = _.select(_.data() >= cutoff)
    a, b = visible.indices().size(), _.indices().size()
    print("  Overall: %d / %d = %.3f" % (a, b, a/max(1,b)))
    visible.setup_binner(n_bins=10)
    visible.completeness(use_binning=True).show(prefix="  ")
    print()

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
