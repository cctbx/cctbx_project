def run(args):
  from rstbx.simage import create
  work_params = create.process_args(
    args=args,
    extra_phil_str="""\
wavelength_2 = None
  .type = float
""")
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
      print "%3d %3d %3d" % h, "(%8.2f, %8.2f)" % tuple(x[:2]), "%.6g" % s
    print
  else:
    work_params.wavelength = work_params.wavelength_2
    i_calc_2, image_info_2 = create.compute(
      work_params=work_params,
      store_miller_index_i_seqs=True,
      store_spots=True,
      store_signals=True)
    assert i_calc.p1_anom.indices().all_eq(i_calc_2.p1_anom.indices())
    lookup_dict = {}
    for ii,i_seq in enumerate(image_info_2.miller_index_i_seqs):
      lookup_dict[i_seq] = ii
    d_spacings = i_calc.p1_anom.d_spacings().data()
    for i_seq,x,s in zip(
          image_info.miller_index_i_seqs,
          image_info.spots,
          image_info.signals):
      ii = lookup_dict.get(i_seq)
      if (ii is not None):
        x = x[:2]
        x2 = image_info_2.spots[ii][:2]
        s2 = image_info_2.signals[ii]
        h = i_calc.p1_anom.indices()[i_seq]
        print "%3d %3d %3d" % h, " %8.5f" % d_spacings[i_seq], \
          "(%7.1f, %7.1f)" % tuple(x), "%.1f" % s
        print "                     ", \
          "(%7.1f, %7.1f)" % tuple(x2), "%.1f" % s2, " %.1f" % (s2-s)
        print "                     ", \
          "(%7.1f, %7.1f)" % tuple([b-a for a,b in zip(x,x2)])

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
