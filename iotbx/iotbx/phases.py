def miller_array_as_phases_phs(self, out):
  assert self.is_complex_array()
  ampl = self.amplitudes().data()
  pdeg = self.phases(deg=True).data()
  for h,a,p in zip(self.indices(),ampl,pdeg):
    print >> out, "%4d%4d%4d" % h + "%8.2f" % a + "%8.2f" % 1 + "%8.2f" % p
