from cctbx import miller
from cctbx.array_family import flex

def miller_array_export_as_cns_hkl(self,
      file_object,
      file_name=None,
      info=[],
      array_names=None):
  out = file_object
  if (file_name): print >> out, "{ file:", file_name, "}"
  if (self.info() is not None):
    print >> out, "{", self.info(), "}"
  for line in info: print >> out, "{", line, "}"
  print >> out, "NREFlections=%d" % self.indices().size()
  if (self.anomalous_flag()):
    print >> out, "ANOMalous=TRUE"
  else:
    print >> out, "ANOMalous=FALSe"
  if (self.sigmas() is not None):
    if (array_names is None): array_names = ["FOBS", "SIGMA"]
    else: assert len(array_names) == 2
    assert isinstance(self.data(), flex.double)
    assert isinstance(self.sigmas(), flex.double)
    if (self.is_xray_intensity_array()):
      f_obs = self.f_sq_as_f()
    else:
      f_obs = self
    nf, ns = array_names
    print >> out, "DECLare NAME=%s DOMAin=RECIprocal TYPE=REAL END" % nf
    print >> out, "DECLare NAME=%s DOMAin=RECIprocal TYPE=REAL END" % ns
    for h,f,s in zip(f_obs.indices(),f_obs.data(),f_obs.sigmas()):
      print >> out, "INDEx %d %d %d" % h, "%s=%.6g %s=%.6g" % (nf,f,ns,s)
  elif (self.is_complex_array()):
    if (array_names is None): array_names = ["F"]
    else: assert len(array_names) == 1
    n = array_names[0]
    print >> out, "DECLare NAME=%s  DOMAin=RECIprocal TYPE=COMPLEX END" % n
    for h,a,p in zip(self.indices(),
                     flex.abs(self.data()),
                     flex.arg(self.data(), True)):
      print >> out, "INDEx %d %d %d" % h, "%s=%.6g %.6g" % (n,a,p)
  else:
    if (array_names is None): array_names = ["DATA"]
    else: assert len(array_names) == 1
    if (isinstance(self.data(), flex.double)):
      print >> out, \
        "DECLare NAME=%s  DOMAin=RECIprocal TYPE=REAL END" % array_names[0]
      fmt = "%.6g"
    elif (   isinstance(self.data(), flex.int)
          or isinstance(self.data(), flex.bool)):
      print >> out, \
        "DECLare NAME=%s  DOMAin=RECIprocal TYPE=INTEger END" % array_names[0]
      fmt = "%d"
    else:
      raise RuntimeError, \
        "Cannot write array type %s to CNS reflection file" % type(self.data())
    fmt = array_names[0] + "=" + fmt
    for h,d in zip(self.indices(),self.data()):
      print >> out, "INDEx %d %d %d" % h, fmt % d

miller.array.export_as_cns_hkl = miller_array_export_as_cns_hkl
