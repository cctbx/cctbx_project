from cctbx import miller
from cctbx.array_family import flex

def miller_array_export_as_cns_hkl(self, file_object, file_name=None, info=[]):
  f = file_object
  if (file_name): print >> f, "{ file:", file_name, "}"
  if (self.info() != None):
    print >> f, "{", self.info(), "}"
  for line in info: print >> f, "{", line, "}"
  print >> f, "NREFlections=%d" % self.indices().size()
  if (self.anomalous_flag()):
    print >> f, "ANOMalous=TRUE"
  else:
    print >> f, "ANOMalous=FALSe"
  if (self.sigmas() != None):
    assert isinstance(self.data(), flex.double)
    assert isinstance(self.sigmas(), flex.double)
    if (isinstance(self, miller.intensity_array)):
      f_obs = self.f_sq_as_f()
    else:
      f_obs = self
    print >> f, "DECLare NAME=FOBS  DOMAin=RECIprocal TYPE=REAL END"
    print >> f, "DECLare NAME=SIGMA DOMAin=RECIprocal TYPE=REAL END"
    for i,h in f_obs.indices().items():
      print >> f, "INDEx %d %d %d" % h,
      print >> f, "FOBS=%.6g" % f_obs.data()[i],
      print >> f, "SIGMA=%.6g" % f_obs.sigmas()[i]
  elif (self.is_complex()):
    print >> f, "DECLare NAME=F  DOMAin=RECIprocal TYPE=COMPLEX END"
    a = flex.abs(self.data())
    p = flex.arg(self.data(), 0001)
    for i,h in self.indices().items():
      print >> f, "INDEx %d %d %d" % h,
      print >> f, "F=%.6g %.6g" % (a[i], p[i])
  else:
    if (isinstance(self, flex.double)):
      print >> f, "DECLare NAME=DATA  DOMAin=RECIprocal TYPE=REAL END"
      fmt = "%.6g"
    elif (isinstance(self, flex.int) or isinstance(self, flex.bool)):
      print >> f, "DECLare NAME=DATA  DOMAin=RECIprocal TYPE=INTEger END"
      fmt = "%d"
    else:
      raise RuntimeError, \
        "Cannot write array type %s to CNS reflection file" % type(self.data())
    for i,h in self.indices().items():
      print >> f, "INDEx %d %d %d" % h,
      print >> f, ("DATA="+fmt) % self.data()[i]

miller.array.export_as_cns_hkl = miller_array_export_as_cns_hkl
