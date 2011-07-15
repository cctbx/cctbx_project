from iotbx.cns.crystal_symmetry_utils import crystal_symmetry_as_sg_uc
from cctbx.array_family import flex

def crystal_symmetry_as_cns_comments(crystal_symmetry, out):
  if (   crystal_symmetry.unit_cell() is not None
      or crystal_symmetry.space_group_info() is not None):
    print >> out, "{ %s }" % crystal_symmetry_as_sg_uc(
      crystal_symmetry=crystal_symmetry)

def export_as_cns_hkl(self,
      file_object,
      file_name,
      info,
      array_names,
      r_free_flags):
  out = file_object
  if (file_name): print >> out, "{ file:", file_name, "}"
  if (self.info() is not None):
    print >> out, "{", self.info(), "}"
  crystal_symmetry_as_cns_comments(crystal_symmetry=self, out=out)
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
    if (r_free_flags is None):
      for h,f,s in zip(f_obs.indices(),f_obs.data(),f_obs.sigmas()):
        print >> out, "INDEx %d %d %d" % h, "%s= %.6g %s= %.6g" % (nf,f,ns,s)
    else:
      assert r_free_flags.indices().all_eq(f_obs.indices())
      print >> out, "DECLare NAME=TEST DOMAin=RECIprocal TYPE=INTE END"
      for h,f,s,t in zip(f_obs.indices(),f_obs.data(),f_obs.sigmas(),
                         r_free_flags.data()):
        print >> out, "INDEx %d %d %d" % h, "%s= %.6g %s= %.6g" % (nf,f,ns,s),\
          "TEST= %d" % int(t)
  elif (self.is_complex_array()):
    if (array_names is None): array_names = ["F"]
    else: assert len(array_names) == 1
    assert r_free_flags is None
    n = array_names[0]
    print >> out, "DECLare NAME=%s  DOMAin=RECIprocal TYPE=COMPLEX END" % n
    for h,a,p in zip(self.indices(),
                     flex.abs(self.data()),
                     flex.arg(self.data(), True)):
      print >> out, "INDEx %d %d %d" % h, "%s= %.6g %.6g" % (n,a,p)
  elif (self.is_hendrickson_lattman_array()):
    if (array_names is None): array_names = ["PA", "PB", "PC", "PD"]
    else: assert len(array_names) == 4
    assert r_free_flags is None
    for i in range(4):
      print >> out, "DECLare NAME=%s  DOMAin=RECIprocal TYPE=REAL END" %(
        array_names[i])
    print >> out, "GROUp  TYPE=HL"
    for i in range(4):
      print >> out, "     OBJEct=%s" %(array_names[i])
    print >> out, "END"
    for h,hl in zip(self.indices(), self.data()):
      print >> out, "INDEx %d %d %d" % h,
      print >> out, "%s= %.6g" % (array_names[0], hl[0]),
      print >> out, "%s= %.6g" % (array_names[1], hl[1]),
      print >> out, "%s= %.6g" % (array_names[2], hl[2]),
      print >> out, "%s= %.6g" % (array_names[3], hl[3])
  else:
    if (array_names is None): array_names = ["DATA"]
    else: assert len(array_names) == 1
    assert r_free_flags is None
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
    fmt = array_names[0] + "= " + fmt
    for h,d in zip(self.indices(),self.data()):
      print >> out, "INDEx %d %d %d" % h, fmt % d
