from __future__ import absolute_import, division, print_function
from iotbx.cns.crystal_symmetry_utils import crystal_symmetry_as_sg_uc
from cctbx.array_family import flex
from six.moves import range
from six.moves import zip

def crystal_symmetry_as_cns_comments(crystal_symmetry, out):
  if (   crystal_symmetry.unit_cell() is not None
      or crystal_symmetry.space_group_info() is not None):
    print("{ %s }" % crystal_symmetry_as_sg_uc(
      crystal_symmetry=crystal_symmetry), file=out)

def export_as_cns_hkl(self,
      file_object,
      file_name,
      info,
      array_names,
      r_free_flags):
  out = file_object
  if (file_name): print("{ file:", file_name, "}", file=out)
  if (self.info() is not None):
    print("{", self.info(), "}", file=out)
  crystal_symmetry_as_cns_comments(crystal_symmetry=self, out=out)
  for line in info: print("{", line, "}", file=out)
  print("NREFlections=%d" % self.indices().size(), file=out)
  if (self.anomalous_flag()):
    print("ANOMalous=TRUE", file=out)
  else:
    print("ANOMalous=FALSe", file=out)
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
    print("DECLare NAME=%s DOMAin=RECIprocal TYPE=REAL END" % nf, file=out)
    print("DECLare NAME=%s DOMAin=RECIprocal TYPE=REAL END" % ns, file=out)
    if (r_free_flags is None):
      for h,f,s in zip(f_obs.indices(),f_obs.data(),f_obs.sigmas()):
        print("INDEx %d %d %d" % h, "%s= %.6g %s= %.6g" % (nf,f,ns,s), file=out)
    else:
      assert r_free_flags.indices().all_eq(f_obs.indices())
      print("DECLare NAME=TEST DOMAin=RECIprocal TYPE=INTE END", file=out)
      for h,f,s,t in zip(f_obs.indices(),f_obs.data(),f_obs.sigmas(),
                         r_free_flags.data()):
        print("INDEx %d %d %d" % h, "%s= %.6g %s= %.6g" % (nf,f,ns,s),\
          "TEST= %d" % int(t), file=out)
  elif (self.is_complex_array()):
    if (array_names is None): array_names = ["F"]
    else: assert len(array_names) == 1
    assert r_free_flags is None
    n = array_names[0]
    print("DECLare NAME=%s  DOMAin=RECIprocal TYPE=COMPLEX END" % n, file=out)
    for h,a,p in zip(self.indices(),
                     flex.abs(self.data()),
                     flex.arg(self.data(), True)):
      print("INDEx %d %d %d" % h, "%s= %.6g %.6g" % (n,a,p), file=out)
  elif (self.is_hendrickson_lattman_array()):
    if (array_names is None): array_names = ["PA", "PB", "PC", "PD"]
    else: assert len(array_names) == 4
    assert r_free_flags is None
    for i in range(4):
      print("DECLare NAME=%s  DOMAin=RECIprocal TYPE=REAL END" %(
        array_names[i]), file=out)
    print("GROUp  TYPE=HL", file=out)
    for i in range(4):
      print("     OBJEct=%s" %(array_names[i]), file=out)
    print("END", file=out)
    for h,hl in zip(self.indices(), self.data()):
      print("INDEx %d %d %d" % h, end=' ', file=out)
      print("%s= %.6g" % (array_names[0], hl[0]), end=' ', file=out)
      print("%s= %.6g" % (array_names[1], hl[1]), end=' ', file=out)
      print("%s= %.6g" % (array_names[2], hl[2]), end=' ', file=out)
      print("%s= %.6g" % (array_names[3], hl[3]), file=out)
  else:
    if (array_names is None): array_names = ["DATA"]
    else: assert len(array_names) == 1
    assert r_free_flags is None
    if (isinstance(self.data(), flex.double)):
      print("DECLare NAME=%s  DOMAin=RECIprocal TYPE=REAL END" % array_names[0], file=out)
      fmt = "%.6g"
    elif (   isinstance(self.data(), flex.int)
          or isinstance(self.data(), flex.bool)):
      print("DECLare NAME=%s  DOMAin=RECIprocal TYPE=INTEger END" % array_names[0], file=out)
      fmt = "%d"
    else:
      raise RuntimeError("Cannot write array type %s to CNS reflection file" % type(self.data()))
    fmt = array_names[0] + "= " + fmt
    for h,d in zip(self.indices(),self.data()):
      print("INDEx %d %d %d" % h, fmt % d, file=out)
