from cctbx import adptbx
from cctbx.eltbx.caasf import it1992
import sys, os

def check_cns_availability():
  if (not os.environ.has_key("CNS_INST")):
    print "CNS not available."
    sys.exit(1)

def write(file_name, cns_input):
  f = open(file_name, "w")
  for line in cns_input:
    print >> f, line
  f.close()

def topology(scatterers):
  cns_input = []
  a = cns_input.append
  a("topology")
  for scatterer in scatterers:
    name = scatterer.label
    a("  residue %s" % name)
    a("    atom %s mass=1 charge=0 {chemical}type=%s end" % (name, name))
    a("  end")
  a("end")
  a("")
  a("segment")
  a("  name=A")
  for scatterer in scatterers:
    name = scatterer.label
    a("  molecule {res}name=%s number=1 end" % name)
  a("end")
  a("")
  return cns_input

def coordinates(scatterers):
  cns_input = []
  a = cns_input.append
  resid = 1
  for scatterer in scatterers:
    x = scatterer.site
    q = scatterer.occupancy
    assert not scatterer.anisotropic_flag
    b = adptbx.u_as_b(scatterer.u_iso)
    caasf = it1992(scatterer.caasf.label())
    fp = scatterer.fp
    fdp = scatterer.fdp
    for i in xrange(3):
      a("do (%s=%.12g) (resid=%d)" % ("xyz"[i], x[i], resid))
    a("do (q=%.12g) (resid=%d)" % (q, resid))
    a("do (b=%.12g) (resid=%d)" % (b, resid))
    a("xray")
    a("  scatter (chemical=%s)" % (scatterer.label,))
    for i in xrange(4):
      a("    %.6g %.6g" % (caasf.a()[i], caasf.b()[i]))
    a("    %.6g" % (caasf.c(),))
    a("end")
    a("do (scatter_fp=%.12g) (resid=%d)" % (fp, resid))
    a("do (scatter_fdp=%.12g) (resid=%d)" % (fdp, resid))
    a("")
    resid += 1
  a("coordinates orthogonalize end")
  a("show (name) (all)")
  a("show (chemical) (all)")
  a("show (scatter_fp) (all)")
  a("show (scatter_fdp) (all)")
  a("")
  return cns_input

def xray_unit_cell(unit_cell):
  cns_input = []
  a = cns_input.append
  a("xray")
  a("  a=%.12g b=%.12g c=%.12g alpha=%.12g beta=%.12g gamma=%.12g"
    % unit_cell.parameters())
  a("end")
  a("")
  return cns_input

def xray_symmetry(space_group):
  cns_input = []
  a = cns_input.append
  a("xray")
  for s in space_group:
    a("  symm=(" + s.mod_positive().as_xyz() + ")")
  a("end")
  a("")
  return cns_input

def xray_structure(structure):
  cns_input = xray_unit_cell(structure.unit_cell())
  e = cns_input.extend
  e(xray_symmetry(structure.space_group()))
  e(topology(structure.scatterers()))
  e(coordinates(structure.scatterers()))
  return cns_input

def xray_anomalous(flag):
  cns_input = []
  a = cns_input.append
  a("xray")
  if (flag):
    a("  anomalous=true")
  else:
    a("  anomalous=false")
  a("end")
  a("")
  return cns_input

def xray_declare(names, miller_array):
  assert len(names) in (1,2)
  if (len(names) == 2): assert miller_array.sigmas() is not None
  if (miller_array.is_complex()):
    type = "complex"
  else:
    type = "real"
  cns_input = []
  a = cns_input.append
  a("xray")
  a("  declare name=%s domain=reciprocal type=%s end" % (names[0], type))
  if (len(names) == 2):
    a("  declare name=%s domain=reciprocal type=real end" % names[1])
  a("end")
  a("")
  return cns_input

def xray_reflection(names, miller_array):
  assert len(names) in (1,2)
  if (len(names) == 2): assert miller_array.sigmas() is not None
  cns_input = []
  a = cns_input.append
  a("xray")
  a("  nreflections %d" % miller_array.indices().size())
  a("  reflection")
  data = miller_array.data()
  sigmas = miller_array.sigmas()
  for i,h in miller_array.indices().items():
    s = "    index %d %d %d" % h + " %s %s" % (names[0], str(data[i]))
    if (sigmas is not None):
      s += " %s %s" % (names[1], str(sigmas[i]))
    a(s)
  a("  end")
  a("end")
  a("")
  return cns_input

def xray_generate(d_max, d_min):
  cns_input = []
  a = cns_input.append
  a("xray")
  a("  generate %s %s" % (str(d_max), str(d_min)))
  a("end")
  a("")
  return cns_input

def xray_mapresolution(d_min):
  cns_input = []
  a = cns_input.append
  a("xray")
  a("  mapresolution %s" % str(d_min))
  a("end")
  a("")
  return cns_input

def xray_method(method):
  assert method in ["direct", "fft"]
  cns_input = []
  a = cns_input.append
  a("xray")
  a("  method %s" % method)
  a("end")
  a("")
  return cns_input

def xray_predict(reciprocal_space_array):
  cns_input = []
  a = cns_input.append
  a("""\
xray
  declare name=%(reciprocal_space_array)s domain=reciprocal type=complex end
  predict
    mode=reciprocal
    to=%(reciprocal_space_array)s
    selection=(all)
    atomselection=(all)
  end
end
""" % vars())
  return cns_input

def script_predict_methods_comparison(d_min, structure, f_obs):
  cns_input = xray_structure(structure)
  a,e = cns_input.append,cns_input.extend
  e(xray_anomalous(f_obs.anomalous_flag()))
  e(xray_declare(["fobs"], f_obs))
  e(xray_reflection(["fobs"], f_obs))
  e(xray_mapresolution(d_min))
  e(xray_method("direct"))
  e(xray_predict("fcalc_dir"))
  e(xray_method("fft"))
  e(xray_predict("fcalc_fft"))
  a("xray")
  a("  write reflections output=tmp.hkl fcalc_dir fcalc_fft end")
  a("end")
  a("")
  a("stop")
  return cns_input

def xray_gradients():
  return ["""\
xray
  declare name=fcalc domain=reciprocal type=complex end
  declare name=fpart domain=reciprocal type=complex end
  declare name=weight domain=reciprocal type=real end
  do (weight=1) (all)
  do (fpart=0) (all)
  associate reset
  associate fcalc (all)
  target=(resi(amplitude(fobs),(fcalc+fpart), weight))
  dtarget=(dresi(amplitude(fobs),(fcalc+fpart), weight))
  tselection (all)
end

flags exclude * include xref end

energy end

show (dx) (all)
show (dy) (all)
show (dz) (all)
"""]

def script_xray_gradients(d_min, f_obs, structure, method):
  cns_input = xray_structure(structure)
  a,e = cns_input.append,cns_input.extend
  e(xray_anomalous(f_obs.anomalous_flag()))
  e(xray_declare(["fobs"], f_obs))
  e(xray_reflection(["fobs"], f_obs))
  e(xray_mapresolution(d_min))
  e(xray_method(method))
  e(xray_gradients())
  a("stop")
  return cns_input
