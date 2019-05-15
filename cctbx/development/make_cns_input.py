from __future__ import absolute_import, division, print_function
import cctbx.eltbx.xray_scattering
from cctbx import adptbx
from cctbx import eltbx
from cctbx.development import random_structure
from cctbx.development import debug_utils
from libtbx import adopt_init_args
import sys
from six.moves import range

def tst_run_requiring_cns(args, call_back):
  import libtbx.path
  if (libtbx.path.full_command_path(command="cns") is None):
    print("Skipping tests: cns not available.")
  else:
    debug_utils.parse_options_loop_space_groups(
      argv=args, call_back=call_back, show_cpu_times=False)
  from libtbx.utils import format_cpu_times
  print(format_cpu_times())

def write(file_name, cns_input):
  f = open(file_name, "w")
  for line in cns_input:
    print(line, file=f)
  f.close()

def topology(scatterers):
  cns_input = []
  a = cns_input.append
  a("topology")
  for scatterer in scatterers:
    name = scatterer.label.replace(" ","")[-4:]
    a("  residue %s" % name)
    a("    atom %s mass=1 charge=0 {chemical}type=%s end" % (name, name))
    a("  end")
  a("end")
  a("")
  a("segment")
  a("  name=A")
  for scatterer in scatterers:
    name = scatterer.label.replace(" ","")[-4:]
    a("  molecule {res}name=%s number=1 end" % name)
  a("end")
  a("")
  return cns_input

def coordinates(scatterers, xyz_only=False):
  cns_input = []
  a = cns_input.append
  resid = 1
  for scatterer in scatterers:
    x = scatterer.site
    q = scatterer.occupancy
    assert scatterer.flags.use_u_iso_only()
    b = adptbx.u_as_b(scatterer.u_iso)
    gaussian = eltbx.xray_scattering.it1992(scatterer.scattering_type).fetch()
    fp = scatterer.fp
    fdp = scatterer.fdp
    for i in range(3):
      a("do (%s=%.12g) (resid=%d)" % ("xyz"[i], x[i], resid))
    if (not xyz_only):
      a("do (q=%.12g) (resid=%d)" % (q, resid))
      a("do (b=%.12g) (resid=%d)" % (b, resid))
      a("xray")
      a("  scatter (chemical=%s)" % scatterer.label.replace(" ","")[-4:])
      for i in range(4):
        a("    %.6g %.6g" %
          (gaussian.array_of_a()[i], gaussian.array_of_b()[i]))
      a("    %.6g" % (gaussian.c(),))
      a("end")
      a("do (scatter_fp=%.12g) (resid=%d)" % (fp, resid))
      a("do (scatter_fdp=%.12g) (resid=%d)" % (fdp, resid))
    a("")
    resid += 1
  a("coordinates orthogonalize end")
  a("show (name) (all)")
  a("show (chemical) (all)")
  if (not xyz_only):
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
  if (miller_array.is_complex_array()):
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
  for i,h in enumerate(miller_array.indices()):
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

class parameter_nbonds(object):

  def __init__(self, ctofnb=7.5,
                     ctonnb=6.5,
                     cutnb=8.5,
                     irexponent=2,
                     nbxmod=5,
                     rconst=100,
                     repel=0,
                     rexponent=2,
                     tolerance=0.5,
                     wmin=1.5):
    adopt_init_args(self, locals())

class parameter_nonb(object):

  def __init__(self, epsilon=0.1,
                     sigma=1.5,
                     epsilon14=0.1,
                     sigma14=1.3):
    adopt_init_args(self, locals())

def script_vdw_energy(structure, param_nbonds=None, param_nonb=None):
  if (param_nbonds is None): param_nbonds = parameter_nbonds()
  if (param_nonb is None): param_nonb = parameter_nonb()
  cns_input = xray_unit_cell(structure.unit_cell())
  a,e = cns_input.append,cns_input.extend
  e(xray_symmetry(structure.space_group()))
  e(topology(structure.scatterers()))
  e(coordinates(structure.scatterers(), xyz_only=True))
  a("""\
param
  nbonds
    ctofnb=%(ctofnb).6g
    ctonnb=%(ctonnb).6g
    cutnb=%(cutnb).6g
    irexponent=%(irexponent).6g
    nbxmod=%(nbxmod).6g
    rconst=%(rconst).6g
    repel=%(repel).6g
    rexponent=%(rexponent).6g
    tolerance=%(tolerance).6g
    wmin=%(wmin).6g
  end
""" % vars(param_nbonds))
  a("  nonb (all) %(epsilon).6g %(sigma).6g %(epsilon14).6g %(sigma14).6g" %
    vars(param_nonb))
  a("""\

  nbonds
    ?
  end
end

flags exclude * include vdw pvdw end

energy end

show (dx) (all)
show (dy) (all)
show (dz) (all)

distance
  from=(all)
  to=(all)
  cuton=0
  cutoff=%(cutnb).6g
  disp=print
end

stop""" % vars(param_nbonds))
  return cns_input

def exercise(space_group_info,
             anomalous_flag=False,
             d_min=2.,
             verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=("N", "C", "C", "O"),
    random_f_prime_d_min=1.0,
    random_f_double_prime=anomalous_flag,
    random_u_iso=True,
    random_occupancy=True)
  if (0 or verbose):
    structure.show_summary()
  f_calc = structure.structure_factors(
    d_min=d_min,
    algorithm="direct").f_calc()
  for f_obs in [f_calc, abs(f_calc)]:
    script_predict_methods_comparison(d_min, structure, f_obs)
    script_xray_gradients(d_min, f_obs, structure, method="direct")
  script_vdw_energy(structure=structure)

def run_call_back(flags, space_group_info):
  for anomalous_flag in (False, True):
    exercise(space_group_info, anomalous_flag,
             verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
