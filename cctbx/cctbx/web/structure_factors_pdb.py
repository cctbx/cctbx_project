from scitbx.python_utils import complex_math
from cctbx.web import io_utils
from cctbx.web import cgi_utils

def interpret_form_data(form):
  inp = cgi_utils.inp_from_form(form,
    (("algorithm", ""),
     ("form_factors", ""),
     ("d_min", "")))
  inp.coordinates = cgi_utils.coordinates_from_form(form)
  return inp

def run(server_info, inp, status):
  print "<h2>Structure factor calculation</h2>"
  print "<pre>"
  print "Atomic form factors: %s" % inp.form_factors
  print
  print "Minimum d-spacing (highest resolution): %s" % inp.d_min
  if (float(inp.d_min) <= 0.):
    raise ValueError("d-spacing must be greater than zero.")
  print

  structure = io_utils.structure_from_inp_pdb(inp, status)
  structure.scattering_type_registry(
    table=inp.form_factors,
    d_min=float(inp.d_min))
  print "Input model summary:"
  print
  structure.show_summary()
  print
  structure.scattering_type_registry().show()
  print
  algorithm = inp.algorithm
  if (algorithm == "automatic"):
    if (structure.scatterers().size() <= 100):
      algorithm = "direct"
    else:
      algorithm = None
  elif (algorithm not in ["direct", "fft"]):
    algorithm = None
  f_calc_manager = structure.structure_factors(
    anomalous_flag = False,
    algorithm      = algorithm,
    d_min          = float(inp.d_min))
  f_calc = f_calc_manager.f_calc()
  print "Number of Miller indices:", f_calc.indices().size()
  print
  print "Structure factor algorithm:", f_calc_manager.algorithm(verbose=True)
  print
  print "    h     k     l           Fcalc      Pcalc"
  for h,f in zip(f_calc.indices(), f_calc.data()):
    print "%5d %5d %5d %15.6f %8.2f" % (h + complex_math.abs_arg(f, deg=True))
  print

  print "</pre>"
