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
  print "<pre>"
  print "Structure factors calculation\n"
  print "Algorithm: %s\n" % inp.algorithm
  print "Formfactors: %s\n" % inp.form_factors
  print "Minimum d-spacing (highest resolution): %s\n" % inp.d_min
  if (float(inp.d_min) <= 0.):
    raise ValueError, "d-spacing must be greater than zero."
  print

  structure = io_utils.structure_from_inp_pdb(inp, status)
  structure.scattering_dict(table = inp.form_factors,
                            d_min = float(inp.d_min))
  print "Input model summary:\n"
  structure.show_summary()
  structure.scattering_dict().show_summary()
  f_calc = structure.structure_factors(anomalous_flag = False,
                                       algorithm      = inp.algorithm,
                                       d_min          = float(inp.d_min)).f_calc()
  print "\nNumber of Miller indices:", f_calc.indices().size()
  print
  status.in_table = True
  print "    h     k     l           Fcalc      Pcalc"
  for i,h in enumerate(f_calc.indices()):
    print "%5d %5d %5d %15.6f %10.4f" % (
      h + complex_math.abs_arg(f_calc.data()[i], deg=True))
  status.in_table = False
  print

  print "</pre>"
