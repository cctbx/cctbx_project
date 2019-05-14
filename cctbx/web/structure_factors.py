from __future__ import absolute_import, division, print_function
from libtbx import complex_math
from cctbx.web import io_utils
from cctbx.web import cgi_utils

def interpret_form_data(form):
  inp = cgi_utils.inp_from_form(form,
    (("ucparams", "1 1 1 90 90 90"),
     ("sgsymbol", "P1"),
     ("convention", ""),
     ("d_min", "1"),
     ("min_distance_sym_equiv", "0.5"),
     ("algorithm", "automatic"),
     ("coor_type", None)))
  inp.coordinates = cgi_utils.coordinates_from_form(form)
  return inp

def run(server_info, inp, status):
  print("<pre>")
  special_position_settings = io_utils.special_position_settings_from_inp(inp)
  special_position_settings.show_summary()
  print("Minimum distance between symmetrically equivalent sites:", end=' ')
  print(special_position_settings.min_distance_sym_equiv())
  print()

  d_min = float(inp.d_min)
  print("Minimum d-spacing:", d_min)
  if (d_min <= 0.):
    raise ValueError("d-spacing must be greater than zero.")
  print()

  structure = io_utils.structure_from_inp(
    inp, status, special_position_settings)
  algorithm = inp.algorithm
  if (algorithm == "automatic"):
    if (structure.scatterers().size() <= 100):
      algorithm = "direct"
    else:
      algorithm = None
  elif (algorithm not in ["direct", "fft"]):
    algorithm = None
  f_calc_manager = structure.structure_factors(
    anomalous_flag=False, d_min=d_min, algorithm=algorithm)
  f_calc = f_calc_manager.f_calc()
  structure.scattering_type_registry().show()
  print()
  print("Number of Miller indices:", f_calc.indices().size())
  print()
  print("Structure factor algorithm:", f_calc_manager.algorithm(verbose=True))
  print()
  print("</pre><table border=2 cellpadding=2>")
  status.in_table = True
  print("<tr>")
  print("<th>hkl<th>Amplitude<th>Phase")
  for i,h in enumerate(f_calc.indices()):
    print("<tr>")
    print("<td>%3d %3d %3d<td>%.6g<td align=right>%.3f" % (
      h + complex_math.abs_arg(f_calc.data()[i], deg=True)))
  print("</table><pre>")
  status.in_table = False
  print()

  print("</pre>")
