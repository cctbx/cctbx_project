"""f' and f" table lookup given element and wavelength"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.form_factor_query

from cctbx.eltbx import sasaki, henke
import libtbx.phil
from libtbx.utils import Sorry, Usage
import sys

master_params = libtbx.phil.parse("""\
form_factor_query {
  element = None
    .type = str
  wavelength = None
    .type = float
  table = *sasaki henke
    .type = choice(multi = True)
  unit = *angstroms ev kev
    .type = choice
}
""")

def run(args, command_name="phenix.form_factor_query"):
  phil_objects=[]
  argument_interpreter = master_params.command_line_argument_interpreter(
    home_scope="form_factor_query")
  plain_element = None
  plain_wavelength = None
  for arg in args:
    if (arg.find('=') != -1):
      try:
        command_line_params = argument_interpreter.process(arg=arg)
      except KeyboardInterrupt: raise
      except Exception:
        raise Sorry("Unknown keyword: %s" % arg)
      else:
        phil_objects.append(command_line_params)
    else:
      try:
        plain_wavelength = float(arg)
      except KeyboardInterrupt: raise
      except Exception:
        plain_element = arg

  working_params = master_params.fetch(sources=phil_objects)
  params = working_params.extract()
  element = params.form_factor_query.element
  if (element is None):
    element = plain_element
  wavelength = params.form_factor_query.wavelength
  if (wavelength is None):
    wavelength = plain_wavelength
  unit = params.form_factor_query.unit

  if (element is None or wavelength is None):
    raise Usage(
      "%s element=symbol wavelength=float [table=sasaki|henke]"
      " [unit=angstroms|ev|kev]" % command_name)

  for table in ["sasaki", "henke"]:
    if (table not in params.form_factor_query.table): continue
    t = eval(table).table(element)
    print("Information from %s table about %s (Z = %s) at %s %s" % (
      table.capitalize(), t.label(), t.atomic_number(), wavelength,
      {"angstroms": "A"}.get(unit, unit)))
    if (unit == "angstroms"):
      f = t.at_angstrom(wavelength)
    elif (unit == "ev"):
      f = t.at_ev(wavelength)
    elif (unit == "kev"):
      f = t.at_kev(wavelength)
    else:
      raise Sorry("Invalid unit chosen for query: %s" % unit)
    if (f.is_valid_fp()):
      print("fp:  %.6g" % f.fp())
    else:
      print("fp:  unknown")
    if (f.is_valid_fdp()):
      print("fdp: %.6g" % f.fdp())
    else:
      print("fdp: unknown")

if (__name__ == "__main__"):
  run(args=sys.argv[1:])

