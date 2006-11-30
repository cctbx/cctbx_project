import libtbx.phil
import libtbx.phil.command_line
from libtbx.utils import Sorry
from cctbx.eltbx import sasaki,henke
import sys

master_params = libtbx.phil.parse("""\
form_factor_query {
  element = Se
    .type = str
  wavelength = 1.5481
    .type = float
  table = *sasaki henke
    .type = choice(multi = True)
  unit = *angstroms ev kev
    .type = choice
}
""")

def run(args):
  phil_objects=[]
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_params=master_params, home_scope="form_factor_query")
  symbol = None
  wavelength = None
  for arg in args:
    if (arg.find('=') != -1):
      try:
        command_line_params = argument_interpreter.process(arg=arg)
      except:
        raise Sorry("Unknown keyword: %s" % arg)
      else:
        phil_objects.append(command_line_params)
    else:
      try:
        wavelength = float(arg)
      except:
        symbol = arg

  working_params = master_params.fetch(sources=phil_objects)
  params = working_params.extract()

  if (wavelength is None):
    wavelength = params.form_factor_query.wavelength
  if (symbol is None):
    symbol = params.form_factor_query.element
  unit = params.form_factor_query.unit

  if ('sasaki' in params.form_factor_query.table):
    t = sasaki.table(symbol)
    print "Information from %s table about %s (Z = %s) at %s %s" % (
      "Sasaki", t.label(), t.atomic_number(), wavelength,
      params.form_factor_query.unit)
    if (unit == "angstroms"):
      f = t.at_angstrom(wavelength)
    elif (unit == "ev"):
      f = t.at_ev(wavelength)
    elif (unit == "kev"):
      f = t.at_kev(wavelength)
    else:
      raise Sorry("Invalid unit chosen for query")
    if (f.is_valid_fp()):
      print "fp:", f.fp()
      if (f.is_valid_fdp()):
        print "fdp:", f.fdp()

  if ('henke' in params.form_factor_query.table):
    t = henke.table(symbol)
    print "Information from %s table about %s (Z = %s) at %s %s" % (
      "Henke", t.label(), t.atomic_number(), wavelength,
      params.form_factor_query.unit)
    if (unit == "angstroms"):
      f = t.at_angstrom(wavelength)
    elif (unit == "ev"):
      f = t.at_ev(wavelength)
    elif (unit == "kev"):
      f = t.at_kev(wavelength)
    else:
      raise Sorry("Invalid unit chosen for query")
    if (f.is_valid_fp()):
      print "fp:", f.fp()
      if (f.is_valid_fdp()):
        print "fdp:", f.fdp()

if (__name__ == "__main__"):
  run(sys.argv[1:])
