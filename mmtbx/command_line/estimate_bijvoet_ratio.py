
from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry
from math import sqrt
import sys

master_phil_str = """
element = None
  .type = str
  .help = Anomalously scattering atom type
n_sites = None
  .type = int
fdp = None
  .type = float
wavelength = None
  .type = float
  .help = Data collection wavelength (in Angstroms)
energy = None
  .type = float
  .help = Data collection energy (in eV)
n_res = None
  .type = int
  .help = Number of residues in ASU
mw = None
  .type = float
  .help = Molecular weight of ASU
"""


def run(args, out=sys.stdout, params=None):
  import iotbx.phil
  from cctbx.eltbx import chemical_elements, sasaki
  class interpreter(iotbx.phil.process_command_line_with_files):
    def process_other(self, arg):
      if (len(arg) <= 2):
        if (arg.upper() in chemical_elements.proper_upper_list()):
          return iotbx.phil.parse("element=%s" % arg)
  if (params is None):
    cmdline = interpreter(
      args=args,
      master_phil_string=master_phil_str,
      integer_def="n_sites",
      float_def="fdp",
      usage_string="""\
mmtbx.estimate_bijvoet_ratios Se 20 n_res=1000 wavelength=0.9792

Estimate the Bijvoet ratio for a macromolecular X-ray diffraction experiment
given an anomalous scatterer type and expected asymmetric unit contents.""")
    params = cmdline.work.extract()
  validate_params(params)
  fdp = params.fdp
  if (fdp is None):
    if (params.wavelength is not None):
      fdp = sasaki.table(params.element).at_angstrom(params.wavelength).fdp()
      caption = "%s A" % params.wavelength
    else :
      fdp = sasaki.table(params.element).at_ev(params.energy).fdp()
      caption = "%s eV" % params.energy
  if (params.n_res is not None):
    n_atoms = params.n_res * 110 / 15
  else :
    n_atoms = params.mw / 15
  bijvoet_ratio_acentric = sqrt(2 * params.n_sites / n_atoms) * (fdp / 6.7)
  print("Heavy atom type: %s" % params.element, file=out)
  print("Number of sites: %d" % params.n_sites, file=out)
  print("Approx. # of non-H/D atoms: %d" % n_atoms, file=out)
  if (params.fdp is None):
    print("f'' of %s at %s : %6.3f" % (params.element, caption, fdp), file=out)
  else :
    print("f'' (experimental) : %6.3f" % fdp, file=out)
  print("Expected Bijvoet ratio : %4.1f%%" % \
    (bijvoet_ratio_acentric * 100), file=out)
  return bijvoet_ratio_acentric

def validate_params(params):
  from cctbx.eltbx import chemical_elements
  all_elems = chemical_elements.proper_upper_list()
  if (params.element is None):
    raise Sorry("Element symbol not specified.")
  elif (not params.element.upper() in all_elems):
    raise Sorry("Element symbol '%s' not recognized." % params.element)
  if (params.n_sites is None):
    raise Sorry("Number of sites not specified.")
  if ([params.fdp, params.energy, params.wavelength].count(None) != 2):
    print([params.fdp, params.energy, params.wavelength])
    raise Sorry("Please specify either an X-ray wavelength or energy "+
      "(but not both), or the expected f''.")
  if ([params.n_res, params.mw].count(None) != 1):
    raise Sorry("Please specify either the number of residues or the "+
      "approximate molecular weight (but not both).")
  return True

if (__name__ == "__main__"):
  run(sys.argv[1:])
