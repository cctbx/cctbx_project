
"""
Template for writing simple applications which take as input a data file
(amplitudes or intensities plus R-free flags) and a PDB file.  Note that the
script can be modified to also process the covalent geometry of the model by
setting PROCESS_PDB_FILE=True, but this also requires the CCP4 monomer
library.
"""

from __future__ import absolute_import, division, print_function
import mmtbx.command_line
import sys

PROCESS_PDB_FILE = False

def master_phil():
  return mmtbx.command_line.generic_simple_input_phil()

def run(args, out=sys.stdout):
  # this wrapper loads the data and flags (or raises an error if additional
  # input is needed), reads the PDB file, optionally processes the geometry,
  # and creates an fmodel object using the data, flags, and xray.structure
  # object from the PDB file.
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=master_phil(),
    out=out,
    process_pdb_file=PROCESS_PDB_FILE,
    create_fmodel=True,
    prefer_anomalous=False)
  fmodel = cmdline.fmodel
  pdb_hierarchy = cmdline.pdb_hierarchy
  xray_structure = cmdline.xray_structure
  params = cmdline.params
  f_obs = fmodel.f_obs()
  # the fmodel object will already have the bulk solvent correction and
  # scaling performed when created using the above code, so we can immediately
  # use the f_model array.
  f_calc = abs(fmodel.f_model()) # just amplitudes, please
  assert (len(f_calc.indices()) == len(f_obs.indices()))
  from scitbx.array_family import flex
  cc = flex.linear_correlation(f_obs.data(), f_calc.data()).coefficient()
  print("CC(obs-calc): %.3f" % cc, file=out)

if (__name__ == "__main__"):
  run(sys.argv[1:])
