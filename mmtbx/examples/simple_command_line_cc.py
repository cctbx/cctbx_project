
"""
Template for writing simple applications which take as input a data file
(amplitudes or intensities plus R-free flags) and a PDB file.  Note that the
script can be modified to also process the covalent geometry of the model by
setting PROCESS_PDB_FILE=True, but this also requires the CCP4 monomer library.
"""

from __future__ import division
import iotbx.phil
import sys

PROCESS_PDB_FILE = False

master_phil = iotbx.phil.parse("""
include scope mmtbx.utils.cmdline_input_phil_str
""", process_includes=True)

def run (args, out=sys.stdout) :
  import mmtbx.utils
  # this wrapper loads the data and flags (or raises an error if additional
  # input is needed), reads the PDB file, optionally processes the geometry,
  # and creates an fmodel object using the data, flags, and xray.structure
  # object from the PDB file.
  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    update_f_part1_for="map",
    args=args,
    master_phil=master_phil,
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
  print >> out, "CC(obs-calc): %.3f" % cc

if (__name__ == "__main__") :
  run(sys.argv[1:])
