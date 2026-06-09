"""mmtbx.naiad: H-bond-aware placement of hydrogens on water residues."""

from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from mmtbx.programs import water_protonation

if __name__ == "__main__":
  run_program(water_protonation.Program)
