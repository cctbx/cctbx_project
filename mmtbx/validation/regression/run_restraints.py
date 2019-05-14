
from __future__ import absolute_import, division, print_function
from mmtbx.monomer_library import pdb_interpretation
from mmtbx.validation import restraints
import sys

def run(args, log=sys.stdout):
  processed_pdb_file = pdb_interpretation.run(
    args=args,
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    log=log)
  xray_structure = processed_pdb_file.xray_structure()
  if xray_structure is None :
    raise Sorry("Could not calculate X-ray structure from this PDB file. "+
      "This is probably due to missing symmetry information (CRYST1 record.")
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies=False)
  chain_proxies = processed_pdb_file.all_chain_proxies
  pdb_hierarchy = chain_proxies.pdb_hierarchy
  result = restraints.combined(
    pdb_hierarchy=chain_proxies.pdb_hierarchy,
    xray_structure=xray_structure,
    geometry_restraints_manager=geometry,
    ignore_hd=True)
  result.show(out=log, prefix="  ")

if (__name__ == "__main__"):
  run(sys.argv[1:])
