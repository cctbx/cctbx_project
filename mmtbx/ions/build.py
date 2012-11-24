"""
Deals with modifying a structure to include unbuilt and misidentified ions.
"""

from __future__ import division
import sys

def find_and_build_ions(manager, ions, debug = True, out = sys.stdout):
  # Build in the identified ions
  for i_seq, final_choices in ions:
    if len(final_choices) < 1:
      continue
    atom = manager.pdb_atoms[i_seq]
    if len(final_choices) > 1:
      # Ambiguous results
      pass
    else:
      final_choice = final_choices[0]
      if debug:
        print >> out, "Modifying %s to become %s%+d" % \
          (atom.id_str(), final_choice.element, final_choice.charge)

      # Modify the atom object
      _change_identity(atom = atom,
                       element = final_choice.element,
                       charge = str(final_choice.charge))
      # Modify the scattering type
      manager.xray_structure.scatterers()[i_seq].scattering_type = \
        "%s%+d" % (final_choice.element, final_choice.charge)

def _change_identity(atom, element, charge):
  atom.element = element
  atom.charge = charge
  atom.parent().resname = element
  atom.name = element
  atom.segid = "ION"
