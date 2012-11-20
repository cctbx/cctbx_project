"""
Deals with modifying a structure to include unbuilt and misidentified ions.
"""

from __future__ import division
import sys
from libtbx import Auto

def find_and_build_ions(manager, out = sys.stdout, debug = True,
                        show_only_map_outliers = True, candidates = Auto):
  ions = manager.analyze_waters(out = out, debug = debug,
                                show_only_map_outliers = show_only_map_outliers,
                                candidates = candidates)
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
      # print i_seq, final_choice
      _change_identity(atom = atom,
                       element = final_choice.element,
                       charge = str(final_choice.charge))
  return ions

def _change_identity(atom, element, charge):
  atom.element = element
  atom.charge = charge
  atom.fetch_labels().resname = element
  atom.name = element
  atom.segid = "ION"
  # print atom.id_str()
