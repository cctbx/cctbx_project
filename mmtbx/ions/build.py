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
  print ions

def add_ion(atom, xray_structure, identity):
  pass
