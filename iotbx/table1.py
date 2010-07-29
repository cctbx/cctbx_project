
from libtbx.utils import Sorry
from libtbx import adopt_init_args
import sys, os

keywords = [
  ("wavelength", "Wavelength"),
  ("d_max_min", "Resolution range"),
  ("space_group", "Space group"),
  ("unit_cell", "Unit cell"),
  ("n_refl_all", "Total reflections"),
  ("n_refl", "Unique reflections"),
  ("mult", "Multiplicity"),
  ("completeness", "Completeness"),
  ("i_over_sigma", "I/sigma(I)"),
  ("wilson_b", "Wilson B-factor"),
  ("r_sym", "R-sym"),
  ("r_work", "R-factor"),
  ("r_free", "R-free"),
  ("n_residues", "Protein residues"),
  ("n_wat", "Water molecules"),
  ("rms_bonds", "RMS(bonds)"),
  ("rms_angles", "RMS(angles)"),
  ("rama_fav", "Ramachandran favored"),
  ("rama_outlier", "Ramachandran outliers"),
  ("adp_protein", "Average B-factor (protein)"),
  ("adp_water", "Average B-factor (water)"),
]

class table1 (object) :
  def __init__ (self, **kwds) :
    adopt_init_args(self, locals())

  def add_outer_shell (self, **kwds) :
    self.outer_shell = table1(**kwds)

def format_table (tables, params, out=sys.stdout) :
  assert (len(tables) > 0)
