from __future__ import generators
import cctbx.array_family.flex

import boost.python
ext = boost.python.import_ext("iotbx_mtz_wrapper_ext")
from iotbx_mtz_wrapper_ext import *
import iotbx_mtz_wrapper_ext as ext

from cctbx import sgtbx
import sys

column_type_legend_source = \
  "http://www.ccp4.ac.uk/dist/html/mtzlib.html#fileformat"
column_type_legend = {
  "H": "index h,k,l",
  "J": "intensity",
  "F": "amplitude",
  "D": "anomalous difference",
  "Q": "standard deviation",
  "G": "F(+) or F(-)",
  "L": "standard deviation",
  "K": "I(+) or I(-)",
  "M": "standard deviation",
  "P": "phase angle in degrees",
  "W": "weight (of some sort)",
  "A": "phase probability coefficients (Hendrickson/Lattmann)",
  "B": "BATCH number",
  "Y": "M/ISYM, packed partial/reject flag and symmetry number",
  "I": "integer",
  "R": "real",
}

class _object(boost.python.injector, ext.object):

  def space_group_info(self):
    return sgtbx.space_group_info(group=self.space_group())

  def columns(self):
    for crystal in self.crystals():
      for dataset in crystal.datasets():
        for column in dataset.columns():
          yield column

  def column_labels(self):
    return [column.label() for column in self.columns()]

  def column_types(self):
    return [column.type() for column in self.columns()]

  def show_summary(self, out=None):
    if (out is None): out = sys.stdout
    print >> out, "Title:", self.title()
    print >> out, "Space group symbol from file:", self.space_group_name()
    self.space_group_info().show_summary(
      f=out, prefix="Space group from matrices: ")
    print >> out, "Number of crystals:", self.n_crystals()
    print >> out, "Number of Miller indices:", self.n_reflections()
    print >> out, "History:"
    for line in self.history():
      print >> out, " ", line.rstrip()
    for i_crystal,crystal in enumerate(self.crystals()):
      print >> out, "Crystal %d:" % (i_crystal+1)
      print >> out, "  Name:", crystal.name()
      print >> out, "  Project:", crystal.project_name()
      crystal.unit_cell().show_parameters(f=out, prefix="  Unit cell: ")
      print >> out, "  Number of datasets:", crystal.n_datasets()
      for i_dataset,dataset in enumerate(crystal.datasets()):
        print >> out, "  Dataset %d:" % (i_dataset+1)
        print >> out, "    Name:", dataset.name()
        print >> out, "    Id:", dataset.id()
        print >> out, "    Wavelength: %.6g" % dataset.wavelength()
        print >> out, "    Number of columns:", dataset.n_columns()
        print >> out, "    Column number, label, observations, type:"
        for i_column,column in enumerate(dataset.columns()):
          n_observations = column.valid_indices().size()
          print >> out, "      %3d, %s, %d/%d=%.2f%%, %s: %s" % (
            i_column+1,
            column.label(),
            n_observations,
            self.n_reflections(),
            100.*n_observations/self.n_reflections(),
            column.type(),
            column_type_legend[column.type()])
