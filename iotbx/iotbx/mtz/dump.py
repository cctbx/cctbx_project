from iotbx import mtz
import sys

def dump(file_name):
  p = mtz.Mtz(file_name)
  print "Title:", p.title()
  print "Space group symbol:", p.SpaceGroup()
  print "Number of crystals:", p.ncrystals()
  print "Number of Miller indices:", p.size()
  print "History:"
  for line in p.history():
    print " ", line
  for i_crystal in xrange(p.ncrystals()):
    print "Crystal %d:" % (i_crystal+1)
    cryst = p.getCrystal(i_crystal)
    print "  Name:", cryst.crystal_name()
    print "  Project:", cryst.project_name()
    print " ",
    cryst.UnitCell().show_parameters()
    print "  Number of datasets:", cryst.ndatasets()
    for i_dataset in xrange(cryst.ndatasets()):
      print "  Dataset %d:" % (i_dataset+1)
      dataset = cryst.getDataset(i_dataset)
      print "    Name:", dataset.dataset_name()
      print "    Wavelength: %.6g" % dataset.wavelength()
      print "    Number of columns:", dataset.ncolumns()
      print "    Column number, label, observations, type:"
      for i_column in xrange(dataset.ncolumns()):
        column = dataset.getColumn(i_column)
        n_oversations = p.valid_indices(column.label()).size()
        print "      %3d, %s, %d/%d=%.2f%%, %s: %s" % (
          i_column+1, column.label(),
          n_oversations, p.size(), 100.*n_oversations/p.size(),
          column.type(), mtz.column_type_legend[column.type()])
