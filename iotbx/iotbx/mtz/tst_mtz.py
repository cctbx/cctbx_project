from iotbx import mtz
import urllib

def get_test_files(file):
  urllib.urlretrieve('http://cci.lbl.gov/build/'+file,file)

def print_uc(uc):
  uctup = tuple(uc)
  for n in uc:
    print "%10.4f"%n,
  print

def comprehensive_mtz(file):
  p = mtz.Mtz(file)    # Instantiate an MTZ file read object
  p.title()                  # title, as a string
  p.SpaceGroup()             # space group, as a string 
  p.UnitCell(0)              # unit cell of the first crystal as af_shared double
  p.ncolumns(0,0)            # column count of crystal 1, dataset 1
  p.columns()                # all column labels in the file, as a list of strings
  p.size()                   # number of reflections in file
  p.ncrystals()              # number of crystals
  p.history()                # history, as a list of strings
 #p.printHeader(3)           # Old-style header printout from ccp4
 #p.printHeaderAdv(4)        # New-style header printout from ccp4
  p.MIx()                    # copy of miller index list, as af_shared miller
  cryst = p.getCrystal(0)    # lightweight object representing the first crystal
  cryst.crystal_name()       # the crystal name as a string
  cryst.project_name()       # the project name as a string
  cryst.UnitCell()           # unit cell associated with this crystal
  cryst.ndatasets()          # number of data sets associated with this crystal
  data = cryst.getDataset(0) # lightweight object representing the first dataset
  data.dataset_name()        # the dataset name as a string
  data.wavelength()          # the dataset wavelength
  data.ncolumns()            # number of columns associated with this dataset
  col = data.getColumn(0)    # lightweight object representing the first column
  col.label()                # the label of this column
  col.type()                 # the ccp4 type of this column, as a one-letter code
  col[0]                     # reference to the first element in this column
  H = p.H                    # a copy of label "H" column in af_shared (double) format
  K = p.K                    #
  L = p.getShared("L")       # a copy of label "L" column in af_shared (double) format
  I = p.getColumn("I")       # lightweight object representing the column labelled "I"
  SIGI = p.getColumn("SIGI") #
  [H[0],K[0],L[0],I[0],SIGI[0]] # column elements can be referenced from either type

def comprehensive_mtz_2(file):
  p = mtz.Mtz(file)          # Instantiate an MTZ file read object
  m=p.MIx()                    # copy of miller index list, as af_shared (miller)
  h=p.HL("HLA","HLB","HLC","HLD") # copy of HL-coefficients, as af_shared (HL)  
  c=p.complex("FP","PHIB")     # complex structure factor, as af_shared (complex)
  #for x in xrange(p.size()):
  #  print m[x],h[x],c[x]
  
def exercise_mtzread(file):
  p = mtz.Mtz(file)
  print "Title:",p.title()
  print "Spacegroup:",p.SpaceGroup()
  print "Unitcell:",; print_uc(p.UnitCell(0))

  print "Columns:",
  p.columns()
  for c in p.columns():
    print "%s "%c,
  print

  print "Number of Reflections:",p.size()
  print "Number of Crystals:",p.ncrystals()
  print "History:",
  for line in p.history():
    print "  ",line
  for x in xrange(p.ncrystals()):
    print "Crystal",x+1
    cryst = p.getCrystal(x)
    print "  Name:",cryst.crystal_name()
    print "  Project:",cryst.project_name()
    print "  Unit Cell:",;print_uc(cryst.UnitCell())
    for d in xrange(cryst.ndatasets()):
      print "  Dataset",d+1,
      data=cryst.getDataset(d)
      print " columns:",data.ncolumns()
      print "    Name:",data.dataset_name()
      print "    Wavelength:",data.wavelength()
      for c in xrange(data.ncolumns()):
        col = data.getColumn(c)
        print "      Column%9d%20s%3s"%(c+1,col.label(),col.type())

def run():
  get_test_files("test.mtz")
  get_test_files("phase_noanom.mtz")
  comprehensive_mtz("test.mtz")
  exercise_mtzread("test.mtz")
  exercise_mtzread("phase_noanom.mtz")
  comprehensive_mtz_2("phase_noanom.mtz")

if __name__=="__main__":
  run()

