import mmtbx.chemical_components
from mmtbx.chemical_components import get_atom_names, get_hydrogen_names
from mmtbx.chemical_components import get_bond_pairs

def excercise(code):
  print '\nAtom names'
  print get_atom_names(code)
  print '\nAlternate atom names'
  print get_atom_names(code, alternate=True)
  print '\nHydrogen names'
  print get_hydrogen_names(code)
  print '\nAlternate hydrogen names'
  print get_hydrogen_names(code, alternate=True)
  print '\nWrapped hydrogen names'
  print get_hydrogen_names(code, wrap=True)
  print '\nBond pairs'
  print get_bond_pairs(code)
  print '\nAlternate name bond pairs'
  print get_bond_pairs(code, alternate=True)

def run():
  if (mmtbx.chemical_components.data_dir is None):
    print "Skipping tests: mmtbx.chemical_components.data_dir not available"
  else:
    for code in ["HOH", "ATP", "hem"]:
      print "\n%s\n%s" % ('_'*80,("%s " % code)*20)
      excercise(code)
  print "OK"

if __name__=="__main__":
  run()
