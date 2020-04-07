from __future__ import division, print_function
from rstbx.sublattice_support.change_basis import sublattice_change_of_basis

def example_application():
  S = sublattice_change_of_basis(6)
  n_items = 0
  for item in S.yield_transformations_descending_modulus():
    #item.show_summary()
    print (item.sublattice_cosets())
    n_items += 1
  assert n_items == 177

if __name__=="__main__":
  user_max_modulus = 6
  S = sublattice_change_of_basis(user_max_modulus)
  S.show_summary()
  example_application()
  print ("OK")
