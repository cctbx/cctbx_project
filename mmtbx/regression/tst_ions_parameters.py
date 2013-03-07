
from __future__ import division

def exercise () :
  from mmtbx.ions.parameters import server, CHARGES
  import iotbx.pdb.hierarchy
  s = server()
  for elem in ["NA","MG","K","CA","MN","FE","NI","CO","CU","ZN","CD"] :
    params = s.get_metal_parameters(elem)
    assert (len(params.allowed_coordinating_atoms) > 0)
    assert (params.charge is not None)
    for elem2 in params.allowed_coordinating_atoms :
      atom1 = iotbx.pdb.hierarchy.atom()
      atom1.name = elem
      atom1.element = elem
      atom1.charge = "+" + str(params.charge)
      atom2 = iotbx.pdb.hierarchy.atom()
      atom2.name = elem2
      atom2.element = elem2
      charge2 = str(CHARGES[elem2])
      if (not "-" in charge2) :
        charge2 = "+" + charge2
      atom2.charge = charge2
      assert (s.get_valence_params(atom1, atom2) != (None, None))
  print "OK"

if (__name__ == "__main__") :
  exercise()
