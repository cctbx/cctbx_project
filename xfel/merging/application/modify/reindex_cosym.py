from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from cctbx.sgtbx import change_of_basis_op
from cctbx import miller

"""The sole purpose here is to reindex a reflection table, once the coset assignments
   for each experiment have been determined"""

def reindex_refl_by_coset(refl,data,symms,uuids,co,anomalous_flag = False, verbose=True):
  if verbose:
    cache_miller = refl["miller_index"].deep_copy()
    cache_asu = refl["miller_index_asymmetric"].deep_copy()

  for icoset,partition in enumerate(co.partitions):
    if icoset==0: continue # no change of basis

    cb_op = change_of_basis_op(partition[0])
    mi_new = cb_op.apply(refl["miller_index"])
    mi_asu_new = mi_new.deep_copy()
    miller.map_to_asu(symms[0].space_group().info().type(), anomalous_flag, mi_asu_new)

    # now select only those expts assigned to that coset
    good_refls = flex.bool(len(refl), False)
    all_expt_id = list(data["experiment"])
    all_coset = list(data["coset"]) # would like to understand how to use pandas rather than Python list
    for iexpt in range(len(symms)):
        iexpt_id = uuids[iexpt]
        this_coset = all_coset[ all_expt_id.index(iexpt_id) ]
        if this_coset == icoset:
          good_refls |= refl["exp_id"] == iexpt_id

    re_millers = mi_new.select(good_refls)
    re_asu = mi_asu_new.select(good_refls)

    refl["miller_index"].set_selected(good_refls, re_millers)
    refl["miller_index_asymmetric"].set_selected(good_refls, re_asu)
  if verbose:
    for x in range(len(refl)):
      print ("%3d"%x,refl["exp_id"][x],
             all_coset[ all_expt_id.index(refl["exp_id"][x]) ],
             "%10s"%(change_of_basis_op(co.partitions[   all_coset[ all_expt_id.index(refl["exp_id"][x]) ]   ][0]).as_hkl()),
             "(%4d%4d%4d)"%(cache_miller[x]), "(%4d%4d%4d)"%(cache_asu[x]),
             "(%4d%4d%4d)"%(refl["miller_index"][x]), "(%4d%4d%4d)"%(refl["miller_index_asymmetric"][x]))

  return refl

if __name__=="__main__":

  import pickle
  with open("refl.pickle","rb") as F:
    refl = pickle.load(F)
    print(refl)
    good_refls = flex.bool(len(refl), False)
    for x in range(0,len(refl),500):
      good_refls[x]=True
    refl = refl.select(good_refls)
    data = pickle.load(F)
    print(data)
    symms = pickle.load(F)
    print(symms)
    uuids = pickle.load(F)
    co = pickle.load(F)
    print(co)
    co.show()
    reindex_refl_by_coset(refl,data,symms,uuids,co)
