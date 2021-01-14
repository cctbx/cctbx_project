from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from cctbx.miller import set as miller_set
from cctbx.sgtbx import change_of_basis_op

"""The sole purpose here is to reindex a reflection table, once the coset assignments
   for each experiment have been determined"""

def reindex_refl_by_coset(refl,data,symms,uuids,co,verbose=False):
  if verbose:
    cache_miller = refl["miller_index"].deep_copy()
    cache_asu = refl["miller_index_asymmetric"].deep_copy()

  for icoset,partition in enumerate(co.partitions):
    if icoset==0: continue # no change of basis
    print("Coset",icoset)
    #from IPython import embed; embed()
    MI = miller_set(crystal_symmetry=symms[0], indices = refl["miller_index"], anomalous_flag=False)
    MI_new = MI.change_basis( change_of_basis_op(partition[0]) )
    MIA_new = MI_new.map_to_asu()
    #MIA_new = MIA_new.as_non_anomalous_set() # highly questionable.  Is this really what we want?

    # now select only those expts assigned to that coset
    good_refls = flex.bool(len(refl), False)
    all_expt_id = list(data["experiment"])
    all_coset = list(data["coset"]) # would like to understand how to use pandas rather than Python list
    for iexpt in range(len(symms)):
        iexpt_id = uuids[iexpt]
        this_coset = all_coset[ all_expt_id.index(iexpt_id) ]
        if this_coset == icoset:
          good_refls |= refl["exp_id"] == iexpt_id

    re_millers = MI_new.indices().select(good_refls)
    re_asu = MIA_new.indices().select(good_refls)

    refl["miller_index"].set_selected(good_refls, re_millers)
    refl["miller_index_asymmetric"].set_selected(good_refls, re_asu)
  if verbose:
    for x in range(len(refl)):
      print (x,refl["exp_id"][x],
             all_coset[ all_expt_id.index(refl["exp_id"][x]) ],
             "(%4d%4d%4d)"%(cache_miller[x]), "(%4d%4d%4d)"%(cache_asu[x]),
             "(%4d%4d%4d)"%(refl["miller_index"][x]), "(%4d%4d%4d)"%(refl["miller_index_asymmetric"][x]))

  return refl

# map_to_asu
if __name__=="__main__":

  import pickle
  with open("refl.pickle","rb") as F:
    refl = pickle.load(F)
    print(refl)
    data = pickle.load(F)
    print(data)
    symms = pickle.load(F)
    print(symms)
    uuids = pickle.load(F)
    co = pickle.load(F)
    print(co)
    reindex_refl_by_coset(refl,data,symms,uuids,co)
