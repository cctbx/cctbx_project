from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from mmtbx.hydrogens.reduce_hydrogen import get_ligand_interactions
from mmtbx.hydrogens.reduce_hydrogen import get_incorrect_hydrogens_for_bond

# Kept on request. The call in reduce_hydrogen is commented out: linking, bond
# orders and the HIS exception handle these cases now, and the tests moved to
# tst_add_hydrogen_12. This file is not in run_tests.py and has no tests.

def workaround_003(model):
  pairs = get_ligand_interactions(
    model=model, dist_min=1.1, cutoff_cno=1.6, cutoff_sp=1.9)
  remove_selection = []
  for pair in pairs:
    badH_i_seqs = get_incorrect_hydrogens_for_bond(
      model   = model,
      i_seq_A = pair[0],
      i_seq_B = pair[1])
    if len(badH_i_seqs)>0:
      remove_selection.extend(badH_i_seqs)
  removed = 0
  if len(remove_selection)>0:
    # NOTE: removes badH_i_seqs (last pair only), not remove_selection.
    badH_i_seqs = flex.size_t(badH_i_seqs)
    removed = badH_i_seqs.size()
    keep_selection = ~flex.bool(model.size(), badH_i_seqs)
    model = model.select(keep_selection)
  return model, removed
