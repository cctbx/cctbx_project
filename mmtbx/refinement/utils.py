from __future__ import absolute_import, division, print_function
from libtbx.utils import null_out

def two_group_selections_per_residue(
      model):
      # all_chain_proxies=None,
      # pdb_hierarchy=None,
      # sidechain_selection=None,
      # backbone_selection=None):
  # if(all_chain_proxies is None):
  #   assert [pdb_hierarchy,sidechain_selection,backbone_selection].count(None)==0
  #   sc = sidechain_selection
  #   bb = backbone_selection
  # else:
  #   assert [pdb_hierarchy,sidechain_selection,backbone_selection].count(None)==3
  #   pdb_hierarchy = all_chain_proxies.pdb_hierarchy
  #   sc = all_chain_proxies.sel_sidechain().iselection()
  #   bb = all_chain_proxies.sel_backbone().iselection()
  model_ = model.deep_copy()
  if not model_.processed():
    model_.set_stop_for_unknowns(value=False)
    model_.set_log(log=null_out())
    model_.process(make_restraints=False)
  sc = model_.sel_sidechain().iselection()
  bb = model_.sel_backbone().iselection()
  ogpr = [rg.atoms().extract_i_seq() for rg in model_.get_hierarchy().residue_groups()]
  result = []
  for g in ogpr:
    s1 = g.intersection(sc)
    s2 = g.intersection(bb)
    s3 = s1.intersection(s2)
    assert s3.size() == 0
    if(s1.size() > 0):
      result.append( s1 )
    if(s2.size() > 0):
      result.append( s2 )
  return result
