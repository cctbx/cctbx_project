from __future__ import absolute_import, division, print_function
from libtbx.utils import null_out
import iotbx.pdb
from scitbx.array_family import flex

def two_group_selections_per_residue(model):
  model_ = model.deep_copy()
  if not model_.processed():
    model_.set_stop_for_unknowns(value=False)
    model_.set_log(log=null_out())
    model_.process(make_restraints=False)
  sc = model_.sel_sidechain().iselection()
  bb = model_.sel_backbone().iselection()
  result = []
  get_class = iotbx.pdb.common_residue_names_get_class
  for rg in model_.get_hierarchy().residue_groups():
    for ag in rg.atom_groups():
      g = ag.atoms().extract_i_seq()
      if(get_class(name = ag.resname) in ["common_rna_dna", "common_amino_acid"]):
        s1 = g.intersection(sc)
        s2 = g.intersection(bb)
        s3 = s1.intersection(s2)
        assert s3.size() == 0
        if(s1.size() > 0): result.append( s1 )
        if(s2.size() > 0): result.append( s2 )
      else:
        result.append(g)
  # Expensive sanity check
  tmp = flex.size_t()
  for g in result:
    tmp.extend(g)
  assert tmp.size() == model.size()
  #
  return result
