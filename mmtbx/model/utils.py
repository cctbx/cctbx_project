from __future__ import absolute_import, division, print_function

import iotbx.pdb.hierarchy
import mmtbx.model

def merge_models(models):
  """
  Merge models together.

  Models must have same number of hierarchy's 'models' there must be no duplicate
  chain/residue combinations. Models must not contain geometry restraints manager

  Args:
      models (iterable): List or tuple with models to be merged.
        They are not changed inside the function.

  Returns new model object.
  """
  hierarchies = []
  h_models_size = []
  for i, m in enumerate(models):
    h_errors = m.get_hierarchy().overall_counts().errors()
    assert len(h_errors) == 0,\
      "Model #%d has erros: %s" % (i+1, h_errors)
    hierarchies.append(m.get_hierarchy())
    h_models_size.append(m.get_hierarchy().models_size())
    assert m.restraints_manager is None, "Model #%d has GRM initialized." % (i+1)
  assert len(set(h_models_size)) == 1,\
    "Cannot merge models that have differing numbers of hierarchy 'models'"
  for i, m in enumerate(models[1:]):
    assert models[0].crystal_symmetry().is_identical_symmetry(m.crystal_symmetry()),\
      "Cannot merge models with different crystal symmetries (#%d)" % (i+2)
  n_h_models = set(h_models_size).pop()
  # check chain ids are different
  for i in range(n_h_models):
    cids = []
    for h in hierarchies:
      cids += [c.id for c in h.models()[i].chains()]
    assert len(cids) == len(set(cids)), "Duplicated chain ids are present."

  # Merge models
  root = iotbx.pdb.hierarchy.root()
  root.pre_allocate_models(n_h_models)
  for i in range(n_h_models):
    new_model = iotbx.pdb.hierarchy.model()
    new_model.id = '%d' % (i+1)
    for h in hierarchies:
      m_copy = h.models()[i].detached_copy()
      new_model.transfer_chains_from_other(m_copy)
    root.append_model(new_model)
  root.atoms().reset_i_seq()
  root.atoms().reset_serial()
  oc = root.overall_counts()
  assert len(oc.errors()) == 0, oc.errors()
  # should we check warnings as well?
  # print(oc.warnings())
  return mmtbx.model.manager(
      model_input=None,
      pdb_hierarchy=root,
      crystal_symmetry=models[0].crystal_symmetry())
