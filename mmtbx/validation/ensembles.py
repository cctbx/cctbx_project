
# TODO tests

from __future__ import absolute_import, division, print_function
from mmtbx.validation import ramalyze
from mmtbx.validation import rotalyze
from libtbx import slots_getstate_setstate_default_initializer
from libtbx.utils import null_out
from libtbx import easy_mp
from libtbx import Auto
import sys
from six.moves import range

class ensemble_validation(object):
  def __init__(self,
      pdb_hierarchies, # XXX these must be single-conformer
      reference_hierarchy=None,
      nproc=Auto,
      log=sys.stdout):
    self.pdb_hierarchies = pdb_hierarchies
    self.selection_caches = [h.atom_selection_cache() for h in pdb_hierarchies]
    if (reference_hierarchy == None):
      reference_hierarchy = self.pdb_hierarchies[0]
    self.residue_ids = []
    self.residue_id_dict = {}
    for chain in reference_hierarchy.only_model().chains():
      if (not chain.is_protein()) : # TODO
        continue
      for residue_group in chain.residue_groups():
        residue = residue_group.only_atom_group()
        id_str = residue.id_str()
        self.residue_ids.append(id_str)
        self.residue_id_dict[id_str] = len(self.residue_ids) - 1
    self.residue_ensembles = [ list([]) for id_str in self.residue_ids ]
    for hierarchy in pdb_hierarchies :
      for chain in hierarchy.only_model().chains():
        if (not chain.is_protein()):
          continue
        for residue_group in chain.residue_groups():
          residue = residue_group.only_atom_group()
          id_str = residue.id_str()
          i_res = self.residue_id_dict.get(id_str)
          assert (i_res is not None)
          self.residue_ensembles[i_res].append(residue)
    self.validations = easy_mp.pool_map(
      fixed_func=self.validate_single_model,
      iterable=range(len(pdb_hierarchies)),
      processes=nproc)
    rama_by_residue = combine_model_validation_results(
      validation_objects=[ rama for rama, rota in self.validations ],
      ensemble_result_class=ramalyze.ramachandran_ensemble,
      residue_ids=self.residue_ids,
      ignore_unexpected_residues=False,
      log=log)
    rota_by_residue = combine_model_validation_results(
      validation_objects=[ rota for rama, rota in self.validations ],
      ensemble_result_class=rotalyze.rotamer_ensemble,
      residue_ids=self.residue_ids,
      ignore_unexpected_residues=False,
      log=log)
    assert len(rama_by_residue)==len(rota_by_residue)==len(self.residue_ids)
    self.residue_data = []
    for i_res, id_str in enumerate(self.residue_ids):
      residues = self.residue_ensembles[i_res]
      self.residue_data.append(
        residue_analysis(
          id_str=id_str,
          residues=residues,
          rama=rama_by_residue[i_res],
          rota=rota_by_residue[i_res]))

  def validate_single_model(self, i_model):
    hierarchy = self.pdb_hierarchies[i_model]
    rama_validation = ramalyze.ramalyze(
      pdb_hierarchy=hierarchy,
      outliers_only=False)
    rota_validation = rotalyze.rotalyze(
      pdb_hierarchy=hierarchy,
      outliers_only=False)
    return rama_validation, rota_validation

  def show(self, out=sys.stdout, prefix=""):
    for residue in self.residue_data :
      residue.show(out=out, prefix=prefix)

  def get_validated_residues_in_selection(self,
      selection,
      require_n_residues=None,
      log=None):
    if (log is None) : log = null_out()
    results = []
    i_model = 0
    while (i_model < len(self.pdb_hierarchies)):
      hierarchy = self.pdb_hierarchies[i_model]
      sel_cache = self.selection_caches[i_model]
      rama, rota = self.validations[i_model]
      i_model += 1
      isel = sel_cache.selection(selection).iselection()
      hierarchy_sel = hierarchy.select(isel)
      residue_groups = hierarchy_sel.only_model().only_chain().residue_groups()
      if (require_n_residues is not None):
        if (len(residue_groups) != require_n_residues):
          results.append(None)
          continue
      reject = False
      for residue_group in residue_groups :
        atom_group = residue_group.only_atom_group()
        rama_result = rama.find_atom_group(other=atom_group)
        rota_result = rota.find_atom_group(other=atom_group)
        assert (not None in [rama_result, rota_result]), atom_group.id_str()
        if (rama_result.is_outlier()) or (rota_result.is_outlier()):
          reject = True
          break
      if (reject):
        results.append(None)
      else :
        results.append(hierarchy_sel)
    return results

class residue_analysis(slots_getstate_setstate_default_initializer):
  __slots__ = ["id_str", "rama", "rota", "residues"]
  def show(self, out=sys.stdout, prefix=""):
    if (self.rota is not None):
      print(prefix + self.rota.as_string(), file=out)

def combine_model_validation_results(
    validation_objects,
    ensemble_result_class,
    residue_ids=None,
    get_key_callback=None,
    ignore_unexpected_residues=True,
    log=sys.stderr):
  def default_get_key_callback(result, i_model):
    return result.atom_group_id_str()
  if (get_key_callback is None):
    get_key_callback = default_get_key_callback
  else :
    assert hasattr(get_key_callback, "__call__")
  if (len(residue_ids) is None):
    residue_ids = []
    for result in validation_objects[0].results :
      residue_ids.append(get_key_callback(result=result, i_model=0))
  results_dict = dict([ (id_str, list([])) for id_str in residue_ids ])
  for i_model, validation_object in enumerate(validation_objects):
    for result in validation_object.results :
      id_str = get_key_callback(result=result, i_model=i_model)
      if (not id_str in results_dict):
        msg = "The residue ID '%s' was not found in the expected list."%id_str
        if (ignore_unexpected_residues):
          print("  " + msg, file=log)
        else :
          raise RuntimeError(msg)
      else :
        results_dict[id_str].append(result)
  results = []
  for id_str in residue_ids :
    rr = results_dict[id_str]
    if (len(rr) == 0):
      results.append(None)
    else :
      results.append(ensemble_result_class(rr))
  return results
