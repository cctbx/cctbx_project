# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate
import mmtbx.model
import iotbx.pdb
import iotbx.phil
from scitbx.array_family import flex

from mmtbx.wwpdb.rcsb_web_services import reference_chain_search
from mmtbx.wwpdb.rcsb_entry_request import get_info
from iotbx.pdb.fetch import fetch

from libtbx.utils import Sorry, urlopen
from libtbx import Auto
import os

# =============================================================================

def fetch_model(model_id, model_info):
  """presently only experimental ones

  Args:
      model_id (_type_): _description_

  Raises:
      Sorry: _description_

  Returns:
      _type_: _description
  """
  try:
    if model_info.is_computational():
      url = model_info.get_source_url()
      m_data = urlopen(url)
    else:
      m_data = fetch(id=model_id, entity='model_cif')
  except RuntimeError as e:
    if str(e).find("Couldn't download") >= 0:
      pass
    else:
      raise e
  lines = m_data.read().splitlines()
  m = mmtbx.model.manager(
      model_input = iotbx.pdb.input(
          lines=flex.std_string(lines),
          source_info=None))
  return m

class Program(ProgramTemplate):
  program_name = 'find_reference'
  description = '''
mmtbx.find_reference: Tool to find reference models for supplied model

Usage examples:
  mmtbx.find_reference model.pdb
  '''

  datatypes = ['model', 'phil']
  data_manager_options = ['model_skip_expand_with_mtrix']

  master_phil_str = """\
    identity_cutoff = 0.9
      .type=float
    include_computed_models=False
      .type = bool
    output_phil_format = *refine rsr
      .type = choice
"""

  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
    model = self.data_manager.get_model()
    self.data_manager.has_models(expected_n=1, exact_count=True, raise_sorry=True)
    if model.get_hierarchy().models_size() != 1:
      raise Sorry("Multi-model files are not supported.")
    print ('Inputs OK', file=self.logger)

  # ---------------------------------------------------------------------------

  def write_homologue_summary(self, model_infos, first_n=10):
    print("    pdb id, resolution, rwork, rfree, first %d items" % first_n, file=self.logger)
    print("    pdb id, pLDDT, first %d items" % first_n, file=self.logger)
    for mi in model_infos[:first_n]:
      if not mi.is_computational():
        print("    %s: %s, %s, %s" % (
            mi.get_pdb_id(),
            mi.get_resolution(),
            mi.get_rwork(),
            mi.get_rfree()), file=self.logger)
      else:
        print("    %s: %.2f" % (
            mi.get_pdb_id(),
            mi.get_plddt()), file=self.logger)

  def pick_homologue(self, search_results):
    """pick the best homologue
    """
    model_infos = get_info([a.split('.')[0].replace(":", "_") for a in search_results])
    self.write_homologue_summary(model_infos, first_n=10)
    # try to return first computational instead:
    # for sr, mi in zip(search_results, model_infos):
    #   if mi.is_computational():
    #     return sr, mi
    return search_results[0], model_infos[0]

  def run(self):
    self.result = []
    model = self.data_manager.get_model()
    for chain in model.get_hierarchy().only_model().chains():
      # find homology
      print("Working with chain %s" % chain.id, file=self.logger)
      if not chain.is_protein():
        continue
      seq = chain.as_padded_sequence()
      ref_search_result = reference_chain_search(sequence=seq,
                             identity_cutoff=self.params.identity_cutoff,
                             include_csm=self.params.include_computed_models)
      print('  Chain sequence:', seq)
      print('  Reference search results:',ref_search_result, file=self.logger)
      homology_model, model_info = self.pick_homologue(ref_search_result)
      print('  Picked homology model: %s' % homology_model, file=self.logger)
      if homology_model:
        ref_model_id, ref_label_id = homology_model.split('.')
        ref_m = fetch_model(ref_model_id, model_info)
        ref_m.add_crystal_symmetry_if_necessary()
        label_auth_dict = ref_m.get_model_input().label_to_auth_asym_id_dictionary()
        ref_chain_id = label_auth_dict.get(ref_label_id, ref_label_id)
        # print("  ref label --> chain: '%s' --> '%s' " % (ref_label_id, ref_chain_id), file=self.logger)
        # print('  ref_m atom size: ', ref_m.get_hierarchy().atoms_size())
        selection = ref_m.selection("chain %s and protein" % ref_chain_id)
        cutted_ref_m = ref_m.select(selection)
        atoms_removed = cutted_ref_m.remove_alternative_conformations(always_keep_one_conformer=True)
        # print('  alt conf atoms_removed', atoms_removed)
        cutted_ref_m = cutted_ref_m.remove_hydrogens()
        print('  Final atom size of reference chain: ', cutted_ref_m.get_hierarchy().atoms_size(), file=self.logger)
        fname = 'reference_for_chain_%s_%s.%s' % (chain.id, ref_model_id, ref_chain_id)
        ref_fname = self.data_manager.write_model_file(
            model_str=cutted_ref_m,
            filename=fname,
            overwrite=True)
        chain_result = (chain.id, ref_chain_id, ref_fname)
        print("  Result appended with: ", chain_result, file=self.logger)
        self.result.append(chain_result)
    print("SELF RESULT", self.result, file=self.logger)
    print(self.get_reference_model_param_file_content(output_phil_format=self.params.output_phil_format))

    fn = "%s.eff" % self.get_default_output_filename()
    if not self.params.output.prefix and not self.params.output.suffix:
      fn = "%s.eff" % self.get_default_output_filename(
          prefix='%s_references' % (os.path.basename(self.data_manager.get_default_model_name())[:-4]),
          serial=Auto)
    print("Writing parameters to: ", fn, file=self.logger)
    self.data_manager.write_phil_file(
      self.get_reference_model_param_file_content(output_phil_format=self.params.output_phil_format),
      filename=fn)

  def get_reference_model_param_file_content(self, output_phil_format='refine'):
    result = ""
    if not self.result:
      return result
    files_text = "\n".join("  file=%s" %i[2] for i in self.result)
    groups_text = ""
    for item in self.result:
      groups_text += """
  reference_group {
    reference = chain %s
    selection = chain %s
    file_name = %s
  }""" % (item[1], item[0], item[2])
    result = ''
    if output_phil_format == 'refine':
      result = 'refinement.'
    result +="""\
reference_model {
  enabled=True
%s
%s
}
""" % (files_text, groups_text)
    return result

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.result

