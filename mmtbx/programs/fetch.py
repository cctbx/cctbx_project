"""Automatically retrieve data from the PDB via the RCSB \
    web server"""
from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate
from libtbx import easy_run
from libtbx.utils import Sorry
from iotbx.pdb.fetch import valid_pdb_id, fetch_and_write
from mmtbx.wwpdb import rcsb_web_services
import os

master_phil_str = """
fetch
  .caption = Automatically retrieve data from the PDB via the RCSB \
    web server.  If you intend to re-refine or re-build the structure we \
    recommend creating a new project, but this is not required.  Note that \
    you may also use this tool to generate an MTZ file from the mmCIF \
    structure factors (if available), but the options \
    are more limited than what is available in the phenix.cif_as_mtz.
  .style = auto_align caption_img:icons/custom/pdb_import64.png \
    caption_width:400
{
  pdb_ids = None
    .type = strings
    .short_caption = PDB ID(s)
    .input_size = 400
    .style = bold
  action = *model data sequence half_maps all
    .type = choice(multi=True)
    .caption = model_file(s) data_file(s) sequence half_maps
  convert_to_mtz = False
    .type = bool
    .caption = Try to convert X-ray data to mtz format
  mirror = *rcsb pdbe pdbj
    .type = choice
    .caption = RCSB PDBe PDBj
    .short_caption = Mirror site
    .style = bold
}
"""

class Program(ProgramTemplate):
  description = """
  Fetch model, data, sequence files. Optionally convert to sf data to mtz format.
  Usage:
  Get only model:
    iotbx.fetch_pdb 1yjp
  Get model and data(xray or cryo-em):
    iotbx.fetch_pdb 1yjp action=model+data
  Get everything:
    iotbx.fetch_pdb 1yjp action=all
"""
  datatypes = ['phil']
  master_phil_str = master_phil_str

  def custom_init(self):
    # store output filenames and error messages
    self.output_filenames = []
    self.errors = []

  def validate(self):
    print('Validating inputs:\n', file=self.logger)
    if self.params.fetch.pdb_ids is  None:
      raise Sorry("Provide pdb id for fetching")
    for pdb_id in self.params.fetch.pdb_ids:
      if not valid_pdb_id(pdb_id):
        raise Sorry("Invalid PDB code: %s" % pdb_id)
    for a in self.params.fetch.action:
      if a not in ['model', 'data', 'sequence', 'half_maps', 'all']:
        raise Sorry("Unsupported action %s" % a)

  def define_entities_to_fetch(self, emdb_number):
    entities_to_fetch = []
    if 'model' in self.params.fetch.action:
      entities_to_fetch += ['model_pdb', 'model_cif']
    if 'data' in self.params.fetch.action:
      if emdb_number is None:
        entities_to_fetch += ['sf']
      else:
        entities_to_fetch += ['em_map']
    if 'sequence' in self.params.fetch.action:
      entities_to_fetch += ['sequence']
    if 'half_maps' in self.params.fetch.action and emdb_number:
      entities_to_fetch += ['em_half_map_1', 'em_half_map_2']

    if 'all' in self.params.fetch.action:
      entities_to_fetch = ['model_pdb', 'model_cif', 'sequence']
      if emdb_number is not None:
        entities_to_fetch += ['em_map','em_half_map_1', 'em_half_map_2']
      else:
        entities_to_fetch += ['sf']
    return entities_to_fetch

  def run(self):
    for pdb_id in self.params.fetch.pdb_ids:
      emdb_number = None
      emdb_ids = rcsb_web_services.get_emdb_id_for_pdb_id(pdb_id)
      if emdb_ids is not None:
        emdb_number = emdb_ids[0].split('-')[1]
      print("Fetching: PDB ID: %s, EMDB ID: %s" % (pdb_id, emdb_number), file=self.logger)
      entities_to_fetch = self.define_entities_to_fetch(emdb_number)
      for e in entities_to_fetch:
        fn = fetch_and_write(
            id=pdb_id,
            entity=e,
            mirror=self.params.fetch.mirror,
            emdb_number=emdb_number,
            log=self.logger)
        if e == 'sf' and fn is not None and self.params.fetch.convert_to_mtz:
          # here we have successfully have structure factor file and
          # need to convert it to mtz.
          cmd = "phenix.cif_as_mtz %s --merge --map_to_asu --extend_flags --ignore_bad_sigmas" % fn
          easy_run.call(cmd)
          if os.path.isfile("%s-sf.mtz" % pdb_id):
            new_filename = "%s.mtz" % pdb_id
            os.rename("%s-sf.mtz" % pdb_id, new_filename)
            print("Converted structure factors saved to %s.mtz" % pdb_id, file=self.logger)
            self.output_filenames.append(new_filename)
          if (not os.path.isfile("%s.mtz" % pdb_id)):
            error_msg = "MTZ conversion failed - try running phenix.cif_as_mtz " \
              + "manually (and check %s-sf.cif for format errors)." % pdb_id
            print(error_msg, file=self.logger)
            self.errors.append(error_msg)
        elif fn is not None:
          self.output_filenames.append(fn)

  def get_results(self):
    return self.output_filenames, self.errors
