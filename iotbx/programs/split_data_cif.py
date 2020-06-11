# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate

from mmtbx.command_line.cif_as_mtz import extract
import iotbx.cif
import iotbx.cif.model
from cctbx import crystal
from libtbx import Auto
from libtbx.utils import Sorry
import os

class Program(ProgramTemplate):
  description = '''
iotbx.split_data_cif: Tool to split cif file that has multiple data blocks

Usage examples:
  iotbx.split_data_cif 5r82-sf.cif
  '''

  datatypes = ['phil', 'miller_array']

  master_phil_str = """\
output_r_free_label = 'R-free-flags'
  .type = str
  .help = MTZ column label to use for R-free flags (default: R-free-flags)
merge_non_unique_under_symmetry = True
  .type = bool
  .help = Merge non-unique data where present
incompatible_flags_to_work_set = True
  .type = bool
  .help = When merging place reflections with incompatible flags into the \
          working set
remove_systematic_absences = True
  .type = bool
  .help = Remove systematic absent reflections
map_to_asu = True
  .type = bool
  .help = Map to asymmetric unit
show_details_if_error = True
  .type = bool
  .help = Show data details for some errors
ignore_bad_sigmas = True
  .type = bool
  .help = Set sigmas to None instead of raising an error when bad sigmas \
          are present
extend_flags = True
  .type = bool
  .help = Extend R-free flags to cover all reflections if necessary
output {
  mtz = True
  cif = True
}
"""

  def validate(self):
    print('Validating inputs:\n', file=self.logger)
    self.data_manager.has_miller_arrays(raise_sorry=True)

  def run(self):
    self.out_cif_names = []
    self.out_mtz_names = []
    fnames = self.data_manager.get_miller_array_names()
    for fname in fnames:
      cif_reader = self.data_manager.get_miller_array(filename=fname).file_content()
      cif_model = cif_reader.model()
      dummy_cs = crystal.symmetry(
            unit_cell=None,
            space_group_info=None)
      for data_block_name in cif_model.keys():
        fn_base = self.get_default_output_filename(
            prefix='%s_%s' % (fname, data_block_name),
            serial=Auto)
        if self.params.output.cif:
          fn = fn_base+'.cif'
          if os.path.isfile(fn) and not self.params.output.overwrite:
            raise Sorry("%s already exists and overwrite is set to False." % fn)
          print ("Writing '%s'" % fn, file=self.logger)
          output_model = iotbx.cif.model.cif(blocks={data_block_name:cif_model[data_block_name]})
          with open(fn, 'w') as f:
            output_model.show(out=f)
          self.out_cif_names.append(fn)
        if self.params.output.mtz:
          miller_arrays = cif_reader.build_miller_arrays(data_block_name=data_block_name)
          mtz_object = extract(
              file_name='',
              crystal_symmetry=dummy_cs,
              wavelength_id=None,
              crystal_id=None,
              show_details_if_error=self.params.show_details_if_error,
              output_r_free_label=self.params.output_r_free_label,
              merge_non_unique_under_symmetry=self.params.merge_non_unique_under_symmetry,
              map_to_asu=self.params.map_to_asu,
              remove_systematic_absences=self.params.remove_systematic_absences,
              all_miller_arrays={data_block_name:miller_arrays},
              incompatible_flags_to_work_set=self.params.incompatible_flags_to_work_set,
              ignore_bad_sigmas=self.params.ignore_bad_sigmas,
              extend_flags=self.params.extend_flags,
              return_as_miller_arrays=False,
              log=self.logger)
          fn = fn_base+'.mtz'
          print ("Writing '%s'" % fn, file=self.logger)
          self.data_manager.write_miller_array_file(mtz_object, filename=fn)
          self.out_mtz_names.append(fn)

  def get_results(self):
    return self.out_cif_names, self.out_mtz_names
