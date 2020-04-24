# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate


class Program(ProgramTemplate):
  description = '''
iotbx.split_data_cif: Tool to split cif file that has multiple data blocks

Usage examples:
  iotbx.split_data_cif 5r82-sf.cif
  '''
  def validate(self):
    print('Validating inputs:\n', file=self.logger)
    print(dir(self.data_manager), file=self.logger)
    print(self.data_manager.get_miller_array_labels(), file=self.logger)
    print(self.data_manager.get_miller_array_names(), file=self.logger)
    print(self.data_manager.get_miller_array_type(), file=self.logger)
    print(self.data_manager.get_miller_array_types(), file=self.logger)


  def run(self):
    pass

  def get_results(self):
    return None
