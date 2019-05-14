#
# Here we define one common way of constructing command-line tools which is
# defined here, so everybody could finally forget how to parse command-line,
# read model file etc.
# They should have common (consistent):
# - way of parsing command-line
# - output help
# - reading files, check and react for e.g. IO errors, wrong formats, etc.
#
# This file will end up somewhere like libtbx/command_line/base_cl_app.py
#

from __future__ import absolute_import, division, print_function

import sys
import iotbx.phil
import iotbx.pdb
import six

input_model_phil_str = """\
input_model_fname = None
  .type = path
"""
input_model_phil = iotbx.phil.parse(input_model_phil_str)


class cl_app(object):
  def __init__(self,
      cl_args):
    self.cl_args = cl_args
    self.new_master_par_str = self.master_params_str
    self.pdbf_def = None
    for k, v in six.iteritems(self.needed_info):
      setattr(self, v, None)
      # adjust params_str for files
      if k == "model":
        self.new_master_par_str = input_model_phil_str + self.new_master_par_str
        self.pdbf_def = "input_model_fname"
      # other options to follow:
      if k == "xray_data":
        pass

  def parse_cl(self):
    print("In parent, running parse_cl")
    print("  ", self.needed_info)
    # processing and parsing cl_args
    self.new_master_par = iotbx.phil.parse(self.new_master_par_str, process_includes=True)
    self.params = self.new_master_par.extract()
    if self.log is None:
      self.log = sys.stdout

    # this should go to command-line processing class, but for now:
    if len(self.cl_args) <= 1 or "--help" in self.cl_args or "-h" in self.cl_args:
      self.show_help()
      return 1
    input_objects = iotbx.phil.process_command_line_with_files(
      args=self.cl_args[1:],
      master_phil=self.new_master_par,
      pdb_file_def=self.pdbf_def,
      )
    self.work_params = input_objects.work.extract()
    # If everything looks fine:
    return 0

  def show_help(self):
    print(self.help_message, file=self.log)
    self.new_master_par.show(self.log)

  def read_and_validate_inputs(self):
    # some generic validation, read model and data
    print("In parent, running validate inputs")
    if self.pdbf_def is not None:
      self.read_model_file()
    return 0

  def read_model_file(self):
    pdb_inp = iotbx.pdb.input(file_name=getattr(self.work_params, self.pdbf_def))
    setattr(self, self.needed_info["model"], pdb_inp.construct_hierarchy())
