#
# This is a template for command-line callable script.
# It is a plain script because nothing else should import it anywhere.
# Ideology:
# Inherit application from base application class (cl_app), where
# majority of common functionality is defined.
# Specify needs of particular application here:
# - help message
# - need for model, data, restraints etc in a simple form (see self.needed_info)
# - Nothing should be validated here, because it will not be accessible in
#     other parts, e.g. GUI, pipeline.
# - Define run() method, call actual functionality with available info
# - output results.
#
# Unique things here are:
# - help message
# - what functionality it should run
# - what is needed for the functionality: model, restraints, data (what kind)?
#
# This kind of files will be in */command_line/<name>.py
#
#

from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME phenix.cl_app_example

import sys
from cutter import cutter
from cl_app import cl_app
import iotbx.phil

# Note that we just including existing scope here, which is wrapped in
# single name "cutter" and will be easily passed to the actual class
master_params_str = """\
include scope cutter.master_params

output_prefix = "outp"
  .type = str
"""

# def master_params():
#   return iotbx.phil.parse(master_params_str, process_includes=True)

class cl_cutter(cl_app):
  def __init__(self, cl_args):
    # Set up particular application needs here
    self.help_message = """\
This tool will cut out range of residues from chain. Usage:
  %s model.pdb chain_id='A' start_resid=10, end_resid=20  """ % cl_args[0]
    self.log = sys.stdout

    self.master_params_str = master_params_str
    self.needed_info = {
      "model" : "rmodel", # model-predefined word to indicate the meaning of parameter, rmodel-desired name of class attribute
    }

    # and then just call parent's __init__ to take care of the rest:
    super(cl_cutter, self).__init__(
      cl_args=cl_args)

  def run(self):
    # call parent's function to parse CL. Examine results -> fast failing for
    # typos in parameters etc, IOErrors.
    r = self.parse_cl()
    if r != 0:
      return
    # call parent's function to actually construct objects -> not so fast,
    # could be doing restraints etc.
    r = self.read_and_validate_inputs()
    if r != 0:
      return

    # Now we should have all necessary parsed objects in internal structures
    # like pdb_h etc
    print("Parsed model in form of pdb_hierarchy:", self.rmodel, file=self.log)

    #
    # So do the job like we are in pipeline
    # If needed, this could be wrapped in try...except to catch errors.
    c = cutter(self.rmodel, self.work_params.cutter, self.log) # super-fast init
    # validation if the work can be done.
    # The outcome is to be determined.
    c.validate_inputs()
    c.run() # actually work, may be time-consuming, maybe call-backs are needed
    result = c.get_results()
    # implementation of run() could be long and complicated, so let's do
    # results preparation separately, so everybody could understand what's there
    # at a glance.
    # If something is not there - go to the class and implement new function like:
    # extra_result_value = c.get_extra_result_value()
    # Or use cutting-edge technology:
    # inherit from cutter, overload get_results() and adjust them for your needs.

    #
    # That's it. We got result, now we just need to display it
    # and the job of CL-tool is over.
    print("Result, in form of pdb_hierarchy", result.answer, file=self.log)
    result.answer.write_pdb_file(file_name="%s.pdb" % self.work_params.output_prefix)

cutter_app = cl_cutter(
    cl_args=sys.argv)
cutter_app.run()
