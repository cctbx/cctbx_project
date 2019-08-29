#
# This is functionality implementation. It should be easily callable from inside
# other, bigger functions/pipelines.
#
# This file is not supposed to be executed. So no if __name__ == "__main__" here.
#
# It should not care about:
# - CL parsing
# - Basic validation (like file is not there)
# - No filenames available here at all - everything is handled outside.
# - In a rare cases, when absolutely necessary it is still OK to write
#   _INTERMEDIATE_ files to disc here.
# - Results should be in internal objects so they could be used differently
#   by CL-tool, GUI or calling script (pipeline).
#
#
# It should get:
# - parsed and filled out phil parameters
# - objects to work with (model, data, GRM, etc) as internal structures,
#   preferably consistent between each other, passed as arguments
#
# It should do:
# - Specific quck validation: e.g. cannot work with model with alternative
#   conformations to determine beforehand if the job can be done and fail
#   quckly and in uniform way. The result of such validation is up for
#   discussion (raising exception, returning value, etc).
# - The job
#
#
# This kind of files sit in regular places according to functionality:
# cctbx/..., mmtbx/..., etc
#

from __future__ import absolute_import, division, print_function
import iotbx.phil
from libtbx import group_args
from libtbx.utils import Sorry


# Parameters are wrapped with unique top-level name so they could be
# easily imported in higher-level procedures and passed from there
# to this functionality
master_params_str = """
cutter {
  start_resid = None
    .type = int
  end_resid = None
    .type = int
  chain_id = None
    .type = str
}
"""


def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

class cutter():
  """ This class is for cutting out pieces of chain out of hierarchy provided
  chain id, start resid and end resid."""
  def __init__(self, pdb_h, params, log):
    # just basic setup.
    self.pdb_h = pdb_h
    self.params = params
    self.log = log

  def validate_inputs(self):
    # Very specific validation goes here:
    # How validation fails need to be determined.
    # Sorry, return value, assertions allowed?
    if self.params.start_resid is None or self.params.end_resid is None:
      raise Sorry("Need to specify start_resid and start_resid params")
    if len(self.params.chain_id) > 3:
      raise Sorry("Do not support long chain ids yet.")

  def run(self):
    # here goes actual work
    asc = self.pdb_h.atom_selection_cache()
    sel = asc.selection("chain %s and resid %d through %d" % (
        self.params.chain_id, self.params.start_resid, self.params.end_resid))
    self.answer = self.pdb_h.select(sel)

  def get_results(self):
    # small, easy to understand function preparing results. This way other
    # developers will not have to read run() to figure out what is returning.
    # On top of that, this function could be overloaded if necessary!
    return group_args(
        answer = self.answer)

