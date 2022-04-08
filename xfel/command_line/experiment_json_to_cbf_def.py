from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME cxi.experiment_json_to_cbf_def

# Script to convert the output from a joint refinement using dials.refine to a CSPAD
# cbf header file. Note hardcoded distance of 100 isn't relevant for just a cbf header

from dials.util import show_mail_on_error
from dials.util.options import ArgumentParser
from dials.util.options import flatten_experiments
from xfel.cftbx.detector.cspad_cbf_tbx import write_cspad_cbf, map_detector_to_basis_dict
from libtbx import phil

phil_scope = phil.parse("""
  output_def_file = refined_detector.def
    .type = str
    .help = Name of output .def file
""")


class Script(object):
  def __init__(self):
    # Create the parser
    self.parser = ArgumentParser(
      phil = phil_scope,
      read_experiments = True,
      check_format = False)

  def run(self):
    params, options = self.parser.parse_args(show_diff_phil=True)
    experiments = flatten_experiments(params.input.experiments)

    detector = experiments[0].detector

    metro = map_detector_to_basis_dict(detector)
    write_cspad_cbf(None, metro, 'cbf', None, params.output_def_file, None, detector.hierarchy().get_distance(), header_only=True)

    print("Done")

if __name__ == '__main__':
  with show_mail_on_error():
    script = Script()
    script.run()

