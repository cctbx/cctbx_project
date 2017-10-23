#!/usr/bin/env python
#
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.process

from __future__ import division

help_message = '''

See dials.stills_process.

'''

from libtbx.phil import parse
control_phil_str = '''
  verbosity = 1
    .type = int(value_min=0)
    .help = "The verbosity level"

  input {
    template = None
      .type = str
      .help = "The image sweep template"
      .multiple = True
   }

  output {
    datablock_filename = datablock.json
      .type = str
      .help = "The filename for output datablock"

    strong_filename = strong.pickle
      .type = str
      .help = "The filename for strong reflections from spot finder output."

    indexed_filename = indexed.pickle
      .type = str
      .help = "The filename for indexed reflections."

    refined_experiments_filename = refined_experiments.json
      .type = str
      .help = "The filename for saving refined experimental models"

    integrated_filename = integrated.pickle
      .type = str
      .help = "The filename for final integrated reflections."

    profile_filename = profile.phil
      .type = str
      .help = "The filename for output reflection profile parameters"

    integration_pickle = int-%d-%s.pickle
      .type = str
      .help = Filename for cctbx.xfel-style integration pickle files
  }
'''

delete_shoeboxes_override_str = '''
  integration {
    debug {
      delete_shoeboxes = True
        .type = bool
        .help = "Delete shoeboxes immediately before saving files. This option"
                "in combination with debug.output=True enables intermediate"
                "processing steps to make use of shoeboxes."
    }
  }
'''

from dials.command_line.stills_process import Processor, dials_phil_str, program_defaults_phil_str, Script as DialsScript

phil_scope = parse(control_phil_str + dials_phil_str, process_includes=True).fetch(parse(program_defaults_phil_str))
phil_scope = phil_scope.fetch(parse(delete_shoeboxes_override_str))

class Script(DialsScript, Processor):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    self.reference_detector = None

    # The script usage
    usage = "usage: %s [options] [param.phil] datablock.json" % libtbx.env.dispatcher_name

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message,
      read_datablocks=True,
      read_datablocks_from_images=True)

  def run(self):
    '''Execute the script.'''
    from dxtbx.datablock import DataBlockTemplateImporter
    from dials.util.options import flatten_datablocks
    from dials.util import log
    from logging import info
    from time import time
    from libtbx.utils import Abort

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=False)
    datablocks = flatten_datablocks(params.input.datablock)

    # Check we have some filenames
    if len(datablocks) == 0:

      # Check if a template has been set and print help if not, otherwise try to
      # import the images based on the template input
      if len(params.input.template) == 0:
        self.parser.print_help()
        exit(0)
      else:
        importer = DataBlockTemplateImporter(
          params.input.template,
          options.verbose)
        datablocks = importer.datablocks

    # Save the options
    self.options = options
    self.params = params
    self.load_reference_geometry()

    st = time()

    # Import stuff
    if len(datablocks) == 0:
      raise Abort('No datablocks specified')
    elif len(datablocks) > 1:
      raise Abort('Only 1 datablock can be processed at a time.')
    datablock = datablocks[0]

    if self.reference_detector is not None:
      for imageset in datablock.extract_imagesets():
        imageset.set_detector(self.reference_detector)

    # Configure logging
    log.config(
      params.verbosity,
      info='dials.process.log',
      debug='dials.process.debug.log')

    # Log the diff phil
    diff_phil = self.parser.diff_phil.as_str()
    if diff_phil is not '':
      info('The following parameters have been modified:\n')
      info(diff_phil)

    if self.params.output.datablock_filename:
      from dxtbx.datablock import DataBlockDumper
      dump = DataBlockDumper(datablock)
      dump.as_json(self.params.output.datablock_filename)

    # Do the processing
    observed = self.find_spots(datablock)
    experiments, indexed = self.index(datablock, observed)
    experiments, indexed = self.refine(experiments, indexed)
    integrated = self.integrate(experiments, indexed)

    # Total Time
    info("")
    info("Total Time Taken = %f seconds" % (time() - st))

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
