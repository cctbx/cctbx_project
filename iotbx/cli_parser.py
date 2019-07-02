from __future__ import absolute_import, division, print_function
from six.moves import range

'''
Standard command-line parser for CCTBX programs

The CCTBXParser class will read files and process PHIL parameters from the
command-line as well as have standard command-line flags for showing the
PHIL scope and citations for a program.

'''

import argparse, getpass, logging, os, sys, time

import iotbx.phil
import libtbx.phil

from iotbx.data_manager import DataManager, data_manager_type
from iotbx.file_reader import any_file
from libtbx import citations
from libtbx.program_template import ProgramTemplate
from libtbx.str_utils import wordwrap
from libtbx.utils import multi_out, show_times, Sorry

# =============================================================================
def run_program(program_class=None, custom_process_arguments=None,
                args=None, logger=None):
  '''
  Function for running programs using CCTBXParser and the program template

  :param program_class:  ProgramTemplate type (required)
  :param custom_process_arguments:
                         Custom function to parse unknown arguments (optional)
  :param args:           list of command-line arguments (optional)
  :param logger:         logger (e.g. multi_out) for output (optional)
  :rtype:                whatever is returned from program_class.get_results()
  '''

  assert (program_class is not None)

  if (args is None):
    args = sys.argv[1:]

  # create logger
  if (logger is None):
    logger = multi_out()
    logger.register('stdout', sys.stdout)

  # start timer
  t = show_times(out=logger)

  # create parser
  parser = CCTBXParser(program_class=program_class,
                       custom_process_arguments=custom_process_arguments,
                       logger=logger)
  namespace = parser.parse_args(args)

  # start program
  print('Starting job', file=logger)
  print('='*79, file=logger)
  task = program_class(parser.data_manager, parser.working_phil.extract(),
                       master_phil=parser.master_phil,
                       logger=logger)

  # custom constructor (optional)
  task.custom_init()

  # validate inputs
  task.validate()

  # run program
  task.run()

  # clean up (optional)
  task.clean_up()

  # stop timer
  print('', file=logger)
  print('='*79, file=logger)
  print('Job complete', file=logger)
  t()

  return task.get_results()

# =============================================================================
class ParserBase(argparse.ArgumentParser):

  def __init__(self, parse_files=True, parse_phil=True, parse_dir=False,
               *args, **kwargs):
    super(ParserBase, self).__init__(*args, **kwargs)

    # store options
    self.parse_files = parse_files
    self.parse_phil = parse_phil
    self.parse_dir = parse_dir

    # add default behavior for positional arguments
    if (self.parse_files):
      self.add_argument('files', nargs='*', help='Input file(s) (e.g. model.cif)',
                        action=ParsePositionalArgumentsAction)
    if (self.parse_phil):
      self.add_argument('phil', nargs='*', help='Parameter(s) (e.g. d_min=2.0)',
                        action=ParsePositionalArgumentsAction)
    if (self.parse_dir):
      self.add_argument('dir', nargs='*', help='Input directory',
                        action=ParsePositionalArgumentsAction)

# =============================================================================
class ParsePositionalArgumentsAction(argparse.Action):
  '''
  This action is a first pass for command-line arguments. It does basic checks
  to see if an argument is a file, a directory, or a phil parameter (contains
  an equals sign). Command-line switches (options beginning with "-") are
  handled by default actions in the parser
  '''

  def __call__(self, parser, namespace, values, option_string=None):
    # figure out what options are set
    parse_files = hasattr(namespace, 'files')
    parse_phil = hasattr(namespace, 'phil')
    parse_dir = hasattr(namespace, 'dir')

    # get previous values or define default
    if ( parse_files and (getattr(namespace, 'files') is not None) ):
      files = namespace.files
    else:
      files = list()
    if ( parse_phil and (getattr(namespace, 'phil') is not None) ):
      phil = namespace.phil
    else:
      phil = list()
    if ( parse_dir and (getattr(namespace, 'dir') is not None) ):
      directory = namespace.dir
    else:
      directory = list()
    if ( hasattr(namespace, 'unknown') and
         (getattr(namespace, 'unknown') is not None) ):
      unknown = namespace.unknown
    else:
      unknown = list()

    # separate values
    for value in values:
      if (os.path.isfile(value)):
        files.append(value)
      elif ( (isinstance(value, str)) and
             ('=' in value) ):
        phil.append(value)
      elif (os.path.isdir(value)):
        directory.append(value)
      else:
        unknown.append(value)

    # update options
    if (parse_files):
      setattr(namespace, 'files', files)
    else:
      unknown.extend(files)
    if (parse_phil):
      setattr(namespace, 'phil', phil)
    else:
      unknown.extend(phil)
    if (parse_dir):
      setattr(namespace, 'dir', directory)
    else:
      unknown.extend(directory)

    # store unknown values for custom processing (if available)
    setattr(namespace, 'unknown', unknown)

# =============================================================================
class CCTBXParser(ParserBase):

  def __init__(self, program_class, custom_process_arguments=None,
               logger=None, *args, **kwargs):
    '''
    '''
    # program name
    self.prog = os.getenv('LIBTBX_DISPATCHER_NAME')
    if (self.prog is None):
      self.prog = sys.argv[0]
    self.prefix = self.prog.split('.')[-1]

    # PHIL filenames
    self.data_filename = self.prefix + '_data.eff'
    self.modified_filename = self.prefix + '_modified.eff'
    self.all_filename = self.prefix + '_all.eff'

    # terminal width
    self.text_width = 79

    # print header
    border = '-' * self.text_width
    description = border + program_class.description + border
    epilog = border + program_class.epilog
    super(CCTBXParser, self).__init__(
      prog=self.prog, description=description, epilog=epilog,
      formatter_class=argparse.RawDescriptionHelpFormatter,
      *args, **kwargs)

    # default values
    self.program_class = program_class
    self.custom_process_arguments = custom_process_arguments
    self.logger = logger
    if (self.logger is None):
      self.logger = logging.getLogger('main')
    self.data_manager = DataManager(datatypes=program_class.datatypes,
                                    logger=self.logger)

    # add PHIL converters if available
    if (len(program_class.phil_converters) > 0):
      iotbx.phil.default_converter_registry = \
        libtbx.phil.extended_converter_registry(
          additional_converters=program_class.phil_converters,
          base_registry=iotbx.phil.default_converter_registry)

    # set up master and working PHIL scopes
    self.master_phil = iotbx.phil.parse(
      program_class.master_phil_str, process_includes=True)
    required_output_phil = iotbx.phil.parse(ProgramTemplate.output_phil_str)
    self.master_phil.adopt_scope(required_output_phil)
    self.working_phil = None

    self.add_default_options()

  # ---------------------------------------------------------------------------
  def add_default_options(self):
    '''
    '''
    # --show-defaults by itself is set to 0
    # --show-defaults=n sets it to n and it can only be {0, 1, 2, 3}
    self.add_argument(
      '--show-defaults', '--show_defaults',
      nargs='?', const=0, type=int, choices=list(range(0,4)),
      help='show default parameters with expert level (default=0)')

    # --attributes-level by itself is set to 1
    # --attributes-level=n sets it to n and it can only be {1, 2, 3}
    self.add_argument(
      '--attributes-level', '--attributes_level',
      nargs='?', const=1, type=int, choices=list(range(0,4)),
      help='show parameters with attributes (default=0)'
    )

    # --write-data
    # switch for writing only DataManager PHIL parameters
    self.add_argument(
      '--write-data', '--write_data', action='store_true',
      help='write DataManager PHIL parameters to file (%s)' % \
      self.data_filename
    )

    # --write-modified
    # switch for writing only modified PHIL parameters
    self.add_argument(
      '--write-modified', '--write_modified', action='store_true',
      help='write modifed PHIL parameters to file (%s)' % \
      self.modified_filename
    )

    # --write-all
    # switch for writing all PHIL parameters
    self.add_argument(
      '--write-all', '--write_all', action='store_true',
      help='write all (modified + default + data) PHIL parameters to file (%s)' %
      self.all_filename
    )

    # --overwrite
    # switch for overwriting files, takes precedence over PHIL definition
    self.add_argument(
      '--overwrite', action='store_true', default=False,
      help='overwrite files, this overrides the output.overwrite PHIL parameter'
    )

    # --citations will use the default format
    # --citations=<format> will use the specified format
    self.add_argument(
      '--citations',
      nargs='?', const='default', type=str, choices=['default', 'cell', 'iucr'],
      help='show citation(s) for program in different formats')

  # ---------------------------------------------------------------------------
  def parse_args(self, args):
    '''
    '''
    # default behavior with no arguments
    if (len(args) == 0):
      self.print_help()
      self.exit()

    # parse arguments
    self.namespace = super(CCTBXParser, self).parse_args(args)

    # process command-line options
    if (self.namespace.attributes_level is not None):
      if (self.namespace.show_defaults is None):
        self.error('--attributes-level requires --show-defaults to be set')
    if (self.namespace.show_defaults is not None):
      self.master_phil.show(expert_level=self.namespace.show_defaults,
                            attributes_level=self.namespace.attributes_level,
                            out=self.logger)
      self.exit()
    if (self.namespace.citations is not None):
      self.show_citations()
      self.exit()

    # start program header
    print ('Starting %s' % self.prog, file=self.logger)
    print('on %s by %s' % (time.asctime(), getpass.getuser()), file=self.logger)
    print('='*self.text_width, file=self.logger)
    print('', file=self.logger)

    # process files
    if (self.parse_files):
      self.process_files(self.namespace.files)

    # process phil and phil files
    if (self.parse_phil):
      self.process_phil(self.namespace.phil)

    # process directories
    if (self.parse_dir):
      self.process_dir(self.namespace.dir)

    # custom processing of arguments (if available)
    # the function for custom processing of arguments should take a CCTBXParser
    # object as its argument. The function should modify the
    # CCTBXParser.namespace.unknown, CCTBXParser.data_manager,
    # CCTBXParser.working_phil, or other CCTBXParser members.
    # At the end of the function, CCTBXParser.working_phil should have the final
    # libtbx.phil.scope object (not libtbx.phil.scope_extract) for the program
    # A libtbx.phil.scope_extract object can be converted into a
    # libtbx.phil.scope object by
    #    CCTBXParser.master_phil.format(python_object=<scope_extract>)
    if (self.custom_process_arguments is not None):
      self.custom_process_arguments(self)
      assert(isinstance(self.working_phil, libtbx.phil.scope))

    # post processing after all arguments are parsed
    self.post_process()

    # show final PHIL parameters
    self.show_phil_summary()

    return self.namespace

  # ---------------------------------------------------------------------------
  def process_files(self, file_list, message = 'Processing files:'):
    '''
    Second pass to process files. The first pass already checked that these
    files exist. There may be conditions where the file is deleted in the time
    between the first pass and calling this function.

    Use iotbx.file_reader.any_file to process files.
    Will need updating to work with mmtbx.model.manager class more efficiently
    '''
    print(message, file=self.logger)
    print('-'*self.text_width, file=self.logger)
    print('', file=self.logger)
    printed_something = False

    unused_files = list()

    for filename in file_list:
      a = any_file(filename)
      process_function = 'process_%s_file' % data_manager_type.get(a.file_type)
      if (hasattr(self.data_manager, process_function)):
        getattr(self.data_manager, process_function)(filename)
        print('  Found %s, %s' % (data_manager_type[a.file_type], filename),
              file=self.logger)
        printed_something = True
      else:
        unused_files.append(filename)

    # show unrecognized files
    if (len(unused_files) > 0):
      if (printed_something):
        print('', file=self.logger)
      print('  Files not used by program:', file=self.logger)
      print('  --------------------------', file=self.logger)
      for filename in unused_files:
        print('  %s' % filename, file=self.logger)
      printed_something = True

    # process PHIL files for DataManager scope in order from command-line
    # files are appended and the default is not overridden
    # files from the command-line take precedence
    phil_names = self.data_manager.get_phil_names()
    for name in phil_names:
      phil = self.data_manager.get_phil(name)
      if (hasattr(phil.extract(), 'data_manager')):
        self.data_manager.load_phil_scope(phil)

    if (not printed_something):
      print('  No files found', file=self.logger)

    print('', file=self.logger)

  # ---------------------------------------------------------------------------
  def process_phil(self, phil_list):
    ''''
    Process PHIL arguments
    Also checks PHIL arguments (command line and files) for parameters that
    specify files (.type = path)
    '''
    print('Processing PHIL parameters:', file=self.logger)
    print('-'*self.text_width, file=self.logger)
    print('', file=self.logger)

    printed_something = False

    data_sources = list()
    sources = list()
    unused_phil = list()

    # PHIL files are processed in order from command-line
    if (self.data_manager.has_phils()):
      phil_names = self.data_manager.get_phil_names()
      phil = list()
      print('  Adding PHIL files:', file=self.logger)
      print('  ------------------', file=self.logger)
      for name in phil_names:
        # remove DataManager scope since input files are already loaded
        phil_scope = self.data_manager.get_phil(name)
        for phil_object in phil_scope.objects:
          if (phil_object.name == 'data_manager'):
            phil_scope.objects.remove(phil_object)
        phil.append(phil_scope)
        print('    %s' % name, file=self.logger)
      data_sources.extend(phil)
      print('', file=self.logger)
      printed_something = True

    # command-line PHIL arguments override any previous settings and are
    # processed in given order
    def custom_processor(arg):
      unused_phil.append(arg)
      return True

    if (len(phil_list) > 0):
      interpreter = self.master_phil.command_line_argument_interpreter(
        home_scope='')
      print('  Adding command-line PHIL:', file=self.logger)
      print('  -------------------------', file=self.logger)
      for phil in phil_list:
        print('    %s' % phil, file=self.logger)
      print('', file=self.logger)
      printed_something = True
      working = interpreter.process_args(
        phil_list, custom_processor=custom_processor)
      if (len(working) > 0):
        sources.extend(working)
    if (self.namespace.overwrite):  # override overwrite if True
      sources.append(iotbx.phil.parse('output.overwrite=True'))
    if ( (len(data_sources) + len(sources)) > 0):
      self.working_phil, more_unused_phil = self.master_phil.fetch(
        sources=data_sources + sources, track_unused_definitions=True)
      unused_phil.extend(more_unused_phil)
    else:
      self.working_phil = self.master_phil.fetch()

    # show unrecognized parameters and abort
    if (len(unused_phil) > 0):
      print('  Unrecognized PHIL parameters:', file=self.logger)
      print('  -----------------------------', file=self.logger)
      for phil in unused_phil:
        print('    %s' % phil, file=self.logger)
      print('', file=self.logger)
      error_message = 'Some PHIL parameters are not recognized by %s.\n' % \
                      self.prog
      error_message += wordwrap('Please run this program with the --show-defaults option to see what parameters are available.', max_chars=self.text_width) + '\n'
      error_message += wordwrap('PHIL parameters in files should be fully specified (e.g. "output.overwrite" instead of just "overwrite")', max_chars=self.text_width) + '\n'
      raise Sorry(error_message)

    # process input phil for file/directory defintions and add to DataManager
    # Note: if a PHIL file is input as a PHIL parameter, the contents of the
    # file will NOT be parsed and validated. The PHIL file should be provided
    # as a command-line argument. This is mostly for finding data files
    # defined by PHIL parameters that should be added to the DataManager
    diff_phil = self.master_phil.fetch_diff(self.working_phil)
    paths = self.check_phil_for_paths(diff_phil)
    if (len(paths) > 0):
      files = list()
      dirs = list()
      for path in paths:
        if (path is not None):
          if (os.path.isfile(path)):
            files.append(path)
          elif (os.path.isdir(path)):
            dirs.append(path)
      if (self.parse_files):
        self.process_files(files, message='Processing files from PHIL:')
      if (self.parse_dir):
        self.process_dir(dirs, message='Processing directories from PHIL:')

    if (not printed_something):
      print('  No PHIL parameters found', file=self.logger)
      print('', file=self.logger)

  # ---------------------------------------------------------------------------
  def process_dir(self, dir_list, message = 'Processing directories:'):
    '''
    '''
    print(message, file=self.logger)
    print('-'*self.text_width, file=self.logger)
    print('', file=self.logger)

  # ---------------------------------------------------------------------------
  def check_phil_for_paths(self, phil_scope):
    '''
    Recursively check PHIL scope if there is a 'path' type.
    Returns the paths (empty list means no paths were found)
    '''
    paths = list()
    if (phil_scope.is_definition):
      if (phil_scope.type.phil_type == 'path'):
        if phil_scope.style is not None and 'new_file' in phil_scope.style:
          pass
        else:
          paths.append(phil_scope.extract())
    elif (phil_scope.is_scope):
      for phil_object in phil_scope.objects:
        paths.extend(self.check_phil_for_paths(phil_object))
    return paths

  # ---------------------------------------------------------------------------
  def post_process(self):
    '''
    Post processing of inputs after all arguments are parsed
    '''

    working_phil_extract = self.working_phil.extract()

    # update default model with program pdb interpretation scope
    if (self.data_manager.supports('model') and
        (self.data_manager.get_default_model_name() is not None)):
      self.data_manager.update_pdb_interpretation_for_model(
        self.data_manager.get_default_model_name(), working_phil_extract)

  # ---------------------------------------------------------------------------
  def show_phil_summary(self):
    '''
    Show final, modified PHIL parameters after all processing is complete
    Also, write phil scopes based on command-line flags
    '''

    overwrite = (self.namespace.overwrite or \
                 self.working_phil.extract().output.overwrite)

    # check for any remaining unknown arguments
    if (len(self.namespace.unknown) > 0):
      error_message = 'The following arguments are not recognized:\n'
      for value in self.namespace.unknown:
        error_message += '  %s\n' % value
      raise Sorry(error_message)

    # get differences
    try:
      data_diff = self.data_manager.master_phil.fetch_diff(
        self.data_manager.export_phil_scope())
    except RuntimeError as err:
      raise Sorry(err)
    try:
      phil_diff = self.master_phil.fetch_diff(self.working_phil)
    except RuntimeError as err:
      raise Sorry(err)
    data_is_different = (len(data_diff.as_str()) > 0)
    phil_is_different = (len(phil_diff.as_str()) > 0)
    is_different = data_is_different or phil_is_different

    # show final processed phil scope
    print('Final processed PHIL parameters:', file=self.logger)
    print('-'*self.text_width, file=self.logger)
    if (is_different):
      data_diff.show(prefix='  ', out=self.logger)
      phil_diff.show(prefix='  ', out=self.logger)
    else:
      print('  All parameters are set to their defaults', file=self.logger)
    print('', file=self.logger)

    # write scopes if requested
    if (self.namespace.write_data or self.namespace.write_modified or
        self.namespace.write_all):
      print('Writing program PHIL file(s):', file=self.logger)

    # write DataManager scope
    if (self.namespace.write_data):
      if (data_is_different):
        self.data_manager.write_phil_file(
          self.data_manager.export_phil_scope().as_str(),
          filename=self.data_filename,
          overwrite=overwrite)
        print('  Input file PHIL written to %s.' % self.data_filename,
              file=self.logger)
      else:
        print('  No input file PHIL to write', file=self.logger)

    # write differences
    if (self.namespace.write_modified):
      if (phil_is_different):
        self.data_manager.write_phil_file(
          phil_diff.as_str(), filename=self.modified_filename,
          overwrite=overwrite)
        print('  Modified PHIL parameters written to %s.' %
              self.modified_filename, file=self.logger)
      else:
        print('  No PHIL modifications to write', file=self.logger)

    # write all parameters (DataManager + Program)
    if (self.namespace.write_all):
      all_phil = self.data_manager.export_phil_scope().as_str()
      all_phil += self.working_phil.as_str(expert_level=3)
      self.data_manager.write_phil_file(
        all_phil, filename=self.all_filename, overwrite=overwrite)
      print('  All PHIL parameters written to %s.' % self.all_filename,
            file=self.logger)

    print('', file=self.logger)

  # ---------------------------------------------------------------------------
  def show_citations(self):
    # build list of program-specific citations
    program_citations = list()
    if (self.program_class.citations is not None):
      class_citations = citations.master_citation_phil.fetch(
        source=self.program_class.citations).extract()
      for citation in class_citations.citation:
        program_citations.append(citation)
    for article_id in self.program_class.known_article_ids:
      citation = citations.citations_db.get(article_id)
      if (citation is not None):
        program_citations.append(citation)
      else:
        raise Sorry('"%s" not found citations database' % article_id)

    # show program-specific citations and general citation for CCTBX
    if (len(program_citations) > 0):
      print('Citation(s) for %s:' % self.prog, file=self.logger)
      print('-'*self.text_width, file=self.logger)
      print('', file=self.logger)
      for citation in program_citations:
        citations.show_citation(citation, out=self.logger,
                                format=self.namespace.citations)
    self.program_class.show_template_citation(
      text_width=self.text_width, logger=self.logger,
      citation_format=self.namespace.citations)

# =============================================================================
# end
