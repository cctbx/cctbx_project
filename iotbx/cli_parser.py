from __future__ import division, print_function

'''
Standard command-line parser for CCTBX programs

The CCTBXParser class will read files and process PHIL parameters from the
command-line as well as have standard command-line flags for showing the
PHIL scope and citations for a program.

'''

import argparse, getpass, logging, os, sys, time

import iotbx.phil

from iotbx.file_reader import any_file
from libtbx import citations
from libtbx.data_manager import DataManager
from libtbx.utils import Sorry

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

    files = list()
    phil = list()
    directory = list()
    unknown = list()

    if ( parse_files and (getattr(namespace, 'files') is not None) ):
      files = namespace.files
    if ( parse_phil and (getattr(namespace, 'phil') is not None) ):
      phil = namespace.phil
    if ( parse_dir and (getattr(namespace, 'dir') is not None) ):
      directory = namespace.dir

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
    if (parse_phil):
      setattr(namespace, 'phil', phil)
    if (parse_dir):
      setattr(namespace, 'dir', directory)

    # raise Sorry with unknown values
    if (len(unknown) > 0):
      error_message = 'The following arguments are not recognized'
      if (parse_files or parse_dir or parse_phil):
        error_message += ' ('
      if (parse_files):
        error_message += ' file '
      if (parse_dir):
        error_message += ' directory '
      if (parse_phil):
        error_message += ' PHIL '
      if (parse_files or parse_dir or parse_phil):
        error_message += ')'
      error_message += ':\n'
      for value in unknown:
        error_message += '  %s\n' % value
      raise Sorry(error_message)

# =============================================================================
class CCTBXParser(ParserBase):

  def __init__(self, program_class, logger=None, *args, **kwargs):
    '''
    '''
    self.prog = os.getenv('LIBTBX_DISPATCHER_NAME')
    if (self.prog is None):
      self.prog = sys.argv[0]
    self.prefix = self.prog.split('.')[-1]

    # PHIL filenames
    self.data_filename = self.prefix + '_data.eff'
    self.modified_filename = self.prefix + '_modified.eff'
    self.all_filename = self.prefix + '_all.eff'

    border = '-' * 79
    description = border + program_class.description + border
    epilog = border + program_class.epilog
    super(CCTBXParser, self).__init__(
      prog=self.prog, description=description, epilog=epilog,
      formatter_class=argparse.RawDescriptionHelpFormatter,
      *args, **kwargs)

    self.program_class = program_class
    self.logger = logger
    if (self.logger is None):
      self.logger = logging.getLogger('main')
    self.data_manager = DataManager(datatypes=program_class.datatypes)
    self.master_phil = iotbx.phil.parse(program_class.master_phil_str,
                                        process_includes=True)
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
      nargs='?', const=0, type=int, choices=range(0,4),
      help='show default parameters with expert level (default=0)')

    # --attributes-level by itself is set to 1
    # --attributes-level=n sets it to n and it can only be {1, 2, 3}
    self.add_argument(
      '--attributes-level', '--attributes_level',
      nargs='?', const=1, type=int, choices=range(0,4),
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
    print('='*79, file=self.logger)
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

    if (self.working_phil is None):
      self.working_phil = self.master_phil

    return self.namespace

  # ---------------------------------------------------------------------------
  def process_files(self, file_list):
    '''
    Second pass to process files. The first pass already checked that these
    files exist. There may be conditions where the file is deleted in the time
    between the first pass and calling this function.

    Use iotbx.file_reader.any_file to process files.
    Will need updating to work with mmtbx.model.manager class more efficiently
    May be rolled into DataManager class
    '''
    print('Processing files:', file=self.logger)
    print('-'*79, file=self.logger)
    printed_something = False

    unused_files = list()

    for filename in file_list:
      a = any_file(filename)
      # models
      if (a.file_type == 'pdb'):
        self.data_manager.process_model_file(filename)
        print('  Found model, %s' % filename, file=self.logger)
        printed_something = True
      # sequences
      elif (a.file_type == 'seq'):
        self.data_manager.add_sequence(filename, a.file_object)
        print('  Found sequence, %s' % filename, file=self.logger)
        printed_something = True
      elif (a.file_type == 'phil'):
        self.data_manager.add_phil(filename, a.file_object)
        print('  Found PHIL, %s' % filename)
        printed_something = True
      # more file types to come!
      else:
        unused_files.append(filename)

    # show unrecognized files
    if (len(unused_files) > 0):
      print('  Unrecognized files:', file=self.logger)
      print('  -------------------', file=self.logger)
      for filename in unused_files:
        print('  %s' % filename, file=self.logger)
      printed_something = True

    # process PHIL files for DataManager scope in ascending order
    # files are appended and the default is not overridden
    # files from the command-line take precedence
    phil_names = self.data_manager.get_phil_names()
    phil_names.sort()
    for name in phil_names:
      phil = self.data_manager.get_phil(name)
      if (hasattr(phil.extract(), 'data_manager')):
        self.data_manager.load_phil_scope(phil)

    if (not printed_something):
      print('  No files found', file=self.logger)

    if (self.namespace.write_data):
      with open(self.data_filename, 'w') as f:
        self.data_manager.export_phil_scope().show(out=f)

    print('', file=self.logger)

  # ---------------------------------------------------------------------------
  def process_phil(self, phil_list):
    ''''
    Process PHIL arguments
    Currently only handles command-line changes
    Will add inclusion of PHIL files from data_manager first, then
    command-line options
    '''
    print('Processing PHIL parameters:')
    print('-'*79, file=self.logger)
    printed_something = False

    data_sources = list()
    sources = list()
    unused_phil = list()

    # PHIL files are processed in ascending order
    if (self.data_manager.has_phils()):
      phil_names = self.data_manager.get_phil_names()
      phil_names.sort()
      phil = list()
      for name in phil_names:
        phil.append(self.data_manager.get_phil(name))
      data_sources.extend(phil)

    # command-line PHIL arguments override any previous settings
    def custom_processor(arg):
      unused_phil.append(arg)
      return True

    interpreter = self.master_phil.command_line_argument_interpreter(
      home_scope='')
    working = interpreter.process_args(
      phil_list, custom_processor=custom_processor)
    if (len(working) > 0):
      sources.extend(working)
    self.working_phil = self.master_phil.fetch(
      sources=data_sources + sources)

    # show differences
    if (len(sources) > 0):
      phil_diff = self.master_phil.fetch_diff(self.working_phil)
      print('  Non-default PHIL parameters:', file=self.logger)
      print('  ----------------------------', file=self.logger)
      phil_diff.show(prefix='  ', out=self.logger)
      printed_something = True

      # write differences (no DataManager scope)
      if (self.namespace.write_modified):
        with open(self.modified_filename, 'w') as f:
          phil_diff.show(out=f)

    # show unrecognized parameters
    if (len(unused_phil) > 0):
      print('  Unrecognized PHIL parameters:', file=self.logger)
      print('  -----------------------------', file=self.logger)
      for phil in unused_phil:
        print('  %s' % phil, file=self.logger)
      printed_something = True

    if (not printed_something):
      print('  No PHIL parameters found', file=self.logger)

    # write all parameters (DataManager + Program)
    if (self.namespace.write_all):
      with open(self.all_filename, 'w') as f:
        self.data_manager.export_phil_scope().show(out=f)
        self.working_phil.show(expert_level=3, out=f)

    print('', file=self.logger)

  # ---------------------------------------------------------------------------
  def process_dir(self, dir_list):
    '''
    '''
    print('Processing directories:')
    print('-'*79, file=self.logger)
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
      print('Citation(s) for this program:', file=self.logger)
      print('-'*79, file=self.logger)
      print('', file=self.logger)
      for citation in program_citations:
        citations.show_citation(citation, out=self.logger,
                                format=self.namespace.citations)
    self.show_cctbx_citation()

  # ---------------------------------------------------------------------------
  def show_cctbx_citation(self):
    print('\nGeneral citation for CCTBX:', file=self.logger)
    print('-'*79, file=self.logger)
    print('', file=self.logger)
    citations.show_citation(citations.citations_db['cctbx'], out=self.logger,
                            format=self.namespace.citations)

# =============================================================================
# end
