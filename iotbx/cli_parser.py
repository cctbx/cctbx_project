'''
Standard command-line parser for CCTBX programs

The CCTBXParser class will read files and process PHIL parameters from the
command-line as well as have standard command-line flags for showing the
PHIL scope and citations for a program.

'''
from __future__ import absolute_import, division, print_function

import argparse, getpass, logging, os, sys, textwrap, time

from six.moves import cStringIO as StringIO

import iotbx.phil
import libtbx.phil

from iotbx.data_manager import DataManager, data_manager_type
from iotbx.file_reader import any_file
from libtbx import citations
from libtbx.program_template import ProgramTemplate
from libtbx.str_utils import wordwrap
from libtbx.utils import multi_out, null_out, show_times, Sorry

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
    if self.parse_files:
      self.add_argument('files', nargs='*', help='Input file(s) (e.g. model.cif)',
                        action=ParsePositionalArgumentsAction)
    if self.parse_phil:
      self.add_argument('phil', nargs='*', help='Parameter(s) (e.g. d_min=2.0)',
                        action=ParsePositionalArgumentsAction)
    if self.parse_dir:
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
    if parse_files and getattr(namespace, 'files') is not None:
      files = namespace.files
    else:
      files = []
    if parse_phil and getattr(namespace, 'phil') is not None:
      phil = namespace.phil
    else:
      phil = []
    if parse_dir and getattr(namespace, 'dir') is not None:
      directory = namespace.dir
    else:
      directory = []
    if hasattr(namespace, 'unknown') \
      and getattr(namespace, 'unknown') is not None:
      unknown = namespace.unknown
    else:
      unknown = []

    # separate values
    for value in values:
      if os.path.isfile(value):
        files.append(value)
      elif isinstance(value, str) and '=' in value:
        phil.append(value)
      elif os.path.isdir(value):
        directory.append(value)
      else:
        unknown.append(value)

    # update options
    if parse_files:
      setattr(namespace, 'files', files)
    else:
      unknown.extend(files)
    if parse_phil:
      setattr(namespace, 'phil', phil)
    else:
      unknown.extend(phil)
    if parse_dir:
      setattr(namespace, 'dir', directory)
    else:
      unknown.extend(directory)

    # store unknown values for custom processing (if available)
    setattr(namespace, 'unknown', unknown)

# =============================================================================
class CCTBXParser(ParserBase):

  def __init__(self, program_class, custom_process_arguments=None,
               unused_phil_raises_sorry=True, logger=None, *args, **kwargs):
    '''
    '''
    # program name
    # Order of precedence:
    # 1) ProgramTemplate.program_name
    # 2) LIBTBX_DISPATCHER_NAME
    # 3) Calling command
    if hasattr(sys, 'argv') and sys.argv:
      self.prog = os.getenv('LIBTBX_DISPATCHER_NAME', sys.argv[0])
    else:
      self.prog = 'unknown.unknown'
    if program_class.program_name is not None:
      self.prog = program_class.program_name
    if program_class.program_name is None:
      program_class.program_name = self.prog
    self.prefix = self.prog.split('.')[-1]
    # Windows dispatchers may have the .bat extension
    if sys.platform == 'win32' and self.prefix.lower() == 'bat':
      self.prefix = self.prog.split('.')[-2]

    # PHIL filenames
    self.data_filename = self.prefix + '_data.eff'
    self.defaults_filename = self.prefix + '_defaults.eff'
    self.modified_filename = self.prefix + '_modified.eff'
    self.all_filename = self.prefix + '_all.eff'

    # JSON filename
    self.json_filename = self.prefix + '_result.json'

    # terminal width
    self.text_width = 79

    # DataManager diff
    self.data_manager_diff = None

    # print header
    border = '-' * self.text_width
    description = border + textwrap.dedent(program_class.description) + border
    epilog = border + program_class.epilog
    super(CCTBXParser, self).__init__(
      prog=self.prog, description=description, epilog=epilog,
      formatter_class=argparse.RawDescriptionHelpFormatter,
      *args, **kwargs)

    # default values
    self.program_class = program_class
    self.custom_process_arguments = custom_process_arguments
    self.unused_phil_raises_sorry = unused_phil_raises_sorry
    self.unused_phil = []
    self.logger = logger
    if self.logger is None:
      self.logger = logging.getLogger('main')
    self.data_manager = DataManager(
      datatypes=program_class.datatypes,
      custom_options=program_class.data_manager_options,
      custom_master_phil_str=program_class.data_manager_custom_master_phil_str,
      logger=self.logger)
    self.namespace = None

    # add PHIL converters if available
    if len(program_class.phil_converters) > 0:
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

    # add master PHIL to description
    extra_description = '\n\nPHIL arguments:\n'
    phil_str = ''
    if self.program_class.show_data_manager_scope_by_default:
      phil_str += self.data_manager.master_phil.as_str(
        prefix='  ', expert_level=3, attributes_level=0)
    phil_str += self.master_phil.as_str(
      prefix='  ', expert_level=3, attributes_level=0)
    if len(phil_str) > 0:
      extra_description += phil_str
      self.description += extra_description

    self.add_default_options()

  # ---------------------------------------------------------------------------
  def add_default_options(self):
    '''
    '''
    # --show-defaults by itself is set to 0
    # --show-defaults=n sets it to n and it can only be {0, 1, 2, 3}
    self.add_argument(
      '--show-defaults', '--show_defaults',
      nargs='?', const=3, type=int, choices=range(0,4),
      help='show default parameters with expert level (default=3)'
    )

    # --attributes-level by itself is set to 0
    # --attributes-level=n sets it to n and it can only be {0, 1, 2, 3}
    self.add_argument(
      '--attributes-level', '--attributes_level',
      nargs='?', const=0, type=int, choices=list(range(0,4)),
      help='show parameters with extra attributes (default=0)'
    )

    # --write-data
    # switch for writing only DataManager PHIL parameters
    self.add_argument(
      '--write-data', '--write_data', action='store_true',
      help='write DataManager PHIL parameters to file (%s)' %
      self.data_filename
    )

    # --write-defaults
    # switch for writing all the default PHIL parameters
    self.add_argument(
      '--write-defaults', '--write_defaults', action='store_true',
      help='write default PHIL parameters to file (%s)' %
      self.defaults_filename
    )

    # --write-modified
    # switch for writing only modified PHIL parameters
    self.add_argument(
      '--write-modified', '--write_modified', action='store_true',
      help='write modifed PHIL parameters to file (%s)' %
      self.modified_filename
    )

    # --write-all
    # switch for writing all PHIL parameters
    self.add_argument(
      '--write-all', '--write_all', action='store_true',
      help='write all (modified + default + data) PHIL parameters to file (%s)' %
      self.all_filename
    )

    # --diff-params
    # switch for writing the differences between the input PHIL files and
    # the current defaults
    self.add_argument(
      '--diff-params', '--diff_params', action='store_true',
      help='similar to --write-modified, but stops program execution after writing and always overwrites'
    )

    # --json
    # return JSON output from program
    self.add_argument(
      '--json', action='store_true',
      help='''\
writes or overwrites the JSON output for the program to file (%s).
Use --json-filename to specify a different filename for the output.''' %
      self.json_filename,
    )

    # --json-filename
    # set a non-default filename for JSON output
    self.add_argument(
      '--json-filename', '--json_filename', action='store',
      type=str, default=None,
      help='''\
optionally specify a filename for JSON output. If a filename is provided,
the .json extension will be added automatically if it does not already exist.
Also, specifying this flag implies that --json is also specified.'''
    )

    # --check-current-dir
    # if a file does exist in the path specified in the PHIL, check the current directory
    self.add_argument(
      '--check-current-dir', '--check_current_dir', action='store_true',
      help='check current directory for a file if the file path specified in a PHIL file does not work'
    )

    # --overwrite
    # switch for overwriting files, takes precedence over PHIL definition
    self.add_argument(
      '--overwrite', action='store_true',
      help='overwrite files, this overrides the output.overwrite PHIL parameter'
    )

    # --profile
    # enable profiling output
    # the output file is hardcoded to "profile.out" to avoid parsing confusion
    # e.g. --profile model.pdb should pass model.pdb to the program, not dump
    # the profiling stats to model.pdb.
    self.add_argument(
      '--profile', action='store_true',
      help='enable profiling and outputs statistics to a file (profile.out).'
    )

    # --dry-run
    # proceeds until the validate step
    self.add_argument(
      '--dry-run', '--dry_run', action='store_true',
      help='performs basic validation the input arguments, but does not run the program'
    )

    # --citations will use the default format
    # --citations=<format> will use the specified format
    self.add_argument(
      '--citations',
      nargs='?', const='default', type=str, choices=['default', 'cell', 'iucr'],
      help='show citation(s) for program in different formats'
    )

    # --quiet
    # suppress output
    self.add_argument(
      '--quiet', action='store_true',
      help='suppress output to terminal'
    )

    # --version
    # returns the program version
    self.add_argument(
      '--version', '-v', action='store_true',
      help='show version information'
    )

  # ---------------------------------------------------------------------------
  def parse_args(self, args, skip_help = False):
    '''
    '''
    # default behavior with no arguments
    if (len(args) == 0) and (not skip_help):
      self.print_help()
      self.exit()

    # parse arguments
    if sys.version_info >= (3, 7):
      # https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser.parse_intermixed_args
      # https://bugs.python.org/issue9338
      # https://bugs.python.org/issue15112
      self.namespace = super(CCTBXParser, self).parse_intermixed_args(args)
    else:
      self.namespace = super(CCTBXParser, self).parse_args(args)

    # process command-line options
    if self.namespace.attributes_level is not None:
      if self.namespace.show_defaults is None:
        self.error('--attributes-level requires --show-defaults to be set')
    if self.namespace.show_defaults is not None:
      if self.namespace.attributes_level is None:
        self.namespace.attributes_level = 0
      if self.program_class.show_data_manager_scope_by_default:
        self.data_manager.master_phil.show(
          expert_level=self.namespace.show_defaults,
          attributes_level=self.namespace.attributes_level,
          out=self.logger)
      self.master_phil.show(expert_level=self.namespace.show_defaults,
                            attributes_level=self.namespace.attributes_level,
                            out=self.logger)
      self.exit()
    if self.namespace.citations is not None:
      self.show_citations()
      self.exit()
    if self.namespace.version:
      print(self.program_class.get_version(), file=self.logger)
      self.exit()

    # start program header
    print ('Starting %s' % self.prog, file=self.logger)
    print('on %s by %s' % (time.asctime(), getpass.getuser()), file=self.logger)
    print('='*self.text_width, file=self.logger)
    print('', file=self.logger)
    self.logger.flush()

    # process files
    if self.parse_files:
      self.process_files(self.namespace.files)

    # process phil and phil files
    if self.parse_phil:
      self.process_phil(self.namespace.phil)

    # process directories
    if self.parse_dir:
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
    if self.custom_process_arguments is not None:
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

    unused_files = []

    for filename in file_list:
      a = any_file(filename)
      process_function = 'process_%s_file' % data_manager_type.get(a.file_type)
      if hasattr(self.data_manager, process_function):
        getattr(self.data_manager, process_function)(filename)
        print('  Found %s, %s' % (data_manager_type[a.file_type], filename),
              file=self.logger)
        printed_something = True
      else:
        unused_files.append(filename)

    # show unrecognized files
    if len(unused_files) > 0:
      if printed_something:
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
      if hasattr(phil.extract(), 'data_manager'):
        phil = self._update_phil_paths_to_cwd(phil)
        self.data_manager.load_phil_scope(phil, process_files=not self.namespace.diff_params)

    if not printed_something:
      print('  No files found', file=self.logger)

    print('', file=self.logger)

  # ---------------------------------------------------------------------------
  def raise_Sorry_for_unused_phil(self):
    '''
    Convenience function for aborting when there are unused PHIL
    parameters. This function is useful when custom PHIL handling is
    necessary.
    '''
    if len(self.unused_phil) > 0 and self.unused_phil_raises_sorry:
      advice = ''
      print('  Unrecognized PHIL parameters:', file=self.logger)
      print('  -----------------------------', file=self.logger)
      for phil in self.unused_phil:
        print('    %s' % phil, file=self.logger)
        if str(phil).find('.qi.')>-1:
          advice = 'Consider setting a QM package using PHENIX_MOPAC, PHENIX_ORCA or similar.'
      print('', file=self.logger)
      error_message = 'Some PHIL parameters are not recognized by %s.\n' % \
                      self.prog
      error_message += wordwrap('Please run this program with the --show-defaults option to see what parameters are available.', max_chars=self.text_width) + '\n'
      error_message += wordwrap('PHIL parameters in files should be fully specified (e.g. "output.overwrite" instead of just "overwrite")', max_chars=self.text_width) + '\n'
      if advice:
        error_message += wordwrap(advice, max_chars=self.text_width) + '\n'
      if self.unused_phil_raises_sorry and not self.namespace.diff_params:
        raise Sorry(error_message)

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

    # DataManager PHIL
    data_manager_data_sources = []  # from files
    data_manager_sources = []       # from command line

    # program PHIL
    data_sources = []
    sources = []

    # PHIL files are processed in order from command-line
    if self.data_manager.has_phils():
      phil_names = self.data_manager.get_phil_names()
      print('  Adding PHIL files:', file=self.logger)
      print('  ------------------', file=self.logger)
      for name in phil_names:
        # remove DataManager scope since input files are already loaded
        phil_scope = self.data_manager.get_phil(name)
        for phil_object in phil_scope.objects:
          if phil_object.name == 'data_manager':
            phil_scope.objects.remove(phil_object)
            data_manager_data_sources.append(phil_scope.customized_copy(objects=[phil_object]))
        data_sources.append(phil_scope)
        print('    %s' % name, file=self.logger)
      print('', file=self.logger)
      printed_something = True

    # command-line PHIL arguments override any previous settings and are
    # processed in given order
    if len(phil_list) > 0:
      interpreter = self.master_phil.command_line_argument_interpreter(
        assume_when_ambiguous=False)
      data_manager_interpreter = self.data_manager.master_phil.command_line_argument_interpreter(
        assume_when_ambiguous=False)
      print('  Adding command-line PHIL:', file=self.logger)
      print('  -------------------------', file=self.logger)
      for phil in phil_list:
        print('    %s' % phil, file=self.logger)
      print('', file=self.logger)
      printed_something = True
      # check each parameter
      for phil in phil_list:
        processed_arg = None
        data_manager_processed_arg = None
        try:
          processed_arg = interpreter.process_arg(arg=phil)
        except Sorry as e:
          if e.__str__().startswith('Unknown'):
            # check if it is a DataManager parameter
            try:
              data_manager_processed_arg = data_manager_interpreter.process_arg(arg=phil)
            except Sorry as e2:
              if e2.__str__().startswith('Unknown'):
                self.unused_phil.append(phil)
              else:
                raise
          else:
            raise
        if processed_arg is not None:
          sources.append(processed_arg)
        if data_manager_processed_arg is not None:
          data_manager_sources.append(data_manager_processed_arg)

      # collect multiple DataManager label specifications from command-line
      #   data_manager.miller_array.labels.name
      #   data_manager.miller_array.user_selected_labels
      def _get_last_object(phil_scope):
        if phil_scope.is_definition:
          return phil_scope
        return _get_last_object(phil_scope.objects[0])

      label_objects = []
      type_objects = []
      for source in data_manager_sources:
        phil_object = _get_last_object(source)
        if phil_object.full_path() in [
          'data_manager.miller_array.labels.name',
          'data_manager.miller_array.user_selected_labels',
          'data_manager.map_coefficients.labels.name',
          'data_manager.map_coefficients.user_selected_labels']:
          label_objects.append(source)
        if phil_object.full_path() == 'data_manager.model.type':
          type_objects.append(source)

      if len(label_objects) > 0:
        print('  Found labels in command-line', file=self.logger)
        print('  ----------------------------', file=self.logger)

        # combine labels of the same datatype
        datatypes = ['miller_array', 'map_coefficients']
        for datatype in datatypes:
          for parent_label_object in label_objects:
            parent_object = _get_last_object(parent_label_object)  # name / user_selected_labels
            pps = parent_object.primary_parent_scope  # miller_array
            if parent_object.primary_parent_scope.name == 'datatype':
              break
          if parent_object.full_path() == 'data_manager.%s.labels.name' % datatype:
            pps = parent_object.primary_parent_scope.primary_parent_scope
          for child_label_object in label_objects:
            if child_label_object is not parent_label_object \
              and child_label_object in data_manager_data_sources \
              and child_label_object in label_objects:
              data_manager_sources.remove(child_label_object)
              phil_object = _get_last_object(child_label_object)
              if phil_object.full_path() == 'data_manager.%s.labels.name' % datatype:
                phil_object = phil_object.primary_parent_scope
              pps.objects.append(phil_object)
              phil_object.primary_parent_scope = pps

        # finally combine labels into one scope
        label_phil = label_objects[0]
        if len(label_objects) == len(datatypes):
          for label_object in label_objects[1:]:
            label_phil.adopt_scope(label_object)
            if label_object in data_manager_sources:
              data_manager_sources.remove(label_object)

        print('', file=self.logger)

        print('  Combined labels PHIL', file=self.logger)
        print('  --------------------', file=self.logger)
        tmp_working_phil = self.data_manager.master_phil.fetch_diff(label_phil)
        print(tmp_working_phil.as_str(prefix='    '), file=self.logger)
        print('', file=self.logger)

      # set model types for each model
      # if model type is specified, a model type needs to be specified
      # for each model even if it is the default.
      if len(type_objects) > 0:
        model_names = self.data_manager.get_model_names()
        if len(type_objects) != len(model_names):
          raise Sorry('Please specify exactly one "model.type" for each model.')
        print('  Matching model type PHIL:', file=self.logger)
        print('  -------------------------', file=self.logger)
        for model_name, type_object in zip(model_names, type_objects):
          type_phil = self.data_manager.master_phil.fetch(source=type_object)
          type_extract = type_phil.extract()
          model_extract = type_extract.data_manager.model[0]
          model_extract.file = model_name
          print('    %s %s' % (model_extract.file, model_extract.type), file=self.logger)
          new_type_object = self.data_manager.master_phil.format(python_object=type_extract)
          for object in new_type_object.objects[0].objects:
            if object.name == 'model':
              e = object.extract()
              if e.file == model_name:
                type_object.objects[0].objects = [object]
        print('', file=self.logger)

    if self.namespace.overwrite:  # override overwrite if True
      sources.append(iotbx.phil.parse('output.overwrite=True'))

    # process program parameters
    skip_incompatible_objects = False
    if self.namespace.diff_params:
      skip_incompatible_objects = True
    if len(data_sources) + len(sources) > 0:
      self.working_phil, more_unused_phil = self.master_phil.fetch(
        sources=data_sources + sources,
        track_unused_definitions=True,
        skip_incompatible_objects=skip_incompatible_objects)
      self.unused_phil.extend(more_unused_phil)
    elif self.working_phil is None:
      self.working_phil = self.master_phil.fetch()

    # process DataManager parameters
    if len(data_manager_data_sources) + len(data_manager_sources) > 0:
      diff_phil, more_unused_phil = self.data_manager.master_phil.fetch_diff(
        sources=data_manager_data_sources + data_manager_sources, track_unused_definitions=True)
      self.unused_phil.extend(more_unused_phil)
      # load remaining files and final fmodel parameters
      diff_phil = self._update_phil_paths_to_cwd(diff_phil)
      self.data_manager.load_phil_scope(diff_phil, process_files=not self.namespace.diff_params)
      self.data_manager_diff = diff_phil

    # show unrecognized parameters and abort
    self.raise_Sorry_for_unused_phil()

    # process input phil for file/directory defintions and add to DataManager
    # Note: if a PHIL file is input as a PHIL parameter, the contents of the
    # file will NOT be parsed and validated. The PHIL file should be provided
    # as a command-line argument. This is mostly for finding data files
    # defined by PHIL parameters that should be added to the DataManager
    diff_phil = self.master_phil.fetch_diff(self.working_phil)
    diff_phil = self._update_phil_paths_to_cwd(diff_phil)
    paths = self.check_phil_for_paths(diff_phil)
    if len(paths) > 0:
      files = set()
      dirs = set()
      for path in paths:
        if path is not None:
          if os.path.isfile(path):
            files.add(path)
          elif os.path.isdir(path):
            dirs.add(path)
      if self.parse_files:
        self.process_files(files, message='Processing files from PHIL:')
      if self.parse_dir:
        self.process_dir(dirs, message='Processing directories from PHIL:')

    if not printed_something:
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
    paths = []
    if phil_scope.is_definition:
      if phil_scope.type.phil_type == 'path':
        if phil_scope.style is not None and 'new_file' in phil_scope.style:
          pass
        else:
          paths.append(phil_scope.extract())
    elif phil_scope.is_scope:
      for phil_object in phil_scope.objects:
        paths.extend(self.check_phil_for_paths(phil_object))
    return paths

  # ---------------------------------------------------------------------------
  def _update_phil_paths_to_cwd(self, phil=None):
    '''
    Convenience function for modifying a phil scope so that file paths
    are reset to use the current working directory if the filename exists.
    '''
    if self.namespace.check_current_dir:
      diff_phil = self.data_manager.master_phil.fetch_diff(phil)
      paths = self.check_phil_for_paths(diff_phil)
      phil_str = diff_phil.as_str()
      for path in paths:
        new_path = os.path.join(os.getcwd(), os.path.basename(path))
        if not os.path.isfile(path) and os.path.isfile(new_path):
          phil_str = phil_str.replace(path, new_path)
      phil = iotbx.phil.parse(phil_str, process_includes=True)
    return phil

  # ---------------------------------------------------------------------------
  def post_process(self):
    '''
    Post processing of inputs after all arguments are parsed
    '''

    working_phil_extract = self.working_phil.extract()

  # ---------------------------------------------------------------------------
  def show_phil_summary(self):
    '''
    Show final, modified PHIL parameters after all processing is complete
    Also, write phil scopes based on command-line flags
    '''

    overwrite = (self.namespace.overwrite or \
                 self.working_phil.extract().output.overwrite)

    # check for any remaining unknown arguments
    if len(self.namespace.unknown) > 0:
      error_message = 'The following arguments are not recognized:\n'
      for value in self.namespace.unknown:
        error_message += '  %s\n' % value
      raise Sorry(error_message)

    # get differences
    try:
      data_diff = self.data_manager.master_phil.fetch_diff(
        self.data_manager.export_phil_scope())
      # keep original DataManager scope when using --diff-params
      if self.namespace.diff_params and self.data_manager_diff is not None:
        data_diff = self.data_manager_diff
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
    if is_different:
      data_diff.show(prefix='  ', out=self.logger)
      phil_diff.show(prefix='  ', out=self.logger)
    else:
      print('  All parameters are set to their defaults', file=self.logger)
    print('', file=self.logger)

    # write scopes if requested
    if self.namespace.write_data or self.namespace.write_defaults \
      or self.namespace.write_modified or self.namespace.write_all:
      print('Writing program PHIL file(s):', file=self.logger)

    if self.data_manager.get_default_output_filename() is None:
      self.data_manager.set_default_output_filename('cctbx_program')

    # write DataManager scope
    if self.namespace.write_data:
      if data_is_different:
        self.data_manager.write_phil_file(
          self.data_manager.export_phil_scope().as_str(),
          filename=self.data_filename,
          overwrite=overwrite)
        print('  Input file PHIL written to %s.' % self.data_filename,
              file=self.logger)
      else:
        print('  No input file PHIL to write', file=self.logger)

    # write all default parameters
    if self.namespace.write_defaults:
      self.data_manager.write_phil_file(self.master_phil.as_str(expert_level=3),
        filename=self.defaults_filename, overwrite=overwrite)
      print('  Default PHIL parameters written to %s.' % self.defaults_filename,
        file=self.logger)

    # write differences
    if self.namespace.write_modified or self.namespace.diff_params:
      if is_different:
        ow = overwrite or self.namespace.diff_params
        self.data_manager.write_phil_file(
          data_diff.as_str() + phil_diff.as_str(), filename=self.modified_filename,
          overwrite=ow)
        print('  Modified PHIL parameters written to %s.' %
              self.modified_filename, file=self.logger)
      else:
        print('  No PHIL modifications to write', file=self.logger)
      if self.namespace.diff_params:
        sys.exit(0)

    # write all parameters (DataManager + Program)
    if self.namespace.write_all:
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
    program_citations = []
    if self.program_class.citations is not None:
      class_citations = citations.master_citation_phil.fetch(
        source=self.program_class.citations).extract()
      for citation in class_citations.citation:
        program_citations.append(citation)
    for article_id in self.program_class.known_article_ids:
      citation = citations.citations_db.get(article_id)
      if citation is not None:
        program_citations.append(citation)
      else:
        raise Sorry('"%s" not found citations database' % article_id)

    # show program-specific citations and general citation for CCTBX
    if len(program_citations) > 0:
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
def run_pyside_check():
  '''
  Function for checking if PySide2 is available
  '''
  try:
    import PySide2
  except ImportError:
    msg = '''
------------------------------------------------------------------------
To run this GUI, PySide2 is required. To install with conda run

  conda install pyside2

or with pip run

  pip install pyside2

------------------------------------------------------------------------
'''
    raise Sorry(msg)

# =============================================================================
def run_program(program_class=None, parser_class=CCTBXParser, custom_process_arguments=None,
                unused_phil_raises_sorry=True, args=None, json=False, logger=None,
                hide_parsing_output=False, allow_default_args=False):
  '''
  Function for running programs using CCTBXParser and the program template

  Parameters
  ----------
  program_class: ProgramTemplate
    The class defining the program. It must be a subclass of ProgramTemplate
  parser_class: CCTBXParser
    The parser class to use for parsing. It must be the CCTBXParser or a subclass
  custom_process_arguments: function(parser)
    Custom function to parse unknown arguments (optional)
  unused_phil_raises_sorry: bool
    If False, any unused PHIL parameters are kept for parsing later
  args: list
    List of command-line arguments (optional)
  json: bool
    If True, get_results_as_JSON is called for the return value instead of get_results
  logger: multi_out
    For logging output (optional)
  hide_parsing_output: bool
    If True, hides the output from parsing the command-line arguments

  Returns
  -------
    Whatever is returned from program_class.get_results() or program_class.get_results_as_JSON()
  '''

  assert program_class is not None

  if args is None:
    args = sys.argv[1:]
  elif allow_default_args and type(args) in [list, tuple]:
    args = sys.argv[1:] + args

  # start profiling
  pr = None
  if '--profile' in args:
    import cProfile
    pr = cProfile.Profile()
    pr.enable()

  # keep output in quiet mode
  if '--quiet' in args:
    logger = multi_out()
    logger.register('parser_log', StringIO())

  # create logger
  if logger is None:
    logger = multi_out()
    logger.register('stdout', sys.stdout)
    logger.register('parser_log', StringIO())

  program_logger = logger
  if hide_parsing_output:
    logger = multi_out()
    logger.register('stderr', sys.stderr)
    logger.register('parser_log', StringIO())

  # start timer
  t = show_times(out=logger)

  # create parser
  parser = parser_class(program_class=program_class,
                        custom_process_arguments=custom_process_arguments,
                        unused_phil_raises_sorry=unused_phil_raises_sorry,
                        logger=logger)
  namespace = parser.parse_args(args)

  # start program
  if namespace.dry_run:
    print('Starting dry run', file=logger)
  else:
    print('Starting job', file=logger)
  print('='*79, file=logger)
  task = program_class(parser.data_manager, parser.working_phil.extract(),
                       master_phil=parser.master_phil,
                       logger=program_logger)

  # validate inputs
  task.validate()

  # stop if dry_run is set
  if namespace.dry_run:
    print('\nArguments have been validated by the program.\n', file=logger)
    print('='*79, file=logger)
    return

  # run program
  task.run()

  # clean up (optional)
  task.clean_up()

  # dump profiling stats
  if pr is not None:
    pr.disable()
    pr.dump_stats('profile.out')

  # output JSON
  if namespace.json or namespace.json_filename:
    result = task.get_results_as_JSON()
    if result is not None:
      json_filename = parser.json_filename
      if namespace.json_filename is not None:
        json_filename = namespace.json_filename
        if not json_filename.endswith('.json'):
          json_filename += '.json'
      with open(json_filename, 'w') as f:
        f.write(result)
    else:
      print('', file=logger)
      print('!'*79, file=logger)
      print('WARNING: The get_results_as_JSON function has not been defined for this program', file=logger)
      print('!'*79, file=logger)

  # stop timer
  print('', file=logger)
  print('='*79, file=logger)
  print('Job complete', file=logger)
  t()

  # clean up file for quiet mode
  if namespace.quiet:
    logger.close()

  if json:
    result = task.get_results_as_JSON()
  else:
    result = task.get_results()

  return result

# =============================================================================
def get_program_params(run):
  """Tool to get parameters object for a program that runs with
     the program template.
  params: run:  the program template object
  returns: parameters for this program as set up by the program template
  Get the run something like this way:
    from phenix.programs import map_to_model as run
  """

  parser = CCTBXParser(program_class=run.Program,
                       logger=null_out())
  _ = parser.parse_args([], skip_help = True)
  return parser.working_phil.extract()

# =============================================================================
# end
