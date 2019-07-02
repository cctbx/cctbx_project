'''
Standard Program Template for CCTBX Programs

The "program" is the actual task to be performed without any user interfaces.
The user interfaces (command-line and graphical) build the data_manager and
params objects for the program. The "data_manager"" handles all file input
and "params" handles all the program settings. These two objects should have
all relevant information for the program to run.

The required functions break up the calling order into discrete phases

- constructor: minimal set up
- validate: check that the inputs (files and parameters) are valid and consistent
- run: run the actual task
- get_results: return the desired output from the program

The optional functions provide some extra tweaking

- custom_init: called at the end of the constructor, additional initialization
- clean_up: if temporary files are written in the course of running the program,
            this step should remove those files.

Additional functions and class attributes can be defined for doing the actual
task, but the above functions define a consistent interface.

More documentation to come
'''
from __future__ import absolute_import, division, print_function

import libtbx.phil

from libtbx import citations
from libtbx.utils import multi_out

# =============================================================================
class ProgramTemplate(object):
  # Class variables for customizing program

  # description of the program
  description = '''
Program Description
'''

  # datatypes for program
  # see libtbx/data_manager/<datatype>.py for list of supported datatypes
  # default datatypes are set in libtbx/data_manager/__init__.py
  datatypes = None

  # master PHIL string for the program (required)
  master_phil_str = '''
# example
program {
  parameter = None
    .type = bool
}
'''

  # unique citations for the program. list of citation phil extract objects
  # see libtbx/citations.py for the PHIL format.
  citations = None

  # common citations used by the program that exist in libtbx/citations.params
  # list of article_id strings, e.g. ["polder", "elbow"]).
  known_article_ids = list()

  # text shown at the end of the command-line program
  epilog = '''
For additional help, you can contact the developers at cctbxbb@phenix-online.org
or https://github.com/cctbx/cctbx_project

'''

  # ---------------------------------------------------------------------------
  # Reserved phil scope for output
  # this will be automatically added to the master_phil_str
  # you should add your own output phil scope, but these parameters will be
  # automatically added, so no need to redefine
  output_phil_str = '''
output {
  prefix = None
    .type = str
    .help = Prefix string added to automatically generated output filenames
  suffix = None
    .type = str
    .help = Suffix string added to automatically generated output filenames
  serial = 0
    .type = int
    .help = Serial number added to automatically generated output filenames
  overwrite = False
    .type = bool
    .help = Overwrite files when set to True
}
'''

  # ---------------------------------------------------------------------------
  # Advanced features

  # PHIL converters (in a list) for additional PHIL types
  phil_converters = list()

  # ---------------------------------------------------------------------------
  # Function for showing default citation for template
  @staticmethod
  def show_template_citation(text_width=80, logger=None,
                             citation_format='default'):
    assert(logger is not None)

    print('\nGeneral citation for CCTBX:', file=logger)
    print('-'*text_width, file=logger)
    print('', file=logger)
    citations.show_citation(citations.citations_db['cctbx'], out=logger,
                            format=citation_format)

  # ---------------------------------------------------------------------------
  def __init__(self, data_manager, params, master_phil=None, logger=None):
    '''
    Common constructor for all programs

    This is supposed to be lightweight. Custom initialization, if necessary,
    should be handled by the custom_init function. Developers should not need to
    override this function.

    Parameters
    ----------
    data_manager :
      An instance of the DataManager (libtbx/data_manager.py) class containing
      data structures from file input
    params :
      An instance of PHIL
    logger :
      Standard Python logger (from logging module), optional. A logger will be
      created if it is not provided.

    '''

    self.data_manager = data_manager
    self.master_phil = master_phil
    self.params = params
    self.logger = logger

    if (self.logger is None):
      self.logger = multi_out()

    # master_phil should be provided by CCTBXParser or GUI because of
    # potential PHIL extensions
    if (self.master_phil is None):
      self.master_phil = libtbx.phil.parse(
        self.master_phil_str, process_includes=True)

    self.custom_init()

    # set DataManager defaults
    if self.data_manager is not None:
      self.data_manager.set_default_output_filename(
        self.get_default_output_filename())
      try:
        self.data_manager.set_overwrite(self.params.output.overwrite)
      except AttributeError:
        pass

  def header(self, text):
    print("-"*79, file=self.logger)
    print(text, file=self.logger)
    print("*"*len(text), file=self.logger)

  # ---------------------------------------------------------------------------
  def custom_init(self):
    '''
    Optional initialization step

    Developers should override this function if additional initialization is
    needed. There should be no arguments because all necessary information
    should be in self.data_manager (file input) and self.params (phil parameters)

    Parameters
    ----------
    None
    '''
    pass

  # ---------------------------------------------------------------------------
  def validate(self):
    '''

    '''
    raise NotImplementedError('The "validate" function is required.')

  # ---------------------------------------------------------------------------
  def run(self):
    '''

    '''
    raise NotImplementedError('The "run" function is required.')

  # ---------------------------------------------------------------------------
  def clean_up(self):
    '''

    '''
    pass

  # ---------------------------------------------------------------------------
  def get_results(self):
    '''

    '''
    return None

  # ---------------------------------------------------------------------------
  def get_results_as_JSON(self):
    '''

    '''
    return None

  # ---------------------------------------------------------------------------
  def get_program_phil(self, diff=False):
    '''
    Function for getting the PHIL extract of the Program

    Parameters
    ----------
    diff: bool
      When set to True, only the differences from the master PHIL are returned

    Returns
    -------
    params: libtbx.phil.scope
    '''
    working_phil = self.master_phil.format(python_object=self.params)
    if diff:
      working_phil = self.master_phil.fetch_diff(working_phil)
    return working_phil

  # ---------------------------------------------------------------------------
  def get_data_phil(self, diff=False):
    '''
    Function for getting the PHIL scope from the DataManager

    Parameters
    ----------
    diff: bool
      When set to True, only the differences from the master PHIL are returned

    Returns
    -------
    params: libtbx.phil.scope
    '''
    if self.data_manager is None:
      return libtbx.phil.parse('')
    working_phil = self.data_manager.export_phil_scope()
    if diff:
      working_phil = self.data_manager.master_phil.fetch_diff(working_phil)
    return working_phil

  # ---------------------------------------------------------------------------
  def get_program_extract(self, diff=False):
    '''
    Function for getting the PHIL extract of the Program

    Parameters
    ----------
    diff: bool
      When set to True, only the differences from the master PHIL are returned

    Returns
    -------
    params: libtbx.phil.scope_extract
    '''
    return self.get_program_phil(diff=diff).extract()

  # ---------------------------------------------------------------------------
  def get_data_extract(self, diff=False):
    '''
    Function for getting the PHIL extract from the DataManager

    Parameters
    ----------
    diff: bool
      When set to True, only the differences from the master PHIL are returned

    Returns
    -------
    params: libtbx.phil.scope_extract
    '''
    return self.get_data_phil(diff=diff).extract()

  # ---------------------------------------------------------------------------
  def get_program_phil_str(self, diff=False):
    '''
    Function for getting the PHIL string of the Program

    Parameters
    ----------
    diff: bool
      When set to True, only the differences from the master PHIL are returned

    Returns
    -------
    params: str
    '''
    return self.get_program_phil(diff=diff).as_str()

  # ---------------------------------------------------------------------------
  def get_data_phil_str(self, diff=False):
    '''
    Function for getting the PHIL string from the DataManager

    Parameters
    ----------
    diff: bool
      When set to True, only the differences from the master PHIL are returned

    Returns
    -------
    params: str
    '''
    return self.get_data_phil(diff=diff).as_str()

  # ---------------------------------------------------------------------------
  def get_full_phil_str(self, diff=False):
    '''
    Function for getting the full PHIL string of the DataManager and Program

    Parameters
    ----------
    diff: bool
      When set to True, only the differences from the master PHIL are returned

    Returns
    -------
    params: str
    '''
    return self.get_data_phil_str(diff=diff) + self.get_program_phil_str(diff=diff)

  # ---------------------------------------------------------------------------
  def get_default_output_filename(self):
    '''
    Given the output.prefix, output.suffix, and output.serial PHIL parameters,
    return the default output filename

    Parameters
    ----------
    None

    Returns
    -------
    filename: str
      The default output filename without a file extension
    '''

    filename = 'cctbx_program'
    if hasattr(self.params, 'output'):
      if getattr(self.params.output, 'prefix', None) is not None:
        filename = self.params.output.prefix
      if getattr(self.params.output, 'suffix', None) is not None:
        filename += '{suffix}'.format(suffix=self.params.output.suffix)
      if getattr(self.params.output, 'serial', None) is not None:
        filename += '_{serial:03d}'.format(serial=self.params.output.serial)

    return filename

# =============================================================================
