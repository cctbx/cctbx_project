from __future__ import division, print_function
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

import logging

# =============================================================================
class ProgramTemplate(object):

  def __init__(self, data_manager, params, logger=None):
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
    self.params = params
    self.logger = logger

    if (self.logger is None):
      self.logger = logging.getLogger('program')

    self.custom_init()

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
# =============================================================================
