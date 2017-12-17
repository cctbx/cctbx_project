---------------------------------------------------------------------
libtbx.program_template: Standard Program Template for CCTBX Programs
---------------------------------------------------------------------

Program Template
================

The "program" is the actual task to be performed without any user interfaces.
The user interfaces (command-line and graphical) build the "data_manager" and
"params" objects for the program. The "data_manager" handles all file input and
"params" handles all the program settings.

The required methods break up the calling order into discrete phases

- constructor: minimal set up
- validate: check that the inputs (files and parameters) are valid and consistent
- run: run the actual task
- get_results: return the desired output from the program

The optional methods provide some extra tweaking

- custom_init: called at the end of the constructor, additional initialization
- clean_up: if temporary files are written in the course of running the program,
            this step should remove those files.

.. autoclass:: libtbx.program_template.ProgramTemplate
               :members: __init__, custom_init, validate, run, clean_up,
                         get_results
