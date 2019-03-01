from __future__ import absolute_import, division, print_function
import libtbx
import libtbx.env_config
import os
libtbx.env = libtbx.env_config.unpickle()
libtbx.env.set_os_environ_all_dist()
libtbx.env.dispatcher_name = os.environ.get("LIBTBX_DISPATCHER_NAME")
if not libtbx.env.dispatcher_name:
  # Attempt to identify dispatcher name if LIBTBX_DISPATCHER_NAME is not set
  try:
    import inspect
    _frame = inspect.stack()[1]
    _module = inspect.getmodule(_frame[0])
    _sourcefile = os.path.realpath(_module.__file__)
    for _dist_path in [os.path.realpath(_d) for _d in libtbx.env.dist_paths()]:
      if _sourcefile.startswith(_dist_path):
        _command = '.'.join([ _dist_path.split(os.path.sep)[-1] ] + list(filter(lambda x: x and x != 'command_line', _sourcefile[len(_dist_path):].split(os.path.sep))))
        for _ext in ('.pyo', '.pyc', '.py'):
          if _command.endswith(_ext):
            _command = _command[:-len(_ext)]
        libtbx.env.dispatcher_name = _command
        break
  except Exception: # intentional
    pass # If anything goes wrong - give up silently

libtbx.env.full_testing = os.environ.get("LIBTBX_FULL_TESTING") is not None
