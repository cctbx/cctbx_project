from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.merge

from xfel.merging.application.phil.phil import phil_scope

help_message = '''
Merge xfel data.
'''

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    # The script usage
    import libtbx.load_env
    self.usage = "usage: %s [options] [param.phil] " % libtbx.env.dispatcher_name
    self.parser = None

    '''Initialise the script.'''
    from dials.util.options import OptionParser
    # Create the parser
    self.parser = OptionParser(
      usage=self.usage,
      phil=phil_scope,
      epilog=help_message)

    # Parse the command line. quick_parse is required for MPI compatibility
    params, options = self.parser.parse_args(show_diff_phil=True,quick_parse=True)
    self.params = params
    self.options = options

  def run(self):
    from xfel.merging import application
    import importlib

    # Create the workers using the factories
    workers = []
    for step in ['input']:
      factory = importlib.import_module('xfel.merging.application.'+step+'.factory')
      workers.extend(factory.factory.from_parameters(self.params))

    # Perform phil validation up front
    for worker in workers:
      worker.validate()

    # Do the work
    experiments = reflections = None
    while(workers):
      worker = workers.pop(0)
      experiments, reflections = worker.run(experiments, reflections)

    print ('Done')

if __name__ == '__main__':
  script = Script()
  result = script.run()
  print ("OK")
