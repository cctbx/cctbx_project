from __future__ import division

from libtbx.queuing_system_utils import processing

class JobFactory(object):
  """
  Creator for Queue.Job objects
  """

  def __init__(self,
    platform,
    name = "libtbx_python",
    command = None,
    asynchronous = True,
    use_target_file = True,
    include = None,
    poller = None,
    handler = None,
    save_error = True,
    ):

    from libtbx.queuing_system_utils.processing import transfer
    self.qinterface = processing.INTERFACE_FOR[ platform ][0](
      name = name,
      command = command,
      asynchronous = asynchronous,
      input = transfer.TemporaryFile if use_target_file else transfer.Stdin,
      include = include,
      poller = poller,
      handler = handler,
      save_error = save_error,
      display_stderr = False,
      )


  def __call__(self, target, args = (), kwargs = {}):

    return self.qinterface.Job( target = target, args = args, kwargs = kwargs)
