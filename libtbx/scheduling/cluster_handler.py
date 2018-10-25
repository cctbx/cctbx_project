from __future__ import absolute_import, division, print_function

from libtbx.queuing_system_utils import processing
from libtbx.queuing_system_utils.processing import errors

from libtbx.scheduling import SetupError

class Job(processing.Job):
  """
  Customisation of the processing.Job class that translates errors into SetupError
  """

  def __exception_converting_call(self, name):

    try:
      return getattr( super( Job, self ), name )()

    except errors.BatchQueueError as e:
      raise SetupError("Queue error: %s" % e)


  @property
  def exitcode(self):

    try:
      return super( Job, self ).exitcode

    except errors.BatchQueueError as e:
      raise SetupError("Queue error: %s" % e)


  def start(self):

    return self.__exception_converting_call( name = "start" )


  def is_alive(self):

    return self.__exception_converting_call( name = "is_alive" )


  def join(self):

    return self.__exception_converting_call( name = "join" )


  def terminate(self):

    return self.__exception_converting_call( name = "terminate" )


class JobFactory(object):
  """
  Creator for Job objects

  Note this bypasses the Qinteface's own Job method
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

    try:
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

    except errors.BatchQueueError as e:
      raise SetupError("Queue error: %s" % e)


  def __call__(self, target, args = (), kwargs = {}):

    return Job(
      qinterface = self.qinterface,
      target = target,
      args = args,
      kwargs = kwargs,
      )
