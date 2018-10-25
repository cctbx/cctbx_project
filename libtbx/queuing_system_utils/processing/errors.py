from __future__ import absolute_import, division, print_function

class ProcessingError(Exception):
  """
  Base class
  """


# Job control and information
class BatchQueueError(ProcessingError):
  """
  Error with the underlying batch queuing system
  """


class ExecutableError(BatchQueueError):
  """
  There was a problem with executing the requested program
  """


class AbnormalExitError(BatchQueueError):
  """
  The requested program returned an error status code
  """


# Output processing errors
class ExtractionError(ProcessingError):
  """
  Output from the underlying program is not what is expected
  """
