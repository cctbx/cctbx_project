
import logging
from libtbx.mpi4py import MPI

# The below code is from https://gist.github.com/sixy6e/ed35ea88ba0627e0f7dfdf115a3bf4d1


class MPIIOStream(object):

  def __init__(self, filename, comm, mode):

    self._file = MPI.File.Open(comm, filename, mode)
    self._file.Set_atomicity(True)

  def write(self, msg):
    # if for some reason we don't have a unicode string...
    try:
      msg = msg.encode()
    except AttributeError:
      pass
    self._file.Write_shared(msg)

  def sync(self):
    self._file.Sync()

  def close(self):
    self.sync()
    self._file.Close()


class MPIFileHandler(logging.StreamHandler):

  def __init__(self, filename,
               mode=MPI.MODE_WRONLY|MPI.MODE_CREATE, comm=MPI.COMM_WORLD):
    self.filename = filename
    self.mode = mode
    self.comm = comm

    super(MPIFileHandler, self).__init__(self._open())

  def _open(self):
    stream = MPIIOStream(self.filename, self.comm, self.mode)
    return stream

  def close(self):
    if self.stream:
      self.stream.close()
      self.stream = None

  def emit(self, record):
    """
    Emit a record.
    We have to override emit, as the logging.StreamHandler has 2 calls
    to 'write'. The first for the message, and the second for the
    terminator. This posed a problem for mpi, where a second process
    could call 'write' in between these two calls and create a
    conjoined log message.
    """
    msg = self.format(record)
    self.stream.write('{}{}'.format(msg, self.terminator))
    self.flush()


def main():
  """
  A sample test run.
  """
  comm = MPI.COMM_WORLD
  logger = logging.getLogger("node[%i]"%comm.rank)
  # logger = logging.getLogger("func-status") # another name example

  logger.setLevel(logging.DEBUG)

  mpi_handler = MPIFileHandler("test.log")
  formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
  mpi_handler.setFormatter(formatter)
  logger.addHandler(mpi_handler)

  # sample log levels
  logger.debug('debug message')
  logger.info('info message')
  logger.warning('warn message')
  logger.error('error message')
  logger.critical('critical message')


if __name__ == "__main__":
  main()
