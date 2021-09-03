
import os
import logging
import socket
HOST = socket.gethostname()
from libtbx.mpi4py import MPI
COMM = MPI.COMM_WORLD
if not hasattr(COMM, "rank"):
    COMM.rank = 0
    COMM.size = 1

# gather all hostnames and create sub-communicators for all processes on a given host
unique_hosts = COMM.gather(HOST)
HOST_MAP = None
if COMM.rank==0:
    HOST_MAP = {HOST:i for i,HOST in enumerate(set(unique_hosts))}
HOST_MAP = COMM.bcast(HOST_MAP)
HOST_COMM = COMM.Split(color=HOST_MAP[HOST])


LEVELS = {"low": logging.WARNING, "normal": logging.INFO, "high": logging.DEBUG}
DETAILED_FORMAT = 'RANK%d:%s | ' % (COMM.rank, HOST) + '%(asctime)s | %(filename)s:%(funcName)s >>  %(message)s'
SIMPLE_FORMAT = "RANK%d" % COMM.rank + " | %(message)s"

from simtbx.diffBragg import utils


def _make_logger(loggername, filename=None, formatter=None, level=None, overwrite=False,
                propagate=False):
    """
    :param loggername: name for logger used throughout application
    :param filename: if logging to a file, this should be set
    :param formatter: logging Formatter object
    :param level: logging level number
    :param overwrite: whether to overwrite the existing logfile, if logging to files
    :param propagate: propagate messages to root logger, can result in duplicate messages
    :return: return the logger instance, however, once created here, the instance is accessible
            from other modules via the command logging.getLogger(loggername)
    """
    logger = logging.getLogger(loggername)
    if filename is not None:
        if COMM.rank == 0 and overwrite and os.path.exists(filename):
            os.remove(filename)
        COMM.barrier()
        handler = MPIFileHandler(filename, comm=HOST_COMM)
    else:
        handler = logging.StreamHandler()
    if formatter is not None:
        handler.setFormatter(formatter)
    if level is not None:
        logger.setLevel(level)
    logger.addHandler(handler)
    logger.propagate = propagate
    return logger


def setup_logging_from_params(params):
    """params: PHIL params, see simtbx.diffBragg.hopper phil string"""
    if params.logging.logfiles:
        if COMM.rank == 0:
            utils.safe_makedirs(params.outdir)
        COMM.barrier()
        main_level = LEVELS[params.logging.logfiles_level]
        main_logger = _make_logger("main",
                                  os.path.join(params.outdir, HOST+"-"+params.logging.logname),
                                  level=main_level,
                                  overwrite=params.logging.overwrite,
                                  formatter=logging.Formatter(DETAILED_FORMAT))

        _make_logger("profile",
                    os.path.join(params.outdir, HOST+"-"+params.profile_name),
                    level=logging.INFO,
                    overwrite=params.logging.overwrite,
                    formatter=logging.Formatter(SIMPLE_FORMAT))

        # for convenience we add a console logger for rank0 so we can optionally see output
        # even when logging to files
        if COMM.rank == 0:
            console = logging.StreamHandler()
            console.setFormatter(logging.Formatter(SIMPLE_FORMAT))
            console.setLevel(LEVELS[params.logging.rank0_level])
            main_logger.addHandler(console)

    else:
        if COMM.rank == 0:
            level = LEVELS[params.logging.rank0_level]
        else:
            level = LEVELS[params.logging.other_ranks_level]
        _make_logger("main", level=level, formatter=logging.Formatter(SIMPLE_FORMAT))
        _make_logger("profile", level=level, formatter=logging.Formatter(SIMPLE_FORMAT))


class MPIIOStream(object):
    # The class is from https://gist.github.com/sixy6e/ed35ea88ba0627e0f7dfdf115a3bf4d1

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
    # The class is from https://gist.github.com/sixy6e/ed35ea88ba0627e0f7dfdf115a3bf4d1

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
    A sample test run from https://gist.github.com/sixy6e/ed35ea88ba0627e0f7dfdf115a3bf4d1
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
