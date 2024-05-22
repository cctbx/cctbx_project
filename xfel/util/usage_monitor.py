from contextlib import ContextDecorator
from enum import Enum
import logging
import platform
import psutil
import subprocess
import threading
import time
from typing import List

from libtbx.mpi4py import MPI
comm = MPI.COMM_WORLD
log_path = f'usage_{comm.rank}.log'
logging.basicConfig(filename=log_path, encoding='utf-8', level=logging.ERROR)
log = logging.getLogger("cctbx.usage")
log.setLevel(logging.DEBUG)


class UsageMonitor(ContextDecorator):
  """This class can be used both as context manager (with) and @decorator"""

  class Detail(Enum):
    off = 'off'
    single = 'single'
    node = 'node'
    rank = 'rank'

  def __init__(self, detail: str = 'node', every: float = 1.0) -> None:
    self.detail = self.Detail(detail)
    self.every = every
    self._gpu_count = None

  def __enter__(self) -> None:
    if comm.rank == 0:
      self.start_logging_deamon()

  def __exit__(self, exc_type, exc_val, exc_tb):
    pass

  @property
  def cpu_memory(self) -> float:
    mem = psutil.virtual_memory()
    return mem.used / mem.total

  @property
  def cpu_usage(self) -> float:
    return psutil.cpu_percent(interval=None)

  @property
  def gpus_memory(self) -> List[float]:
    smi = subprocess.run(['nvidia-smi', '--query-gpu=utilization.memory',
                          '--format=csv,noheader,nounits'],
                         stdout=subprocess.PIPE)
    return [float(u) for u in smi.stdout.decode('utf-8').strip().split()]

  @property
  def gpus_usage(self) -> List[float]:
    smi = subprocess.run(['nvidia-smi', '--query-gpu=utilization.gpu',
                          '--format=csv,noheader,nounits'], stdout=subprocess.PIPE)
    return [float(u) for u in smi.stdout.decode('utf-8').strip().split()]

  @property
  def gpu_count(self):
    if not self._gpu_count:
      self._gpu_count = len(self.gpus_usage)
    return self._gpu_count

  @property
  def node(self) -> str:
    return platform.node()

  def log_usage(self) -> None:
    cpu_usage = self.cpu_usage
    cpu_memory = self.cpu_memory
    gpu_usage = self.gpus_usage
    gpu_memory = self.gpus_memory
    log.info(msg=f'{cpu_usage=}, {cpu_memory=}, {gpu_usage=}, {gpu_memory=}.')

  def logging_deamon(self) -> None:
    """Warning: call as threading daemon only, otherwise will never stop"""
    while True:
      threading.Thread(target=self.log_usage, args=(self, ))
      time.sleep(self.every)

  def start_logging_deamon(self) -> threading.Thread:
    return threading.Thread(target=self.logging_deamon,
                            args=(self,), daemon=True)



