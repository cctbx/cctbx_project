from contextlib import ContextDecorator
from enum import Enum
from functools import cached_property
import logging
import platform
import psutil
import subprocess
import threading
import time
from typing import List

from libtbx.mpi4py import MPI
comm = MPI.COMM_WORLD
node = platform.node()
node_ranks = [r for n, r in comm.allgather((node, comm.rank)) if n == node]


class UsageMonitor(ContextDecorator):
  """This class can be used both as context manager (with) and @decorator"""

  class Detail(Enum):
    rank = 'rank'
    node = 'node'
    single = 'single'
    none = 'none'

  def __init__(self, detail: str = 'node', period: float = 1.0) -> None:
    self.detail = self.Detail(detail)
    self.period: float = period
    self.log: logging.Logger = self.configure_logger()
    self._daemon = None

  def __enter__(self) -> None:
    if comm.rank == 0:
      self.start_logging_daemon()

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
  def daemon(self) -> threading.Thread:
    return self._daemon

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

  @cached_property
  def gpu_count(self) -> int:
    return len(self.gpus_usage)

  @cached_property
  def is_logging(self) -> bool:
    return {
      self.Detail.rank: True,
      self.Detail.node: comm.rank == min(node_ranks),
      self.Detail.single: comm.rank == 0,
      self.Detail.none: False
    }[self.detail]

  @cached_property
  def log_path(self) -> str:
    return {
      self.Detail.rank: f'usage_{comm.rank}.log',
      self.Detail.node: f'usage_{node}.log',
      self.Detail.single: 'usage.log',
      self.Detail.none: ValueError
    }[self.detail]

  def configure_logger(self) -> logging.Logger:
    if self.is_logging:
      logging.basicConfig(filename=self.log_path, encoding='utf-8',
                          level=logging.ERROR)
      log = logging.getLogger("cctbx.usage")
      log.setLevel(logging.DEBUG)
    else:
      log = logging.getLogger("cctbx.usage")
      log.setLevel(logging.CRITICAL + 1)
      log.addHandler(logging.NullHandler())
    return log

  def log_current_usage(self) -> None:
    cpu_usage = self.cpu_usage
    cpu_memory = self.cpu_memory
    gpu_usage = self.gpus_usage
    gpu_memory = self.gpus_memory
    self.log.info(msg=f'{cpu_usage=}, {cpu_memory=}, {gpu_usage=}, {gpu_memory=}.')

  def log_usage_every_period(self) -> None:
    """Warning: call as threading daemon only, otherwise will never stop"""
    while True:
      threading.Thread(target=self.log_current_usage, args=(self,)).start()
      time.sleep(self.period)

  def start_logging_daemon(self) -> None:
    self._daemon = threading.Thread(target=self.log_usage_every_period,
                                    args=(self,), daemon=True)
    self.daemon.start()
