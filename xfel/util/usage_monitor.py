from collections import UserDict
from contextlib import ContextDecorator
from dataclasses import astuple, dataclass
from datetime import datetime, timedelta
from enum import Enum
from functools import cached_property
from glob import glob
import itertools
import logging
from pathlib import Path
import platform
import re
import os
from statistics import mean
import subprocess
import threading
import time
from typing import List, Type, Union
PathLike = Union[str, bytes, os.PathLike]

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from libtbx.mpi4py import MPI
comm = MPI.COMM_WORLD


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LOGGING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class UsageLogManager:
  date_fmt = '%Y-%m-%d %H:%M:%S,uuu'
  fmt = '%(asctime)s - %(message)s'
  formatter = logging.Formatter(fmt=fmt, datefmt=date_fmt)
  line_regex = re.compile(
    r'(\d{2,4}-\d\d-\d\d \d\d:\d\d:\d\d,\d+) - '
    r'UsageStats\(cpu_usage=(-?\d+\.?\d*), ?cpu_memory=(-?\d+\.?\d*), ?'
    r'gpu_usage=(-?\d+\.?\d*), ?gpu_memory=(-?\d+\.?\d*)')

  def __init__(self, logger_name: str) -> None:
    self.log = logging.getLogger(logger_name)

  def get_file_logger(self, file_name: str) -> logging.Logger:
    self.log.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler(filename=file_name)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(self.formatter)
    self.log.addHandler(file_handler)
    return self.log

  def get_null_logger(self) -> logging.Logger:
    self.log.setLevel(logging.CRITICAL + 1)
    self.log.addHandler(logging.NullHandler())
    return self.log


usage_log_manager = UsageLogManager('cctbx.usage')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ STORING USAGE STATS ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class UnitIntervalFloat(float):
  """Convenience wrapper that always returns float of value between 0 and 1"""
  def __new__(cls, value) -> 'UnitIntervalFloat':
    return float.__new__(float, max(min(value, 1.0), 0.0))


@dataclass
class UsageStats:
  cpu_usage: UnitIntervalFloat = 0.0
  cpu_memory: UnitIntervalFloat = 0.0
  gpu_usage: UnitIntervalFloat = 0.0
  gpu_memory: UnitIntervalFloat = 0.0

  @property
  def vector(self) -> np.ndarray:
    return np.array(list(astuple(self)), dtype=float)

  def __add__(self, other: 'UsageStats') -> 'UsageStats':
    return self.__class__(*(self.vector + other.vector))  # noqa

  def __mul__(self, other: float) -> 'UsageStats':
    return self.__class__(*(self.vector * other))  # noqa

  def __truediv__(self, other: float) -> 'UsageStats':
    return self.__class__(*(self.vector / other))  # noqa


class UsageStatsHistory(UserDict):
  """A convenience class to store and manipulate"""

  @classmethod
  def from_file(cls, path: PathLike) -> 'UsageStatsHistory':
    new = cls()
    with open(path, 'r') as file:
      for line in file.readlines():
        if match := UsageLogManager.line_regex.match(line):
          t = datetime.strptime(match.group(1), UsageLogManager.date_fmt)
          us = UsageStats(*match.group(2, 3, 4, 5))
          new[t] = us
    return new

  def get_deltas_in_minutes(self) -> np.ndarray:
    times = np.array(list(self.keys()), dtype='datetime64')
    return (times - times[0]) / np.timedelta64(1, 'm')

  def get_stats(self, key: str) -> np.ndarray:
    return np.array([getattr(u, key) for u in self.values()])


# ~~~~~~~~~~~~~~~~~~~~~~~~~ COLLECTING USAGE STATS ~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class ResourceProbeException(OSError):
  pass


class NoSuitableResourceProbeException(OSError):
  pass


class BaseResourceProbeType(type):
  REGISTRY: dict

  def __new__(mcs, *args, **kwargs):
    new_cls = type.__new__(mcs, *args, **kwargs)
    if kind := getattr(mcs, 'kind', ''):
      mcs.REGISTRY[kind] = new_cls
    return new_cls

  @classmethod
  def discover(mcs) -> Union[Type['BaseCPUResourceProbe'],
                             Type['BaseGPUResourceProbe']]:
    for _, probe_class in mcs.REGISTRY.items():
      try:
        probe_instance = probe_class()
        _ = probe_instance.get_name()
        _ = probe_instance.get_usage_stats()
      except ResourceProbeException:
        continue
      else:
        return probe_class
    raise NoSuitableResourceProbeException


class CPUResourceProbeType(BaseResourceProbeType):
  REGISTRY = {}


class GPUResourceProbeType(BaseResourceProbeType):
  REGISTRY = {}


class BaseCPUResourceProbe(metaclass=CPUResourceProbeType):
  kind = None

  def get_name(self) -> str:  # noqa
    return platform.node()

  def get_usage_stats(self) -> UsageStats:  # noqa
    return NotImplemented


class BaseGPUResourceProbe(metaclass=CPUResourceProbeType):
  kind = None

  def get_name(self) -> str:  # noqa
    return NotImplemented

  def get_usage_stats(self) -> UsageStats:  # noqa
    return NotImplemented


class PsutilCPUResourceProbe(BaseCPUResourceProbe):
  kind = 'psutil'

  def __init__(self):
    try:
      import psutil
    except ImportError as e:
      raise ResourceProbeException from e
    self.process = psutil.Process()

  def get_usage_stats(self) -> UsageStats:
    cpu_usage = UnitIntervalFloat(self.process.cpu_percent(interval=None))
    cpu_memory = UnitIntervalFloat(self.process.memory_percent())
    return UsageStats(cpu_usage=cpu_usage, cpu_memory=cpu_memory)


class NvidiaGPUResourceProbe(metaclass=BaseGPUResourceProbe):
  kind = 'Nvidia'

  def get_name(self) -> str:  # noqa
    args = ['nvidia-smi', '--list-gpus']
    out = subprocess.run(args, stdout=subprocess.PIPE).stdout.decode('utf-8')
    return out.strip()

  def get_usage_stats(self) -> UsageStats:  # noqa
    args = ['nvidia-smi',
            '--query-gpu=utilization.gpu,memory.used,memory.total',
            '--format=csv,noheader,nounits']
    try:
      out = subprocess.run(args, stdout=subprocess.PIPE).stdout.decode('utf-8')
      values = [float(v) for v in out.replace(',', ' ').strip().split()]
      gpu_usage = UnitIntervalFloat(mean(values[::3]))
      gpu_memory = UnitIntervalFloat(mean(values[1::3]) / mean(values[2::3]))
      return UsageStats(gpu_usage=gpu_usage, gpu_memory=gpu_memory)
    except (FileNotFoundError, ValueError) as e:
      raise ResourceProbeException from e


class RankInfo:
  """Helper object that stores info about self and nearby ranks"""
  def __init__(self):
    self.rank = comm.rank
    self.cpu_probe = CPUResourceProbeType.discover()()
    self.gpu_probe = GPUResourceProbeType.discover()()
    self.node = self.cpu_probe.get_name()
    self.gpu = self.gpu_probe.get_name()
    neighborhood = comm.allgather((self.rank, self.node, self.gpu))
    self.same_node_ranks = [r for r, n, g in neighborhood if n == self.node]
    self.same_gpu_ranks = [r for r, n, g in neighborhood if g == self.gpu]
    self.is_gpu_ambassador = self.rank == min(self.same_gpu_ranks)
    self.is_node_ambassador = self.rank == min(self.same_node_ranks)

  def get_rank_usage_stats(self) -> UsageStats:
    cpu_usage_stats = self.cpu_probe.get_usage_stats()
    gpu_usage_stats = self.gpu_probe.get_usage_stats()
    return cpu_usage_stats + gpu_usage_stats

  def get_node_usage_stats(self) -> UsageStats:
    usage_stats = self.get_rank_usage_stats()
    node_comm_color = min(self.same_node_ranks)
    node_comm_rank = self.rank - node_comm_color
    node_comm = comm.Split(node_comm_color, node_comm_rank)
    ri_and_usage = comm.allgather((rank_info, usage_stats))
    node_comm.Free()
    cpu_usage_stats = [u for ri, u in ri_and_usage]
    gpu_usage_stats = [u for ri, u in ri_and_usage if ri.is_gpu_ambassador]
    cpu_usage_stat_sum = sum(cpu_usage_stats, UsageStats())
    gpu_usage_stat_sum = sum(gpu_usage_stats, UsageStats())
    return UsageStats(
      cpu_usage=cpu_usage_stat_sum.cpu_usage / len(cpu_usage_stats),
      cpu_memory=cpu_usage_stat_sum.cpu_memory,
      gpu_usage=gpu_usage_stat_sum.gpu_usage / len(gpu_usage_stats),
      gpu_memory=gpu_usage_stat_sum.gpu_memory / len(gpu_usage_stats),
    )


rank_info = RankInfo()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MONITORING USAGE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


#from xfel.util.usage_monitor import UsageMonitor#\
#
#@UsageMonitor(detail='rank', period=1.0)
#function_to_be_timed()
#
#with UsageMonitor(detail='rank', period=1.0):
#  function_to_be_timed()


class UsageMonitor(ContextDecorator):
  """This class can be used both as context manager (with) and @decorator"""

  class Detail(Enum):
    rank = 'rank'
    node = 'node'
    single = 'single'
    none = 'none'

  def __init__(self, detail: str = 'node', period: float = 5.0) -> None:
    self.detail = self.Detail(detail)
    self.period: float = period  # <5 de-prioritizes sub-procs & they stop...
    self.log: logging.Logger = self.get_logger()
    self.log.info(f'Collecting CPU usage with {rank_info.cpu_probe.kind=}')
    self.log.info(f'Collecting GPU usage with {rank_info.gpu_probe.kind=}')
    self.usage_stats_history = UsageStatsHistory()
    self._daemon = None

  def __enter__(self) -> None:
    self.start_logging_daemon()

  def __exit__(self, exc_type, exc_val, exc_tb):
    self.plot_usage_stats_history()

  @property
  def usage_stats(self) -> UsageStats:
    return {
      self.Detail.rank: rank_info.get_rank_usage_stats,
      self.Detail.node: rank_info.get_node_usage_stats,
      self.Detail.single: rank_info.get_rank_usage_stats,
      self.Detail.none: (lambda _: None)
    }[self.detail]()

  @property
  def daemon(self) -> threading.Thread:
    return self._daemon

  @cached_property
  def is_logging(self) -> bool:
    return {
      self.Detail.rank: True,
      self.Detail.node: rank_info.is_node_ambassador,
      self.Detail.single: rank_info.rank == 0,
      self.Detail.none: False
    }[self.detail]

  @cached_property
  def log_path(self) -> str:
    return {
      self.Detail.rank: f'usage_{rank_info.rank}.log',
      self.Detail.node: f'usage_{rank_info.node}.log',
      self.Detail.single: 'usage.log',
      self.Detail.none: None
    }[self.detail]

  def get_logger(self) -> logging.Logger:
    return usage_log_manager.get_file_logger(self.log_path) \
      if self.is_logging else usage_log_manager.get_null_logger()

  def log_current_usage(self) -> None:
    usage_stats = self.usage_stats
    if self.is_logging:
      self.usage_stats_history[datetime.now()] = usage_stats
    self.log.info(msg=f'{usage_stats}')

  def log_usage_every_period(self) -> None:
    """Warning: call as threading daemon only, otherwise will never stop"""
    start = datetime.now()
    log_iter_counter = itertools.count(start=0, step=1)
    for log_iter in log_iter_counter:
      threading.Thread(target=self.log_current_usage, args=()).start()
      next_iter_time = start + timedelta(seconds=self.period) * (log_iter + 1)
      time.sleep(max((next_iter_time - datetime.now()).total_seconds(), 0))

  def start_logging_daemon(self) -> None:
    self._daemon = threading.Thread(target=self.log_usage_every_period,
                                    args=(), daemon=True)
    self.daemon.start()

  def plot_usage_stats_history(self):
    usage_stats_histories = comm.gather(self.usage_stats_history, root=0)
    if rank_info.rank == 0:
      usage_stats_histories = [u for u in usage_stats_histories if u is not None]
      UsageArtist().plot(usage_stats_histories, 'usage.png')


class UsageArtist:
  def __init__(self) -> None:
    self.colormap = plt.get_cmap('tab10')
    self.colormap_period = 10
    self.fig = plt.figure(tight_layout=True)
    gs = GridSpec(4, 1, hspace=0, wspace=0)
    self.ax_cu = self.fig.add_subplot(gs[0, 0])  # cu = (c)pu (u)sage
    self.ax_cm = self.fig.add_subplot(gs[1, 0], sharex=self.ax_cu)
    self.ax_gu = self.fig.add_subplot(gs[2, 0], sharex=self.ax_cu)
    self.ax_gm = self.fig.add_subplot(gs[3, 0], sharex=self.ax_cu)

  def plot(self,
           usage_stats_histories: List[UsageStatsHistory],
           save_path: PathLike = None,
           ) -> None:
    axes = self.ax_cu, self.ax_cm, self.ax_gu, self.ax_gm
    stats = ['cpu_usage', 'cpu_memory', 'gpu_usage', 'gpu_memory']
    labels = ['CPU usage', 'CPU memory', 'GPU usage', 'GPU memory']
    for i, ush in enumerate(usage_stats_histories):
      minutes = ush.get_deltas_in_minutes()
      color = self.colormap(i % self.colormap_period)
      for ax, stat in zip(axes, stats):
        ax.plot(minutes, ush.get_stats(stat), color=color)
    for ax, label in zip(axes, labels):
      ax.set_ylabel(label)
    self.fig.set_size_inches(12, 10)
    if save_path:
      self.fig.savefig('usage.png')
    else:
      plt.show()


def plot_usage_from_logs(log_glob: PathLike, save_path: PathLike) -> None:
  """Manually plot if UsageMonitor terminates unexpectedly i.e. due to OOM"""
  histories = [UsageStatsHistory.from_file(Path(p)) for p in glob(log_glob)]
  UsageArtist().plot(usage_stats_histories=histories, save_path=save_path)
