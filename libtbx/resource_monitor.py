from __future__ import division

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
from typing import Iterable, Type, Union

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

from libtbx.mpi4py import MPI


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ TYPING AND UTILITY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


comm = MPI.COMM_WORLD
PathLike = Union[str, bytes, os.PathLike]


def _mkdir_if_missing(path: PathLike) -> None:
  path: Path = Path(path)
  if comm.rank == 0:
    if dir_name := Path(path).parent:
      dir_name.mkdir(parents=True, exist_ok=True)
  comm.barrier()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LOGGING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class ResourceLogManager:
  """Create appropriate resource logger for each rank and control its format"""
  date_fmt = '%Y-%m-%d %H:%M:%S,%f'  # exposed logging default, do not change
  fmt = '%(asctime)s - %(message)s'
  formatter = logging.Formatter(fmt=fmt)
  line_regex = re.compile(
    r'(\d{2,4}-\d\d-\d\d \d\d:\d\d:\d\d,\d+) - '
    r'ResourceStats\(cpu_usage=(-?\d+\.?\d*), ?cpu_memory=(-?\d+\.?\d*), ?'
    r'gpu_usage=(-?\d+\.?\d*), ?gpu_memory=(-?\d+\.?\d*)')

  def __init__(self, logger_name: str) -> None:
    self.log = logging.getLogger(logger_name)

  def get_file_logger(self, file_path: PathLike) -> logging.Logger:
    self.log.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler(filename=file_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(self.formatter)
    self.log.addHandler(file_handler)
    return self.log

  def get_null_logger(self) -> logging.Logger:
    self.log.setLevel(logging.CRITICAL + 1)
    self.log.addHandler(logging.NullHandler())
    return self.log


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STORING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class PerCentFloat(float):
  """Convenience wrapper that clamps float between 0 and 100"""
  def __new__(cls, value) -> 'PerCentFloat':
    return float.__new__(float, max(min(float(value), 100.0), 0.0))


@dataclass
class ResourceStats:
  """Convenient dataclass used to pass info about CPU/GPU usage/memory"""
  cpu_usage: PerCentFloat = 0.0
  cpu_memory: PerCentFloat = 0.0
  gpu_usage: PerCentFloat = 0.0
  gpu_memory: PerCentFloat = 0.0

  def __post_init__(self):
    for attr in ['cpu_usage', 'cpu_memory', 'gpu_usage', 'gpu_memory']:
      setattr(self, attr, PerCentFloat(getattr(self, attr)))

  @property
  def vector(self) -> np.ndarray:
    return np.array(list(astuple(self)), dtype=float)

  def __add__(self, other: 'ResourceStats') -> 'ResourceStats':
    return self.__class__(*(self.vector + other.vector))  # noqa

  def __mul__(self, other: float) -> 'ResourceStats':
    return self.__class__(*(self.vector * other))  # noqa

  def __truediv__(self, other: float) -> 'ResourceStats':
    return self.__class__(*(self.vector / other))  # noqa


try:
  _UD = UserDict[datetime, ResourceStats]
except TypeError:
  _UD = UserDict # Needed for PyVer < 3.9
class ResourceStatsHistory(_UD):
  """Store and easily handle a dictionary of datetime-ResourceStats pairs"""

  @classmethod
  def from_file(cls, path: PathLike) -> 'ResourceStatsHistory':
    new = cls()
    with open(path, 'r') as file:
      for line in file.readlines():
        if match := ResourceLogManager.line_regex.match(line):
          t = datetime.strptime(match.group(1), ResourceLogManager.date_fmt)
          new[t] = ResourceStats(*match.group(2, 3, 4, 5))
    return new

  def get_deltas_array(self) -> np.ndarray:  # result in minutes
    times = np.array(list(self.keys()), dtype=np.datetime64)
    return (times - times[0]) / np.timedelta64(1, 'm') \
        if times.size else np.empty((0,), dtype=np.float64)

  def get_stats_array(self, key: str) -> np.ndarray:
    return np.array([getattr(u, key) for u in self.values()])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PROBING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class ResourceProbeException(OSError):
  """Raised if a `ResourceProbe` fails during `get_resource_stats()`"""


class NoSuitableResourceProbeException(OSError):
  """Raised if all of available `ResourceProbe`s fail `get_resource_stats()`"""


class BaseResourceProbeType(type):
  """Metaclass that implements automatic probe `REGISTRY` and `discover`y"""
  REGISTRY: dict

  def __new__(mcs, *args, **kwargs) -> Union[Type['BaseCPUResourceProbe'],
                                             Type['BaseGPUResourceProbe']]:
    new_cls = type.__new__(mcs, *args, **kwargs)
    if kind := getattr(new_cls, 'kind', ''):
      mcs.REGISTRY[kind] = new_cls
    return new_cls

  @classmethod
  def discover(mcs) -> Union[Type['BaseCPUResourceProbe'],
                             Type['BaseGPUResourceProbe']]:
    for _, probe_class in mcs.REGISTRY.items():
      try:
        probe_instance = probe_class()
        _ = probe_instance.get_resource_stats()
      except ResourceProbeException:
        continue
      else:
        return probe_class
    raise NoSuitableResourceProbeException


class DummyResourceProbe:
  kind = 'dummy'

  def get_name(self) -> str:  # noqa
    return 'unknown'

  def get_resource_stats(self) -> ResourceStats:  # noqa
    return ResourceStats(-1., -1., -1., -1.)


class CPUResourceProbeType(BaseResourceProbeType):
  """Metaclass for CPU resource probes"""
  REGISTRY = {}


class BaseCPUResourceProbe(metaclass=CPUResourceProbeType):
  """
  Base class for CPU resource probes to be inherited by all CPU probes.
  Every subclass should define `kind`, `get_name()`, `get_resource_stats()`.
  """
  kind = None

  def get_name(self) -> str:  # noqa
    return platform.node()


class PsutilCPUResourceProbe(BaseCPUResourceProbe):
  """CPU resource probe that attempts collecting CPU usage via psutil"""
  kind = 'psutil'

  def __init__(self) -> None:
    try:
      import psutil
    except ImportError as e:
      raise ResourceProbeException from e
    self.process = psutil.Process()

  def get_resource_stats(self) -> ResourceStats:
    cpu_usage = PerCentFloat(self.process.cpu_percent(interval=None))
    cpu_memory = PerCentFloat(self.process.memory_percent())
    return ResourceStats(cpu_usage=cpu_usage, cpu_memory=cpu_memory)


class DummyCPUResourceProbe(BaseCPUResourceProbe, DummyResourceProbe):
  """CPU resource probe that reports dummy data when no good probe is found"""
  kind = 'dummy'


class GPUResourceProbeType(BaseResourceProbeType):
  """Metaclass for GPU resource probes"""
  REGISTRY = {}


class BaseGPUResourceProbe(metaclass=GPUResourceProbeType):
  """
  Base class for GPU resource probes to be inherited by all GPU probes
  Every subclass should define `kind`, `get_name()`, `get_resource_stats()`.
  """
  kind = None


class NvidiaGPUResourceProbe(BaseGPUResourceProbe):
  """GPU resource probe that attempts collecting GPU resource via nvidia-smi"""
  kind = 'Nvidia'

  def get_name(self) -> str:  # noqa
    args = ['nvidia-smi', '--list-gpus']
    out = subprocess.run(args, stdout=subprocess.PIPE).stdout.decode('utf-8')
    return out.strip()

  def get_resource_stats(self) -> ResourceStats:  # noqa
    args = ['nvidia-smi',
            '--query-gpu=utilization.gpu,memory.used,memory.total',
            '--format=csv,noheader,nounits']
    try:
      out = subprocess.run(args, stdout=subprocess.PIPE).stdout.decode('utf-8')
      values = [float(v) for v in out.replace(',', ' ').strip().split()]
      gpu_usage = PerCentFloat(mean(values[::3]))
      gpu_memory = PerCentFloat(100. * mean(values[1::3]) / mean(values[2::3]))
      return ResourceStats(gpu_usage=gpu_usage, gpu_memory=gpu_memory)
    except (FileNotFoundError, ValueError) as e:
      raise ResourceProbeException from e


class DummyGPUResourceProbe(BaseGPUResourceProbe, DummyResourceProbe):
  """CPU resource probe that reports dummy data when no good probe is found"""
  kind = 'dummy'


class RankInfo:
  """Manage info about this and neighbor ranks using auto-discovered probes"""
  def __init__(self) -> None:
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

  def get_rank_resource_stats(self) -> ResourceStats:
    cpu_resource_stats = self.cpu_probe.get_resource_stats()
    gpu_stats = self.gpu_probe.get_resource_stats()
    return cpu_resource_stats + gpu_stats


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MONITORING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class ResourceMonitor(ContextDecorator):
  """
  Collect and log information about CPU & GPU resources in decorated processes.
  In each case the information is limited to single rank to avoid mpi comm.
  This class can be used in several ways listed below:

  - As a decorator – to monitor every call of a given function, decorate the
    function definition itself. The monitor will start every time the function
    is called and stop every time the function is terminated:
    ```
    @ResourceMonitor(*args)
    def function_to_be_monitored()
      stuff_to_be_timed()
    ```

  - As a context manager – to monitor a part of the code, place it withing
    a context manager using a `with` statement. The monitor will start at the
    start of the `with` block and terminate at the end of the `with` block:
    ```
    with ResourceMonitor(*args):
      stuff_to_be_timed()
    ```

  - As a standalone instance – for more advanced or customizable application,
    the context manager can be manually instantiated, started, stopped, etc.
    The following would be roughly equivalent to context manager approach:
    ```
    rm = ResourceMonitor(*args):
    try:
      um.start()
      stuff_to_be_timed()
    finally:
      self.stop()
    ```
  """

  class Detail(Enum):
    node = 'node'
    rank = 'rank'
    rank0 = 'rank0'
    none = 'none'

  def __init__(self,
               detail: str = 'rank',
               period: float = 5.0,
               plot: bool = True,
               prefix: str = 'monitor',
               write: bool = True,
               ) -> None:
    self.active: bool = False
    self.daemon: threading.Thread = None
    self.detail: 'ResourceMonitor.Detail' = self.Detail(detail)
    self.period: float = period  # <5 sec. de-prioritizes sub-procs & they stop
    self.plot: bool = plot
    self.prefix: PathLike = prefix if prefix else 'monitor'
    self.write: bool = write
    self.rank_info: RankInfo = RankInfo()
    _mkdir_if_missing(prefix)
    logger_name = 'libtbx.resource_monitor.' + Path(self.prefix).name
    self.log_manager = ResourceLogManager(logger_name)
    self.log: logging.Logger = self.get_logger()
    self.log.info(f'Collecting CPU stats with {self.rank_info.cpu_probe.kind=}')
    self.log.info(f'Collecting GPU stats with {self.rank_info.gpu_probe.kind=}')
    self.resource_stats_history = ResourceStatsHistory()

  def __enter__(self) -> None:
    self.start()

  def __exit__(self, exc_type, exc_val, exc_tb) -> None:
    self.stop()

  @property
  def resource_stats(self) -> ResourceStats:
    return {
      self.Detail.node: self.rank_info.get_rank_resource_stats,
      self.Detail.rank: self.rank_info.get_rank_resource_stats,
      self.Detail.rank0: self.rank_info.get_rank_resource_stats,
      self.Detail.none: (lambda _: None)
    }[self.detail]()

  @cached_property
  def is_logging(self) -> bool:
    return {
      self.Detail.node: self.rank_info.is_node_ambassador,
      self.Detail.rank: True,
      self.Detail.rank0: self.rank_info.rank == 0,
      self.Detail.none: False
    }[self.detail]

  @cached_property
  def log_path(self) -> str:
    return {
      self.Detail.node: f'{self.prefix}_{self.rank_info.node}.log',
      self.Detail.rank: f'{self.prefix}_{self.rank_info.rank}.log',
      self.Detail.rank0: f'{self.prefix}.log',
      self.Detail.none: None
    }[self.detail]

  def get_logger(self) -> logging.Logger:
    return self.log_manager.get_file_logger(self.log_path) \
      if self.is_logging and self.write else self.log_manager.get_null_logger()

  def log_current_resource_stats(self) -> None:
    """Get current usage stats and log + save them"""
    resource_stats = self.resource_stats
    if self.is_logging:
      self.resource_stats_history[datetime.now()] = resource_stats
    self.log.info(msg=f'{resource_stats}')

  def log_resource_stats_every_period(self) -> None:
    """Call in thread only, can be stopped only upon `self.active` = False"""
    start = datetime.now()
    log_iter_counter = itertools.count(start=0, step=1)
    for log_iter in log_iter_counter:
      self.snapshot()
      next_iter_time = start + timedelta(seconds=self.period) * (log_iter + 1)
      time.sleep(max((next_iter_time - datetime.now()).total_seconds(), 0))
      if not self.active:
        break

  def start(self) -> None:
    """Call to start probing resources"""
    self.active = True
    self.daemon = threading.Thread(target=self.log_resource_stats_every_period,
                                   args=(), daemon=True)
    self.daemon.start()

  def snapshot(self) -> None:
    """Call to probe resources in current single moment only"""
    threading.Thread(target=self.log_current_resource_stats, args=()).start()

  def stop(self) -> None:
    """Call to stop probing resources"""
    self.active = False
    if self.plot:
      self.plot_resource_stats_history()

  def plot_resource_stats_history(self) -> None:
    resource_stats_histories = comm.gather(self.resource_stats_history, root=0)
    if self.rank_info.rank == 0:
      resource_stats_histories = [h for h in resource_stats_histories if h]
      rsa = ResourceStatsArtist()
      rsa.plot(resource_stats_histories, f'{self.prefix}.png')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PLOTTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class ResourceStatsArtist:
  def __init__(self) -> None:
    self.colormap = plt.get_cmap('tab10')
    self.colormap_period = 10
    self.fig = plt.figure(tight_layout=True, figsize=(12, 10))
    gs = GridSpec(4, 1, figure=self.fig, hspace=0, wspace=0)
    self.ax_cu = self.fig.add_subplot(gs[0, 0])  # cu = (c)pu (u)sage
    self.ax_cm = self.fig.add_subplot(gs[1, 0], sharex=self.ax_cu)
    self.ax_gu = self.fig.add_subplot(gs[2, 0], sharex=self.ax_cu)
    self.ax_gm = self.fig.add_subplot(gs[3, 0], sharex=self.ax_cu)

  def plot(self,
           resource_stats_histories: Iterable[ResourceStatsHistory],
           save_path: PathLike = None,
           ) -> None:
    axes = self.ax_cu, self.ax_cm, self.ax_gu, self.ax_gm
    stats = ['cpu_usage', 'cpu_memory', 'gpu_usage', 'gpu_memory']
    labels = ['CPU usage', 'CPU memory', 'GPU usage', 'GPU memory']
    for i, rsh in enumerate(resource_stats_histories):
      minutes = rsh.get_deltas_array()
      color = self.colormap(i % self.colormap_period)
      for ax, stat in zip(axes, stats):
        ax.plot(minutes, rsh.get_stats_array(stat), color=color)
    self.ax_gm.set_xlabel('Time since first probe [min]', fontsize=12)
    for ax, label in zip(axes, labels):
      ax.grid(True)
      ax.tick_params(axis='both', which='major', labelsize=12)
      ax.set_ylabel(label + ' [%]', fontsize=12)
    if save_path:
      self.fig.savefig(f'{save_path}')
    else:
      plt.show()


def plot_logs(log_glob: PathLike, save_path: PathLike) -> None:
  """Manually plot if ResourceMonitor terminates unexpectedly i.e.due to OOM"""
  histories = [ResourceStatsHistory.from_file(Path(p)) for p in glob(log_glob)]
  rsa = ResourceStatsArtist()
  rsa.plot(resource_stats_histories=histories, save_path=save_path)
