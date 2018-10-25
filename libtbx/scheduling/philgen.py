from __future__ import absolute_import, division, print_function

from libtbx.scheduling import SetupError

# Module exception
class ConfigurationError(SetupError):
  """
  Error in setting up the PHIL generator
  """

# Multiprocessing queue options
def mp_fifo_queue():

  from libtbx.scheduling import mp_handler
  return mp_handler.fifo_qfactory


def mp_managed_queue():

  from libtbx.scheduling import mp_handler
  return mp_handler.managed_qfactory()


# Cluster queue options
class cluster_file_queue(object):

  def __init__(self, prefix = "tmp", folder = ".", waittime = 0.1, multifile = None):

    self.prefix = prefix
    self.folder = folder
    self.waittime = waittime
    self.multifile = multifile


  def __call__(self, params):

    from libtbx.scheduling import file_queue
    qfac = file_queue.qfactory(
      prefix = self.prefix,
      folder = self.folder,
      waittime = self.waittime,
      )

    if self.multifile is not None:
      mqfac = file_queue.mqfactory(
        count = self.multifile,
        prefix = self.prefix,
        folder = self.folder,
        waittime = self.waittime,
        )

    else:
      mqfac = qfac

    return ( mqfac , qfac, mqfac )


  def __str__(self):

    return ""


class cluster_socket_queue(object):

  def __init__(self, port = 0, keylength = 16):

    self.port = port
    self.keylength = keylength


  def __call__(self, params):

    from libtbx.scheduling import socket_queue
    qfac = socket_queue.QFactory( port = self.port, keylength = self.keylength )

    return ( qfac , qfac, qfac )


  def __str__(self):

    return ""


# Pool lifecycle controller options
def pool_unlimited_lifecycle():

  from libtbx.scheduling import process_pool
  return process_pool.unlimited


class pool_runtime_limit_lifecycle(object):

  def __init__(self, value):

    self.seconds = value


  def __call__(self):

    from libtbx.scheduling import process_pool
    import functools
    return functools.partial( process_pool.runtime_limit, seconds = self.seconds )


class pool_wallclock_limit_lifecycle(object):

  def __init__(self, value):

    self.seconds = value


  def __call__(self):

    from libtbx.scheduling import process_pool
    import functools
    return functools.partial( process_pool.wallclock_limit, seconds = self.seconds )


class pool_jobcount_limit_lifecycle(object):

  def __init__(self, value):

    self.maxjobs = value


  def __call__(self):

    from libtbx.scheduling import process_pool
    import functools
    return functools.partial( process_pool.jobcount_limit, maxjobs = self.maxjobs )


# Pool autoscaling options
def pool_constant_capacity(ncpus):

  from libtbx.scheduling import process_pool
  return process_pool.constant_capacity( capacity = ncpus )


class pool_upscaling_capacity(object):

  def __init__(self, minimum = 0, buffering = 1):

    self.minimum = minimum
    self.buffering = buffering


  def __call__(self, ncpus):

    from libtbx.scheduling import process_pool
    return process_pool.upscaling_capacity(
      minimum = self.minimum,
      maximum = ncpus,
      buffering = self.buffering,
      )


class pool_breathing_capacity(object):

  def __init__(self, minimum = 0, buffering = 1, downsize_step = 2, downsize_delay = 5):

    self.minimum = minimum
    self.buffering = buffering
    self.downsize_step = downsize_step
    self.downsize_delay = downsize_delay


  def __call__(self, ncpus):

    from libtbx.scheduling import process_pool
    return process_pool.breathing_capacity(
      minimum = self.minimum,
      maximum = ncpus,
      buffering = self.buffering,
      downsize_step = self.downsize_step,
      downsize_delay = self.downsize_delay,
      )


# Classes describing aspects of the system
class constraint(object):
  """
  Properties of the call to be parallelized
  """

  def __init__(
    self,
    pickleable_target = False,
    pickleable_return_value = False,
    offload_main_thread = False,
    concurrent_python_execution = True
    ):

    self.pickleable_target = pickleable_target
    self.pickleable_return_value = pickleable_return_value
    self.offload_main_thread = offload_main_thread
    self.concurrent_python_execution = concurrent_python_execution


class ability(object):
  """
  Properties of the parallelization technology
  """

  def __init__(
    self,
    allows_startup_for_nonpickleable_target,
    allows_transfer_for_nonpickleable_object,
    allows_unlimited_capacity,
    allows_concurrent_python_execution,
    ):

    self.allows_startup_for_nonpickleable_target = allows_startup_for_nonpickleable_target
    self.allows_transfer_for_nonpickleable_object = allows_transfer_for_nonpickleable_object
    self.allows_unlimited_capacity = allows_unlimited_capacity
    self.allows_concurrent_python_execution = allows_concurrent_python_execution


class setting(object):
  """
  Data for setting up a parallel run
  """

  def __init__(
    self,
    jfactory,
    qfactory,
    inqfactory,
    outqfactory,
    lifecycle,
    ):

    self.jfactory = jfactory
    self.qfactory = qfactory
    self.inqfactory = inqfactory
    self.outqfactory = outqfactory
    self.lifecycle = lifecycle


# Technologies
class threading(object):
  """
  Phil generator for threading option
  """

  def __init__(
    self,
    capture_exception = False,
    pool_lifecycle = pool_unlimited_lifecycle,
    ):

    self.capture_exception = capture_exception
    self.pool_lifecycle = pool_lifecycle


  def abilities(self):

    return ability(
      allows_startup_for_nonpickleable_target = True,
      allows_transfer_for_nonpickleable_object = True,
      allows_unlimited_capacity = False,
      allows_concurrent_python_execution = False,
      )


  def jfactory(self):

    if self.capture_exception:
      from libtbx.scheduling import thread_handler
      jfactory = thread_handler.exception_capturing_thread

    else:
      import threading
      jfactory = threading.Thread

    return jfactory


  def qfactory(self):

    from libtbx.scheduling import thread_handler

    return ( thread_handler.qfactory, ) * 3


  def settings(self, params):

    # No user-controlled parameters expected
    ( qfac, inqfac, outqfac ) = self.qfactory()

    return setting(
      jfactory = self.jfactory(),
      qfactory = qfac,
      inqfactory = inqfac,
      outqfactory = outqfac,
      lifecycle = self.pool_lifecycle,
      )


  def phil(self):

    return ""


class multiprocessing(object):
  """
  Phil generator for multiprocessing option
  """

  def __init__(
    self,
    capture_stderr = True,
    qtype = mp_fifo_queue,
    pool_lifecycle = pool_unlimited_lifecycle,
    ):

    self.capture_stderr = capture_stderr
    self.qtype = qtype
    self.pool_lifecycle = pool_lifecycle


  def abilities(self):

    import sys
    return ability(
      allows_startup_for_nonpickleable_target = not sys.platform.startswith( "win32" ),
      allows_transfer_for_nonpickleable_object = False,
      allows_unlimited_capacity = False,
      allows_concurrent_python_execution = True,
      )


  def jfactory(self):

    if self.capture_stderr:
      from libtbx.scheduling import mp_handler
      jfactory = mp_handler.stderr_capturing_process

    else:
      import multiprocessing
      jfactory = multiprocessing.Process

    return jfactory


  def qfactory(self):

    return ( self.qtype(), ) * 3


  def settings(self, params):

    # No user-controlled parameters expected
    ( qfac, inqfac, outqfac ) = self.qfactory()

    return setting(
      jfactory = self.jfactory(),
      qfactory = qfac,
      inqfactory = inqfac,
      outqfactory = outqfac,
      lifecycle = self.pool_lifecycle,
      )


  def phil(self):

    return ""


class cluster(object):
  """
  Phil generator for managed cluster option (processing package)
  """

  def __init__(
    self,
    platforms = [ "sge", "lsf", "pbs", "pbspro", "condor", "slurm" ],
    queue_philgen_for = {
      "file": cluster_file_queue(),
      "network": cluster_socket_queue(),
      },
    default_queue_option = "file",
    lifecycle_factory_for = {
      "runtime": pool_runtime_limit_lifecycle,
      "wallclock": pool_wallclock_limit_lifecycle,
      },
    default_lifecycle_factory = pool_unlimited_lifecycle,
    name = "libtbx_python",
    asynchronous = True,
    use_target_file = True,
    include = None,
    poller = None,
    handler = None,
    capture_stderr = True,
    ):

    assert platforms
    assert queue_philgen_for
    assert default_queue_option in queue_philgen_for

    self.platforms = platforms
    self.queue_philgen_for = queue_philgen_for
    self.default_queue_option = default_queue_option
    self.lifecycle_factory_for = lifecycle_factory_for
    self.default_lifecycle_factory = default_lifecycle_factory

    # These settings are not user controlled
    self.name = name
    self.asynchronous = asynchronous
    self.use_target_file = use_target_file
    self.include = include
    self.poller = poller
    self.handler = handler
    self.capture_stderr = capture_stderr


  def abilities(self):

    return ability(
      allows_startup_for_nonpickleable_target = False,
      allows_transfer_for_nonpickleable_object = False,
      allows_unlimited_capacity = True,
      allows_concurrent_python_execution = True,
      )


  def jfactory(self, platform, command):

    from libtbx.scheduling import cluster_handler

    return cluster_handler.JobFactory(
      platform = platform,
      name = self.name,
      command = command,
      asynchronous = self.asynchronous,
      use_target_file = self.use_target_file,
      include = self.include,
      poller = self.poller,
      handler = self.handler,
      save_error = self.capture_stderr,
      )


  def settings(self, params):

    ( qfac, inqfac, outqfac ) = self.queue_philgen_for[ params.channel.use ](
      params = getattr( params.channel, params.channel.use, None ),
      )

    if params.limit.type is None:
      lifecycle = self.default_lifecycle_factory

    else:
      lifecycle = self.lifecycle_factory_for[ params.limit.type ]( value = params.limit.value )

    return setting(
      jfactory = self.jfactory( platform = params.platform, command = params.command ),
      qfactory = qfac,
      inqfactory = inqfac,
      outqfactory = outqfac,
      lifecycle = lifecycle,
      )


  def phil(self):

    platform = """
platform = %(platform)s
  .help = "Management software platform"
  .type = choice
  .optional = False
""" % {
  "platform": phil_choice( choices = self.platforms, default = self.platforms[0] ),
  }

    command = """
command = None
  .help = "Custom submission command"
  .type = strings
  .optional = False
"""

    limit = """
limit
  .help = "Job limits"
{
  type = %(choices)s
    .help = "Limit type"
    .type = choice
    .optional = True

  value = None
    .help = "Limit value (s)"
    .type = int( value_min = 1 )
    .optional = False
}
""" % {
  "choices": phil_choice( choices = self.lifecycle_factory_for, default = None ),
  }

    queues = """
channel
  .help = "Data channel setup"
{
  use = %(choices)s
    .help = "Channel type"
    .type = choice
    .optional = False

  %(options)s
}
""" % {
  "choices": phil_choice(
    choices = self.queue_philgen_for,
    default = self.default_queue_option,
    ),
  "options": "\n".join( "%s%s" % ( k ,v ) for ( k, v ) in self.queue_philgen_for.items() if str( v ) )
  }


    return """
.help = "Managed cluster parameters"
{
%(platform)s
%(command)s
%(limit)s
%(queues)s
}
""" % {
  "platform": platform,
  "command": command,
  "limit": limit,
  "queues": queues,
  }


# Scheduler backends
class main_thread_backend(object):
  """
  Creates a mainthread.manager

  Note: this is a singleton object
  """

  @staticmethod
  def is_valid_option(constraint):

    if constraint.offload_main_thread:
      return False

    return True


  @staticmethod
  def creator():

    from libtbx.scheduling import mainthread

    return mainthread.creator


class job_scheduler_backend(object):
  """
  Creates a job_scheduler.manager
  """

  def __init__(self, waittime = 0.01):

    self.waittime = waittime


  def is_valid_option(self, ability, constraint):

    if not constraint.pickleable_target and not ability.allows_startup_for_nonpickleable_target:
      return False

    if not constraint.pickleable_return_value and not ability.allows_transfer_for_nonpickleable_object:
      return False

    return True


  def creator(self, technology, ncpus, params):

    from libtbx.scheduling import job_scheduler

    if ncpus is None: # special for unlimited option
      ability = technology.abilities()

      if not ability.allows_unlimited_capacity:
        raise ConfigurationError("option does not allow unlimited capacity (ncpus=None)")

      else:
        capacity = job_scheduler.unlimited

    else:
      capacity = job_scheduler.limited( njobs = ncpus )

    setting = technology.settings( params = params )

    return job_scheduler.creator(
      job_factory = setting.jfactory,
      queue_factory = setting.qfactory,
      capacity = capacity,
      waittime = self.waittime,
      )


class process_pool_backend(object):
  """
  Creates a process_pool.manager
  """

  def __init__(
    self,
    autoscaling = pool_breathing_capacity(),
    waittime = 0.01,
    stalltime = 2,
    idle_timeout = 120,
    ):

    self.autoscaling = autoscaling
    self.waittime = waittime
    self.stalltime = stalltime
    self.idle_timeout = idle_timeout


  def is_valid_option(self, ability, constraint):

    if not constraint.pickleable_target and not ability.allows_transfer_for_nonpickleable_object:
      return False

    if not constraint.pickleable_return_value and not ability.allows_transfer_for_nonpickleable_object:
      return False

    return True


  def creator(self, technology, ncpus, params):

    from libtbx.scheduling import process_pool

    if ncpus is None: # special for unlimited option
      raise ConfigurationError("option does not allow unlimited capacity (ncpus=None)")

    setting = technology.settings( params = params )

    return process_pool.creator(
      job_factory = setting.jfactory,
      inq_factory = setting.inqfactory,
      outq_factory = setting.outqfactory,
      autoscaling = self.autoscaling( ncpus = ncpus ),
      lifecycle = setting.lifecycle(),
      waittime = self.waittime,
      stalltime = self.stalltime,
      idle_timeout = self.idle_timeout,
      )


# Parallelisation options
class parallel_option(object):
  """
  Valid technology/backend combination
  """

  def __init__(self, caption, technology, backend):

    self.caption = caption
    self.technology = technology
    self.backend = backend


  def __call__(self, ncpus, params):

    return self.backend.creator(
      technology = self.technology,
      ncpus = ncpus,
      params = params,
      )


class single_cpu_option(object):
  """
  Valid technology/backend combination forcing single CPU
  """

  def __init__(self, caption, technology, backend):

    self.caption = caption
    self.technology = technology
    self.backend = backend


  def __call__(self, ncpus, params):

    return self.backend.creator(
      technology = self.technology,
      ncpus = 1,
      params = params,
      )


class main_thread_option(object):
  """
  Special option for mainthread execution
  """

  def __init__(self, caption):

    self.caption = caption


  def __call__(self, ncpus, params):

    return main_thread_backend.creator()


class single(object):
  """
  No parallel execution
  """

  def __init__(
    self,
    caption = "single",
    backend = job_scheduler_backend(),
    ):

    self.caption = caption
    self.backend = backend


  def __call__(self, constraint, technologies):

    if main_thread_backend.is_valid_option( constraint = constraint ):
      return main_thread_option( caption = self.caption )

    else:
      for technology in technologies:
        if self.backend.is_valid_option(
          ability = technology.abilities(),
          contraint = constraint,
          ):
          return single_cpu_option(
            caption = self.caption,
            technology = technology,
            backend = self.backend,
            )

      raise ConfigurationError("No valid '%s' option with these constraints" % self.caption)


class multicore(object):
  """
  Parallel execution within a physical machine
  """

  def __init__(
    self,
    caption = "multicore",
    backends = [
      process_pool_backend(),
      job_scheduler_backend(),
      ],
    ):

    self.caption = caption
    self.backends = backends


  def __call__(self, constraint, technologies):

    if constraint.concurrent_python_execution:
      usefuls = [ t for t in technologies if t.abilities().allows_concurrent_python_execution ]

    else:
      usefuls = technologies

    import itertools

    for ( backend, tech ) in itertools.product( self.backends, usefuls ):
      if backend.is_valid_option( ability = tech.abilities(), constraint = constraint ):
        return parallel_option(
          caption = self.caption,
          technology = tech,
          backend = backend,
          )

    raise ConfigurationError("No valid '%s' option with these constraints" % self.caption)


class batch(object):
  """
  Parallel execution via a batch queue
  """

  def __init__(
    self,
    caption = "cluster",
    backend = process_pool_backend(),
    conf_cluster = cluster(),
    ):

    self.caption = caption
    self.backend = backend
    self.conf_cluster = conf_cluster


  def __call__(self, constraint):

    if self.backend.is_valid_option(
      ability = self.conf_cluster.abilities(),
      constraint = constraint,
      ):
      return parallel_option(
        caption = self.caption,
        technology = self.conf_cluster,
        backend = self.backend,
        )

    raise ConfigurationError("No valid '%s' option with these constraints" % self.caption)


  def __str__(self):

    return "%s%s" % ( self.caption, self.conf_cluster.phil() )


class manager(object):
  """
  Creates PHIL scope for setting up a manager
  """

  def __init__(
    self,
    constraint,
    caption = "concurrency",
    single = single(),
    multicore = multicore(),
    extras = [ batch() ],
    conf_threading = threading(),
    conf_multiprocessing = multiprocessing(),
    prefer_mp = False,
    ):

    self.caption = caption
    self.single = single
    self.multicore = multicore

    if prefer_mp:
      technologies = [ conf_multiprocessing, conf_threading ]

    else:
      technologies = [ conf_threading, conf_multiprocessing ]

    self.option_for = {}

    try:
      self.option_for[ single.caption ] = single(
        constraint = constraint,
        technologies = technologies,
        )

    except ConfigurationError:
      raise ConfigurationError("'%s' option not valid with this setup" % single.caption)

    try:
      self.option_for[ multicore.caption ] = multicore(
        constraint = constraint,
        technologies = technologies,
        )

    except ConfigurationError:
      pass

    self.extras = []

    for platform in extras:
      try:
        self.option_for[ platform.caption ] = platform( constraint = constraint )

      except ConfigurationError:
        continue

      self.extras.append( platform )


  def __call__(self, params):

    scope = getattr( params, self.caption )

    if scope.ncpus == 1 and scope.technology == self.multicore.caption:
      technology = self.single.caption

    else:
      technology = scope.technology

    try:
      return self.option_for[ technology ](
        ncpus = scope.ncpus,
        params = getattr( scope, technology, None ),
        )

    except ConfigurationError as e:
      from libtbx.utils import Sorry
      raise Sorry("Error setting up '%s': %s" % ( technology, e ))


  def __str__(self):

    return """
%(caption)s
  .help = "Concurrency parameters"
{
  technology = %(technology)s
    .help = "Parallelisation option to use"
    .type = choice
    .optional = False

  ncpus = 1
    .help = "Maximum number of CPUs"
    .type = int( value_min = 1 )
    .optional = False

  %(others)s
}
""" % {
  "caption": self.caption,
  "technology": phil_choice(
    choices = self.option_for,
    default = self.single.caption,
    ),
  "others": "\n".join( str( p ) for p in self.extras ),
  }


def phil_choice(choices, default):

  if not choices:
    return None

  if default is not None and default not in choices:
    raise ValueError("Default not among choices")

  return " ".join( "*%s" % c if c == default else c for c in choices )
