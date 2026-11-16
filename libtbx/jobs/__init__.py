"""Job lifecycle management with pluggable execution backends.

``JobManager`` owns a table of jobs, a FIFO, and a reaper thread, and
drives every state transition. A ``Backend`` runs the jobs: locally as
subprocesses, or on a batch scheduler such as SLURM.
"""

from libtbx.jobs.backends.base import Backend, PollResult, SubmitError  # noqa: F401
from libtbx.jobs.manager import (  # noqa: F401
  Job, JobManager, JobSpec, TERMINAL_STATES)
