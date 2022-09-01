## Custom merging workers

`cctbx.xfel.merge` supports custom workers located or linked in the location
defined by environment variable `XFEL_CUSTOM_WORKER_PATH`. If this feature is
used, the custom worker directory should be structured analogously to
`cctbx_project/xfel/merging/application`. This example makes a worker called
`annulus` available:
```
$ export XFEL_CUSTOM_WORKER_PATH=~/psii_spread/merging/application
$ tree $XFEL_CUSTOM_WORKER_PATH
/net/cci/dwpaley/psii_spread/merging/application
└── annulus
    ├── annulus_statistics.py
    ├── factory.py
    ├── phil.py
    └── spread_roi.py
```

### Requirements for all merging workers

- Inherit from `xfel.merging.application.worker.worker`

- Implement `__init__` and `__repr__` methods that look like this:
```
class polarization(worker):
  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(polarization, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)
  def __repr__(self):
    return 'Apply polarization correction'
```

- Implement a `run` method that takes `experiments, reflections` and returns
  the same

- Also include a `factory`, see `xfel/merging/application/modify/factory.py`
  for an example.

### Special requirement for custom workers

- The factory typically imports the worker in order to construct and return it.
  A trick is required in the custom factory to temporarily modify the import
  path:
```
# Required boilerplate for custom workers
import sys, inspect, os
current_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.insert(0, current_dir)
# End boilerplate. You may import directly from the directory containing this
# module. It's not necessary to reset the path.
from unit_cell_statistics import unit_cell_statistics
from beam_statistics import beam_statistics
```

### Optional custom phil strings

- Any file `phil.py` under `$XFEL_CUSTOM_WORKER_PATH` will be checked
  for a module-level variable `phil_str` which will be added to the global
  `cctbx.xfel.merge` phil string before parsing.

