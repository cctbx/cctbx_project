from __future__ import absolute_import, division, print_function

def get_memory_usage():
  '''Return memory used by the process in MB'''
  import resource
  import platform
  # getrusage returns kb on linux, bytes on mac
  units_per_mb = 1024
  if platform.system() == "Darwin":
    units_per_mb = 1024*1024
  return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / units_per_mb
