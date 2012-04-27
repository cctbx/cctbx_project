import boost.python
ext = boost.python.import_ext("ipctbx_ext")
from ipctbx_ext import *

def delete_all_segments():
  shmids = list_segment_ids(attached_status=-1, ignore_errors=False)
  result = 0
  for shmid in shmids:
    if (shmctl_rmid(shmid=shmid, ignore_errors=True)):
      result += 1
  return result
