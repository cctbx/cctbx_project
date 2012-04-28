def exercise_1(repetition_factor=5, keep_segments=False):
  import ipctbx
  from scitbx.array_family import flex
  from libtbx.math_utils import iround
  shmids = []
  for i in range(repetition_factor):
    data = flex.int_range(i,2*i+3)
    shmids.append(ipctbx.copy_to_new_segment_int(data=data))
    data = flex.int_range(i,2*i+4).as_double()
    shmids.append(ipctbx.copy_to_new_segment_double(data=data))
  for i,shmid in enumerate(shmids):
    if (i % 2 == 0):
      data = ipctbx.copy_from_segment_int(shmid)
      assert list(data) == range(i//2,i+3)
    else:
      data = ipctbx.copy_from_segment_double(shmid)
      assert [iround(_) for _ in list(data)] == range(i//2,2*(i//2)+4)
  shmids_kernel = ipctbx.list_segment_ids(
    attached_status=1, ignore_errors=False)
  assert set(shmids_kernel).issuperset(set(shmids))
  if (keep_segments):
    print "Keeping segments alive."
  else:
    for shmid in shmids:
      assert ipctbx.shmctl_rmid(shmid=shmid, ignore_errors=False)
  shmids_kernel = ipctbx.list_segment_ids(
    attached_status=-1, ignore_errors=False)
  assert len(set(shmids_kernel).intersection(set(shmids))) == 0

def run(args):
  assert args in [[], ["--keep-segments"]]
  exercise_1(keep_segments=("--keep-segments" in args))
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
