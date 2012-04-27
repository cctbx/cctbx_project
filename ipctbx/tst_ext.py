def exercise_1(
      number_of_segments=5,
      number_of_integers=20,
      keep_segments=False):
  import ipctbx
  shmids = []
  for _ in range(number_of_segments):
    shmids.append(ipctbx.exercise_create_segment(number_of_integers))
  for shmid in shmids:
    matches = ipctbx.exercise_read_segment(shmid, number_of_integers)
    assert matches == number_of_integers
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
