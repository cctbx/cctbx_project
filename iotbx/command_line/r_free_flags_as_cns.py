from iotbx import reflection_file_reader
import iotbx.cns.miller_array
import sys, os

def run(args):
  assert len(args) == 1
  reflection_file = reflection_file_reader.any_reflection_file(
    file_name=args[0])
  assert reflection_file.file_type is not None
  r_free_flags = None
  for miller_array in reflection_file.as_miller_arrays():
    if (miller_array.is_integer_array()):
      assert r_free_flags is None
      r_free_flags = miller_array
  assert r_free_flags is not None
  file_name = os.path.splitext(os.path.basename(args[0]))[0]+".cns"
  assert not os.path.isfile(file_name)
  print "Writing:", file_name
  r_free_flags.export_as_cns_hkl(
    file_object=open(file_name, "w"),
    array_names=["TEST"])

if (__name__ == "__main__"):
  run(sys.argv[1:])
