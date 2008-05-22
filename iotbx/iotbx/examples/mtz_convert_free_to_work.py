import iotbx.mtz
from cctbx.array_family import flex
import sys, os

def run(args, label="R-free-flags", convert_fraction=0.5, random_seed=0):
  assert len(args) == 1
  input_file_name = args[0]
  output_file_name = "less_free_"+os.path.basename(input_file_name)
  print "Reading file:", input_file_name
  mtz_obj = iotbx.mtz.object(file_name=input_file_name)
  column = mtz_obj.get_column(label=label)
  selection_valid = column.selection_valid()
  flags = column.extract_values()
  def get_and_report(what):
    free_indices = ((flags != 0) & selection_valid).iselection()
    work_indices = ((flags == 0) & selection_valid).iselection()
    if (  free_indices.size()
        + work_indices.size() != selection_valid.count(True)):
      raise RuntimeError("""\
Unexpected array of R-free flags:
  Expected: 0 for work reflections, 1 for test reflections.""")
    print what, "number of free reflections:", free_indices.size()
    print what, "number of work reflections:", work_indices.size()
    return free_indices
  free_indices = get_and_report("Input")
  mt = flex.mersenne_twister(seed=random_seed)
  permuted_indices = free_indices.select(
    mt.random_permutation(size=free_indices.size()))
  n_convert = int(permuted_indices.size() * convert_fraction + 0.5)
  print "Number of reflections converted from free to work:", n_convert
  flags.set_selected(permuted_indices[:n_convert], 0)
  get_and_report("Output")
  column.set_values(values=flags, selection_valid=selection_valid)
  print "Writing file:", output_file_name
  mtz_obj.write(file_name=output_file_name)

if (__name__ == "__main__"):
  run(sys.argv[1:])
