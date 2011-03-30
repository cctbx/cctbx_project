from libtbx import easy_pickle
from libtbx.option_parser import option_parser
from scitbx.python_utils import dicts
from cctbx import euclidean_model_matching as emma
import sys

class match_record(object):

  def __init__(self, n_matches, model_size):
    self.n_matches = n_matches
    self.model_size = model_size

  def __repr__(self):
    return "%d/%d" % (self.n_matches, self.model_size)

def run():
  command_line = (option_parser(
    usage="usage: cctbx.euclidean_model_matching [OPTIONS] "
          "reference_structure.pickle structure.pickle",
    description="")
    .option("--tolerance",
            type="float",
            default=3)
    .option("--match_hydrogens", type='bool', default=True)
  ).process(args=sys.argv[1:])
  if len(command_line.args) != 2:
    command_line.parser.print_help()
    sys.exit(1)
  reference_structure = easy_pickle.load(command_line.args[0])
  if (type(reference_structure) in (type([]), type(()))):
    reference_structure = reference_structure[0]
  structures = easy_pickle.load(command_line.args[1])
  if (not type(structures) in (type([]), type(()))):
    structures = [structures]

  if not command_line.options.match_hydrogens:
    reference_structure.select_inplace(
      ~reference_structure.element_selection('H'))
    for structure in structures:
      structure.select_inplace(~structure.element_selection('H'))
  print "Reference model:"
  reference_structure.show_summary()
  print
  reference_model = reference_structure.as_emma_model()

  match_list = []
  match_histogram = dicts.with_default_value(0)
  for structure in structures:
    structure.show_summary()
    if (hasattr(structure, "info")):
      print structure.info
    print
    sys.stdout.flush()
    refined_matches = emma.model_matches(
      reference_model,
      structure.as_emma_model(),
      tolerance=command_line.options.tolerance,
      models_are_diffraction_index_equivalent=False,
      break_if_match_with_no_singles=True).refined_matches
    if (len(refined_matches)):
      refined_matches[0].show()
      m = len(refined_matches[0].pairs)
    else:
      print "No matches"
      m = 0
    match_list.append(match_record(m, structure.scatterers().size()))
    match_histogram[m] += 1
    print
    sys.stdout.flush()
  print "match_list:", match_list
  keys = match_histogram.keys()
  keys.sort()
  keys.reverse()
  print "matches: frequency"
  sum = 0
  for key in keys:
    v = match_histogram[key]
    sum += v
  s = 0
  for key in keys:
    v = match_histogram[key]
    s += v
    print "  %3d: %3d = %5.1f%%, %5.1f%%" % (key, v, 100.*v/sum, 100.*s/sum)
  print
  sys.stdout.flush()

if (__name__ == "__main__"):
  run()
