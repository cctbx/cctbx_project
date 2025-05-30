"""
Extract R-free flags from reflection files.
Writes reflection count, free fraction pairs to file (for plotting).
Also shows free fraction in bins.
"""
from __future__ import absolute_import, division, print_function

from iotbx import reflection_file_utils
from iotbx import reflection_file_reader
import libtbx.phil
from libtbx.str_utils import show_string
from libtbx.utils import Sorry, Usage
import sys, os
from six.moves import zip

master_params = libtbx.phil.parse("""\
r_free_flags_accumulation {
  file_name=None
    .type=path
  label=None
    .type=str
  test_flag_value=None
    .type=int
  disable_suitability_test=False
    .type=bool
}
r_free_flags_accumulation {
  output=None
    .type=path
  plot_header=None
    .type=str
    .multiple=True
}
""")

def run(args, command_name="iotbx.r_free_flags_accumulation"):
  def raise_usage():
    raise Usage("%s reflection_file [label=value]" % command_name)
  if (len(args) == 0 or "--help" in args or "-h" in args):
    raise_usage()
  phil_objects = []
  argument_interpreter = master_params.command_line_argument_interpreter(
    home_scope="r_free_flags_accumulation")
  reflection_files = []
  for arg in args:
    if (os.path.isfile(arg)):
      refl_file = reflection_file_reader.any_reflection_file(
        file_name=arg)
      if (refl_file.file_type() is not None):
        reflection_files.append(refl_file)
        arg = None
    if (arg is not None):
      try: command_line_params = argument_interpreter.process(arg=arg)
      except KeyboardInterrupt: raise
      except Exception: raise Sorry("Unknown file or keyword: %s" % arg)
      else: phil_objects.append(command_line_params)
  params_scope = master_params.fetch(sources=phil_objects).extract()
  params = params_scope.r_free_flags_accumulation
  srv = reflection_file_utils.reflection_file_server(
    reflection_files=reflection_files)
  r_free_flags, test_flag_value = srv.get_r_free_flags(
    file_name=params.file_name,
    label=params.label,
    test_flag_value=params.test_flag_value,
    disable_suitability_test=params.disable_suitability_test,
    parameter_scope="r_free_flags_accumulation")
  params.file_name = r_free_flags.info().source
  params.label = r_free_flags.info().label_string()
  params.test_flag_value = test_flag_value
  if (params.output is None):
    params.output = os.path.basename(params.file_name) \
                  + ".r_free_flags_accumulation"
  working_params = master_params.format(python_object=params_scope)
  working_params.show()
  print()
  print("#phil __OFF__")
  r_free_flags = r_free_flags.array(
    data=r_free_flags.data()==params.test_flag_value)
  r_free_flags.show_r_free_flags_info()
  print()
  accu = r_free_flags \
    .sort(by_value="resolution") \
    .r_free_flags_accumulation()
  print("Writing file: %s" % show_string(params.output))
  print("  1. column: reflection counts, sorted by resolution")
  print("  2. column: number of free reflections / total number of reflections")
  sys.stdout.flush()
  out = open(params.output, "w")
  for line in params.plot_header:
    print(line, file=out)
  for c,f in zip(accu.reflection_counts, accu.free_fractions):
    print(c, f, file=out)
  out.close()
  print()
  sys.stdout.flush()

if (__name__ == "__main__"):
  run(sys.argv[1:])
