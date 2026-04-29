"""Convert Fortran sources to C++"""
from __future__ import absolute_import, division, print_function
import fable.cout

import hashlib
import optparse
import os
import sys

def compute_hexdigest(text):
  m = hashlib.md5()
  m.update(text.encode("utf-8"))
  return m.hexdigest()

def check_fingerprint(file_name):
  with open(file_name) as f:
    lines = f.read().splitlines()
  if (len(lines) == 0): return None
  flds = lines[0].split()
  if (len(flds) < 2 or flds[-2] != "fingerprint"): return None
  orig_hexdigest = flds[-1]
  curr_text = "\n".join(lines[1:])+"\n"
  curr_hexdigest = compute_hexdigest(text=curr_text)
  if (len(orig_hexdigest) != len(curr_hexdigest)): return None
  return (orig_hexdigest == curr_hexdigest)

def write_only_if_safe(file_name, text):
  from libtbx.str_utils import show_string
  if (os.path.exists(file_name)):
    if (not os.path.isfile(file_name)):
      raise RuntimeError(
        "Not a regular file: %s" % show_string(file_name))
    stat = check_fingerprint(file_name=file_name)
    if (stat is None or not stat):
      raise RuntimeError(
        "File appears to be manually modified: %s" % show_string(file_name))
  hexdigest = compute_hexdigest(text=text)
  with open(file_name, "w") as f:
    f.write("// fingerprint %s\n" % hexdigest)
    f.write(text)

class process(object):

  __slots__ = ["options", "dynamic_parameters", "n_calls"]

  def __init__(O, options):
    O.options = options
    if (options.dynamic_parameter is None):
      O.dynamic_parameters = None
    else:
      from fable.cout import dynamic_parameter_props
      from libtbx.utils import Sorry
      O.dynamic_parameters = []
      for opt_dp in options.dynamic_parameter:
        flds = opt_dp.replace("=", " ").split()
        if (len(flds) != 3):
          raise Sorry('Invalid --dynamic-parameter="%s"' % opt_dp)
        if (flds[1] in O.dynamic_parameters):
          raise Sorry('Duplicate --dynamic-parameter="%s"' % opt_dp)
        O.dynamic_parameters.append(dynamic_parameter_props(
          name=flds[1],
          ctype=flds[0],
          default=flds[2]))
    O.n_calls = 0

  def __call__(O, file_names):
    if (O.n_calls != 0):
      print()
    O.n_calls += 1
    opts = O.options
    lines = fable.cout.process(
      file_names=file_names,
      top_procedures=opts.top_procedure,
      include_guard_suffix=opts.include_guard_suffix,
      dynamic_parameters=O.dynamic_parameters,
      fortran_file_comments=opts.fortran_file_comments,
      fem_do_safe=not opts.no_fem_do_safe,
      arr_nd_size_max=opts.arr_nd_size_max,
      inline_all=opts.inline_all,
      common_equivalence_simple=set(opts.common_equivalence_simple.split(",")),
      namespace=opts.namespace,
      separate_cmn_hpp=opts.separate_cmn_hpp,
      number_of_function_files=opts.number_of_function_files,
      debug=opts.debug)
    text = "\n".join(lines)+"\n"
    if (opts.top_procedure is None or not opts.debug):
      sys.stdout.write(text)
    if (len(file_names) != 0 and opts.compile):
      print()
      write_only_if_safe(file_name="fable_cout.cpp", text=text)
      from fable import simple_compilation
      comp_env = simple_compilation.environment()
      out_name = comp_env.build(exe_name=opts.exe_name,
        link=opts.link, file_name_cpp="fable_cout.cpp", show_command=True)
      print()
      if (opts.run):
        from libtbx import easy_run
        cmd = os.path.join(".", out_name)
        if (opts.valgrind):
          cmd = "valgrind " + cmd
        print(cmd)
        easy_run.call(command=cmd)

def run(args):
  import libtbx.load_env
  if (len(args) == 0):
    args = ["--help"]
  elif (args == ["--example"]):
    args = [
      libtbx.env.under_dist(module_name="fable", path="test/valid/sf.f"),
      "--namespace", "example",
      "--run"]
  parser = optparse.OptionParser(usage="%s [options] fortran_file ..." % libtbx.env.dispatcher_name)
  parser.add_option("-?", action="help", help=optparse.SUPPRESS_HELP)
  parser.add_option("--compile", action="store_true", default=False)
  parser.add_option("--link", action="store_true", default=False)
  parser.add_option("--run", action="store_true", default=False)
  parser.add_option("--valgrind", action="store_true", default=False)
  parser.add_option("--each", action="store_true", default=False)
  parser.add_option("--top_procedure", action="append", type="str", metavar="IDENTIFIER")
  parser.add_option("--top-procedure", action="append", type="str", help=optparse.SUPPRESS_HELP)
  parser.add_option("--include_guard_suffix", action="store", type="str", metavar="STRING")
  parser.add_option("--include-guard-suffix", action="store", type="str", help=optparse.SUPPRESS_HELP)
  parser.add_option("--dynamic_parameter", action="append", type="str", metavar="STRING", help='example: --dynamic_parameter="int array_size=100"')
  parser.add_option("--dynamic-parameter", action="append", type="str", help=optparse.SUPPRESS_HELP)
  parser.add_option("--fortran_file_comments", action="store_true", default=False)
  parser.add_option("--fortran-file-comments", action="store_true", help=optparse.SUPPRESS_HELP)
  parser.add_option("--no_fem_do_safe", action="store_true", default=False)
  parser.add_option("--no-fem-do-safe", action="store_true", help=optparse.SUPPRESS_HELP)
  parser.add_option("--arr_nd_size_max", action="store", type="int", default=fable.cout.default_arr_nd_size_max, metavar='INTEGER (default: %d)' % fable.cout.default_arr_nd_size_max)
  parser.add_option("--arr-nd-size-max", action="store", type="int", help=optparse.SUPPRESS_HELP)
  parser.add_option("--inline_all", action="store_true", default=False)
  parser.add_option("--inline-all", action="store_true", help=optparse.SUPPRESS_HELP)
  parser.add_option("--common_equivalence_simple", action="store", type="str", default="", metavar="STRING", help='comma-separated list of common names')
  parser.add_option("--common-equivalence-simple", action="store", type="str", help=optparse.SUPPRESS_HELP)
  parser.add_option("--namespace", action="store", type="str")
  parser.add_option("--separate_cmn_hpp", action="store_true", default=False)
  parser.add_option("--separate-cmn-hpp", action="store_true", help=optparse.SUPPRESS_HELP)
  parser.add_option("--number_of_function_files", action="store", type="int", metavar="INTEGER")
  parser.add_option("--number-of-function-files", action="store", type="int", help=optparse.SUPPRESS_HELP)
  parser.add_option("--example", action="store_true", default=False)
  parser.add_option("--debug", action="store_true", default=False)
  parser.add_option("--exe_name", action="store", type="str")
  parser.add_option("--exe-name", action="store", type="str",help=optparse.SUPPRESS_HELP)

  co, files = parser.parse_args(args)
  if co.valgrind: co.run = True
  if co.run: co.link = True
  if co.link: co.compile = True
  if not co.each:
    process(options=co)(file_names=files)
  else:
    from fable.command_line.read import process_each
    process_each(process=process(options=co), file_names=files)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])

