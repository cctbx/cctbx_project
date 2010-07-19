import fable.cout

def compute_hexdigest(text):
  from libtbx.utils import hashlib_md5
  m = hashlib_md5()
  m.update(text)
  return m.hexdigest()

def check_fingerprint(file_name):
  from libtbx.utils import hashlib_md5
  lines = open(file_name).read().splitlines()
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
  import os
  op = os.path
  if (op.exists(file_name)):
    if (not op.isfile(file_name)):
      raise RuntimeError(
        "Not a regular file: %s" % show_string(file_name))
    stat = check_fingerprint(file_name=file_name)
    if (stat is None or not stat):
      raise RuntimeError(
        "File appears to be manually modified: %s" % show_string(file_name))
  hexdigest = compute_hexdigest(text=text)
  f = open(file_name, "w")
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
    import sys
    if (O.n_calls != 0):
      print
    O.n_calls += 1
    opts = O.options
    lines = fable.cout.process(
      file_names=file_names,
      top_unit_name=opts.top_unit_name,
      include_guard_suffix=opts.include_guard_suffix,
      dynamic_parameters=O.dynamic_parameters,
      fem_do_safe=not opts.no_fem_do_safe,
      arr_nd_size_max=opts.arr_nd_size_max,
      inline_all=opts.inline_all,
      common_equivalence_simple=set(opts.common_equivalence_simple.split(",")),
      namespace=opts.namespace,
      separate_cmn_hpp=opts.separate_cmn_hpp,
      number_of_function_files=opts.number_of_function_files,
      debug=opts.debug)
    text = "\n".join(lines)+"\n"
    if (opts.top_unit_name is None or not opts.debug):
      sys.stdout.write(text)
    if (len(file_names) != 0 and opts.compile):
      print
      write_only_if_safe(file_name="fable_cout.cpp", text=text)
      from fable import simple_compilation
      comp_env = simple_compilation.environment()
      out_name = comp_env.build(
        link=opts.link, file_name_cpp="fable_cout.cpp", show_command=True)
      print
      if (opts.run):
        from libtbx import easy_run
        import os
        op = os.path
        cmd = op.join(".", out_name)
        if (opts.valgrind):
          cmd = "valgrind " + cmd
        print cmd
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
  from libtbx.option_parser import option_parser
  command_line = (option_parser(
    usage="%s [options] fortran_file ..." % libtbx.env.dispatcher_name)
    .option(None, "--compile",
      action="store_true",
      default=False)
    .option(None, "--link",
      action="store_true",
      default=False)
    .option(None, "--run",
      action="store_true",
      default=False)
    .option(None, "--valgrind",
      action="store_true",
      default=False)
    .option(None, "--each",
      action="store_true",
      default=False)
    .option(None, "--top_unit_name",
      action="store",
      type="str",
      metavar="IDENTIFIER")
    .option(None, "--include_guard_suffix",
      action="store",
      type="str",
      metavar="STRING")
    .option(None, "--dynamic_parameter",
      action="append",
      type="str",
      metavar="STRING",
      help='example: --dynamic_parameter="int array_size=100"')
    .option(None, "--no_fem_do_safe",
      action="store_true",
      default=False)
    .option(None, "--arr_nd_size_max",
      action="store",
      type="int",
      default=fable.cout.default_arr_nd_size_max,
      metavar='INTEGER (default: %d)' % fable.cout.default_arr_nd_size_max)
    .option(None, "--inline_all",
      action="store_true",
      default=False)
    .option(None, "--common_equivalence_simple",
      action="store",
      type="str",
      default="",
      metavar="STRING",
      help='comma-separated list of common names')
    .option(None, "--namespace",
      action="store",
      type="str")
    .option(None, "--separate_cmn_hpp",
      action="store_true",
      default=False)
    .option(None, "--number_of_function_files",
      action="store",
      type="int",
      metavar="INTEGER")
    .option(None, "--example",
      action="store_true",
      default=False)
    .option(None, "--debug",
      action="store_true",
      default=False)
  ).process(args=args)
  co = command_line.options
  if (co.valgrind): co.run = True
  if (co.run): co.link = True
  if (co.link): co.compile = True
  if (not co.each):
    process(options=co)(file_names=command_line.args)
  else:
    from fable.command_line.read import process_each
    process_each(process=process(options=co), file_names=command_line.args)

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
