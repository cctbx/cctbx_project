import sys, os
op = os.path

def run(args=[]):
  assert len(args) == 0
  import libtbx.load_env
  src_dir = libtbx.env.under_dist(
    module_name="scitbx", path="lbfgs", test=op.isdir)
  import fable.read
  all_units = fable.read.process(
    file_names=[op.join(src_dir, f) for f in ["sdrive.f", "lbfgs.f"]])
  namespace = "scitbx::lbfgs_fem"
  functions_public = set(["lbfgs", "blockdata_lb2"])
  functions_detail = set(["lb1", "daxpy", "ddot", "mcstep", "mcsrch"])
  functions_program = set(["one_pass"])
  import fable.cout
  functions_hpp = fable.cout.process(
    all_units=all_units,
    namespace=namespace,
    fem_do_safe=False,
    suppress_program=True,
    suppress_common=False,
    suppress_functions=functions_detail.union(functions_program),
    suppress_function_definitions=functions_public)
  functions_cpp = fable.cout.process(
    all_units=all_units,
    namespace=namespace,
    fem_do_safe=False,
    suppress_program=True,
    suppress_common=True,
    suppress_functions=functions_program)
  functions_cpp[0] = "#include <scitbx/lbfgs_fem.hpp>"
  sdrive_cpp = fable.cout.process(
    all_units=all_units,
    namespace=namespace,
    fem_do_safe=False,
    suppress_common=True,
    suppress_functions=functions_detail.union(functions_public))
  sdrive_cpp[0] = functions_cpp[0]
  #
  def make_target_dir(path):
    result = libtbx.env.under_build(path=path)
    if (not op.isdir(result)):
      os.makedirs(result)
      assert op.isdir(result)
    return result
  target_dir = make_target_dir(path="include/scitbx")
  print >> open(op.join(target_dir, "lbfgs_fem.hpp"), "w"), \
    "\n".join(functions_hpp)
  target_dir = make_target_dir(path="scitbx/lbfgs")
  print >> open(op.join(target_dir, "lbfgs_fem.cpp"), "w"), \
    "\n".join(functions_cpp)
  print >> open(op.join(target_dir, "sdrive_fem.cpp"), "w"), \
    "\n".join(sdrive_cpp)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
