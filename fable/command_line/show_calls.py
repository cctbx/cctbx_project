def run(args):
  if (len(args) == 0): args = ["--help"]
  from libtbx.option_parser import option_parser
  import libtbx.load_env
  command_line = (option_parser(
    usage="%s [options] fortran_file ..." % libtbx.env.dispatcher_name)
    .option(None, "--top_procedure",
      action="append",
      type="str")
    .option(None, "--write_graphviz_dot",
      action="store",
      type="str")
  ).process(args=args)
  co = command_line.options
  from fable.read import process
  all_fprocs = process(file_names=command_line.args)
  topological_fprocs = all_fprocs.build_bottom_up_fproc_list_following_calls(
    top_procedures=co.top_procedure)
  dep_cycles = topological_fprocs.dependency_cycles
  if (len(dep_cycles) != 0):
    print "Dependency cycles:", len(dep_cycles)
    for cycle in dep_cycles:
      print " ", " ".join(cycle)
    print
  print "Top-down procedure list:"
  print
  digraph_lhs_rhs = []
  for fproc in reversed(topological_fprocs.bottom_up_list):
    if (fproc.name is None):
      lhs = fproc.fproc_type
      print lhs
    else:
      lhs = fproc.name.value
      print fproc.fproc_type, fproc.name.value
    fwds = set(
      topological_fprocs.forward_uses_by_identifier.get(
        fproc.name.value, []))
    for identifier in sorted(fproc.fdecl_by_identifier.keys()):
      fdecl = fproc.fdecl_by_identifier[identifier]
      if (fdecl.is_fproc_name()): continue
      if (not fdecl.is_user_defined_callable()):
        continue
      called_name = fdecl.id_tok.value
      passed = fproc.externals_passed_by_arg_identifier.get(called_name)
      if (passed is None):
        digraph_lhs_rhs.append((lhs, called_name))
      else:
        called_name += "->" + ",".join(sorted(passed))
        for indirectly_called_name in passed:
          digraph_lhs_rhs.append((lhs, indirectly_called_name))
      if (fdecl.is_function()):
        sz = ""
        if (fdecl.size_tokens is not None):
          if (len(fdecl.size_tokens) == 1
                and fdecl.size_tokens[0].is_integer()):
            sz = "*%s" % fdecl.size_tokens[0].value
          else:
            sz = "*(*)"
        s = "%s (%s%s)" % (called_name, fdecl.data_type.value, sz)
      else:
        s = called_name
      if (called_name in fwds):
        s += " (dependency cycle)"
      print "  %s" % s
  print
  if (co.write_graphviz_dot is not None):
    f = open(co.write_graphviz_dot, "w")
    print >> f, "digraph G {"
    for lhs_rhs in digraph_lhs_rhs:
      print >> f, "  %s -> %s;" % lhs_rhs
    print >> f, "}"
    del f

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
