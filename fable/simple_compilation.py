import os
op = os.path

class environment(object):

  __slots__ = [
    "compiler",
    "obj_suffix",
    "exe_suffix",
    "pch_suffix",
    "compiler_path",
    "gcc_version",
    "fable_dist",
    "tbxx_root",
    "__have_pch"]

  def __init__(O, compiler=None):
    if (os.name == "nt"):
      O.compiler = "cl"
      O.obj_suffix = ".obj"
      O.exe_suffix = ".exe"
      O.pch_suffix = None
    else:
      O.compiler = "g++"
      O.obj_suffix = ".o"
      O.exe_suffix = ""
      O.pch_suffix = ".gch"
    if (compiler is not None):
      O.compiler = compiler
    compiler_from_os_environ = os.environ.get("FABLE_COMPILER")
    if (compiler_from_os_environ is not None):
      O.compiler = compiler_from_os_environ
    from libtbx.path import full_command_path
    O.compiler_path = full_command_path(command=O.compiler+O.exe_suffix)
    import libtbx.load_env
    if (O.compiler == "g++" and O.compiler_path is not None):
      O.gcc_version = libtbx.env_config.get_gcc_version(
        command_name=O.compiler)
    else:
      O.gcc_version = None
    O.fable_dist = libtbx.env.dist_path(module_name="fable")
    if (op.isdir(op.join(O.fable_dist, "tbxx"))):
      O.tbxx_root = None
    else:
      O.tbxx_root = op.dirname(libtbx.env.dist_path(module_name="tbxx"))
    O.__have_pch = False

  def set_have_pch(O):
    O.__have_pch = True

  def assemble_command(O,
        link,
        disable_warnings,
        file_names,
        out_name):
    def quote(s):
      assert s.find('"') < 0
      return '"'+s+'"'
    def quote_list(l):
      return " ".join([quote(s) for s in l])
    qon = quote(out_name)
    if (O.compiler == "cl"):
      if (not link): part = "/c /Fo%s" % qon
      else:          part = "/Fe%s" % qon
      def add_to_include_search_path(path):
        if (path is None): return ""
        return " /I%s" % quote(path)
      result = "%s /nologo /EHsc %s%s%s %s" % (
        O.compiler,
        part,
        add_to_include_search_path(O.fable_dist),
        add_to_include_search_path(O.tbxx_root),
        quote_list(file_names))
    else:
      if (not link): opt_c = "-c "
      else:          opt_c = ""
      if (disable_warnings or O.gcc_version < 30400):
        opt_w = "-w"
      else:
        opt_w = "-Wall -Wno-sign-compare -Winvalid-pch"
      if (out_name.endswith(O.pch_suffix)):
        assert not O.__have_pch
        opt_x = " -x c++-header"
      else:
        opt_x = ""
      if (not O.__have_pch):
        def add_to_include_search_path(path):
          if (path is None): return ""
          return " -I%s" % quote(path)
        opt_i = "%s%s" % (
          add_to_include_search_path(O.fable_dist),
          add_to_include_search_path(O.tbxx_root))
      else:
        opt_i = " -I."
      result = "%s -o %s %s%s -g -O0%s%s %s" % (
        O.compiler, qon, opt_c, opt_w, opt_i, opt_x, quote_list(file_names))
    return result

  def file_name_obj(O, file_name_cpp):
    assert file_name_cpp.endswith(".cpp")
    return file_name_cpp[:-4] + O.obj_suffix

  def file_name_exe(O, exe_root):
    return exe_root + O.exe_suffix

  def compilation_command(O, file_name_cpp, disable_warnings=False):
    return O.assemble_command(
      link=False,
      disable_warnings=disable_warnings,
      file_names=[file_name_cpp],
      out_name=O.file_name_obj(file_name_cpp=file_name_cpp))

  def link_command(O, file_names_obj, exe_root):
    return O.assemble_command(
      link=True,
      disable_warnings=False,
      file_names=file_names_obj,
      out_name=O.file_name_exe(exe_root=exe_root))

  def build(O,
        link,
        file_name_cpp,
        obj_name=None,
        exe_name=None,
        pch_name=None,
        disable_warnings=False,
        show_command=False,
        Error=RuntimeError):
    assert [obj_name, exe_name, pch_name].count(None) >= 2
    if (link):
      out_name = exe_name
      out_suffix = O.exe_suffix
    elif (pch_name is None):
      out_name = obj_name
      out_suffix = O.obj_suffix
    else:
      assert O.pch_suffix is not None
      out_name = pch_name + O.pch_suffix
      out_suffix = None
    if (out_name is None):
      assert file_name_cpp.endswith(".cpp")
      out_name = file_name_cpp[:-4] + out_suffix
    from libtbx.utils import remove_files
    remove_files(out_name)
    cmd = O.assemble_command(
      link=link,
      disable_warnings=disable_warnings,
      file_names=[file_name_cpp],
      out_name=out_name)
    if (show_command):
      print cmd
    from libtbx import easy_run
    buffers = easy_run.fully_buffered(command=cmd)
    if (O.compiler != "cl" or buffers.stderr_lines != [file_name_cpp]):
      buffers.raise_if_errors(Error=Error)
    return out_name
