# $Id$

import sys, os
from string import split, strip, maketrans, translate, join, upper
from make import expand_cf
transl_table_slash_backslash = maketrans("/", "\\")

class file_emitter:

  def __init__(self,obj):
    self.obj = obj
    self.emitter_tuples = []
    if type(obj.prefix_macro) == type("string"):
      for f in obj.files:
        self.emitter_tuples.append((obj.prefix_macro,f))
    elif type(obj.prefix_macro) == type([]):
      for x in xrange(len(obj.prefix_macro)):
        for f in obj.files[x]:
          self.emitter_tuples.append((obj.prefix_macro[x],f))
    self.length = len(self.emitter_tuples)

  def __len__(self):
    return self.length

  def __getitem__(self, num):
    try:
      return self.emitter_tuples[num]
    except:
      raise IndexError


class write_makefiles:

  def __init__(self, subdir, configuration, package):
    self.platform = strip(configuration[0])
    if (not (self.platform in ("tru64_cxx", "unix_gcc", "sun_gcc", "irix_CC",
                               "macosx", "mingw32", "vc60", "win32_mwcc"))):
      stdout = sys.stdout
      sys.stdout = sys.__stdout__
      print "*" * 78
      print 'WARNING: Unsupported platform identifier: "' \
            + self.platform + '"'
      print "If the installation fails you might have to modify the file"
      print "cctbx/build/makefile_generator.py"
      print "If there are error messages, study them carefully."
      print "The installation will proceed in 10 seconds."
      print "*" * 78
      sys.stdout = stdout
      import time
      time.sleep(10)

    self.macros = configuration[1:]
    self.package = package
    # remove empty lines at beginning
    while (len(strip(self.macros[0])) == 0): del self.macros[0]
    self.dependencies()

    if (not hasattr(self, "prefix_macro")):
      self.prefix_macro = upper(self.package.name)
    self.file_dependency_emitter = file_emitter(self)

    if subdir in self.package.tbx_subpkg:
      expand_cf(self.macros, "@(SUBLEVEL1)", os.sep+subdir)
    else:
      expand_cf(self.macros, "@(SUBLEVEL1)", "")

    self.build_shared_libraries = 0
    if (self.platform in ("irix_CC", "unix_gcc", "sun_gcc")):
      self.build_shared_libraries = 1

  def head(self):
    print r"""# Usage:
#
# Unix:
#
#   make softlinks     Create softlinks to source code and tests
#   make               Compile all sources
#   make clean         Remove all object files
#   make unlink        Remove softlinks
#
# Development option:
#   Rename this Makefile to Makefile.nodepend.
#   Run "make -f Makefile.nodepend depend > Makefile"
#   to automatically generate all dependencies.
#
# Windows:
#
#   make copy          Copy source code and tests
#   make               Compile all sources
#   make clean         Remove all object files
#   make del           Remove source code and tests
"""

    self.make_targets = {
      "libraries": [],
      "executables": [],
      "boost_python_modules": [],
      "depend": [],
      "clean": [],
    }
    for m in self.macros: print m
    print
    print "all: compile"
    print

  def format_objects(self, objects):
    if (self.platform in ("vc60", "win32_mwcc")):
      doto = ".obj"
    else:
      doto = ".o"
    s = ""
    for obj in objects:
      s = s + " " + obj + doto
    return s[1:]

  def format_libs(self, libs, macros):
    s = ""
    if (len(libs)):
      if (self.platform in ("vc60", "win32_mwcc")):
        s = s + " $(SP_LIBDIR_WIN)\\lib*.lib"
      else:
        s = s + " -L$(SP_LIBDIR_UNIX)"
        if (self.platform != "irix_CC"):
          for l in libs: s = s + " -l" + l
        else:
          for l in libs: s = s + " -exports -l" + l
    if (len(macros) > 0):
      s = s + " " + join(macros)
    return strip(s)

  def tail(self):

    all = []
    for t in ("libraries", "executables", "boost_python_modules"):
      s = join(self.make_targets[t])
      if (len(s)):
        print "%s: %s" % (t, s)
        all.append(t)
    if (len(all)):
      print "compile:", join(all)
    else:
      print "compile:"
    print

    if (hasattr(self, "make_test")):
      self.make_test()

    doto = self.format_objects(("",))
    print "CPPOPTS=$(STDFIXINC) $(STDOPTS) $(WARNOPTS) $(OPTOPTS) \\"
    print "        ",
    for pk in [upper(x) for x in (self.package.name,)+self.package.supporting]:
      print "$(%sINC) "%(pk),
    print "\\"
    print "        $(BOOSTINC) $(PYINC)"
    print
    print ".SUFFIXES: %s .cpp .c" % (doto,)
    print
    print ".cpp%s:" % (doto,)
    print "\t$(CPP) $(CPPOPTS) -c $*.cpp"
    print
    print ".c%s:" % (doto,)
    print "\t$(CC) $(CPPOPTS) -c $*.c"
    print

    self.make_clean()

    if (self.platform != "vc60"):
      print "depend:"
      if (self.platform in ("mingw32", "win32_mwcc")):
        print "\t@type Makefile.nodepend"
      else:
        print "\t@cat Makefile.nodepend"
      for src in self.make_targets["depend"]:
        print "\t@$(CPP) $(CPPOPTS) $(MAKEDEP) %s.cpp" % (src,)
      print

  def file_management(self):
    print "softlinks:"
    for srcf in self.file_dependency_emitter:
      print "\t-ln -s $(%s_UNIX)/%s ." % srcf
    print
    print "cp:"
    for srcf in self.file_dependency_emitter:
      print "\t-cp $(%s_UNIX)/%s ." % srcf
    print
    print "unlink:"
    for srcf in self.file_dependency_emitter:
      f = split(srcf[1], "/")[-1]
      print "\t-test -L %s && rm %s" % (f, f)
    print
    print "rm:"
    for srcf in self.file_dependency_emitter:
      print "\t-rm " + split(srcf[1], "/")[-1]
    print
    if (self.platform in ("mingw32", "vc60", "win32_mwcc")):
      print "copy:"
      for srcf in self.file_dependency_emitter:
        f = translate(srcf[1], transl_table_slash_backslash)
        print "\t-copy $(%s_WIN)\\%s" % (srcf[0], f)
      print
      print "del:"
      for srcf in self.file_dependency_emitter:
        print "\t-del " + split(srcf[1], "/")[-1]
      print

  def update_depend(self, objects):
    for obj in objects:
      if (not obj in self.make_targets["depend"]):
        self.make_targets["depend"].append(obj)

  def make_library(self, name, objects):
    objstr = self.format_objects(objects)
    if (self.platform == "vc60"):
      lib = "lib" + name + ".lib"
      print "%s: %s" % (lib, objstr)
      print "\t-del %s" % (lib,)
      print "\t$(LD) -lib /nologo /out:%s %s" % (lib, objstr)
    elif (self.platform == "win32_mwcc"):
      lib = "lib" + name + ".lib"
      print "%s: %s" % (lib, objstr)
      print "\t-del %s" % (lib,)
      print "\t$(LD) -library -o %s %s" % (lib, objstr)
    else:
      lib_suffix = ".a"
      if (self.build_shared_libraries):
        lib_suffix = ".so"
      lib = "lib" + name + lib_suffix
      print "%s: %s" % (lib, objstr)
      if (self.platform == "mingw32"):
        print "\t-del %s" % (lib,)
      else:
        print "\trm -f %s" % (lib,)
      if (self.build_shared_libraries):
        print "\t$(CPP) $(LDDLL) -o %s %s" % (lib, objstr)
      elif (self.platform == "tru64_cxx"):
        print "\tcd cxx_repository; \\"
        print "\t  ls -1 > ../%s.input; \\" % (lib,)
        print "\t  ar r ../%s -input ../%s.input" % (lib, lib)
        print "\trm -f %s.input" % (lib,)
        print "\tar r %s %s" % (lib, objstr)
      elif (self.platform == "macosx"):
        print "\tlibtool -static -o %s %s" % (lib, objstr)
      else:
        print "\tar r %s %s" % (lib, objstr)
    if (self.platform  in ("vc60", "mingw32", "win32_mwcc")):
      print "\tcopy %s $(SP_LIBDIR_WIN)" % (lib,)
    else:
      print "\tcp %s $(SP_LIBDIR_UNIX)" % (lib,)
    print
    self.make_targets["libraries"].append(lib)
    self.update_depend(objects)

  def make_executable(self, name, objects_and_libs):
    objects = objects_and_libs[0]
    libs = objects_and_libs[1]
    objstr = self.format_objects(objects)
    libstr = self.format_libs(libs, ("$(LDMATH)",))
    if (not self.platform in ("mingw32", "vc60", "win32_mwcc")):
      nameexe = name
    else:
      nameexe = name + ".exe"
    if (self.platform != "vc60"):
      out = "-o "
    else:
      out = "/out:"
    print "%s: %s" % (nameexe, objstr)
    print "\t$(LD) $(LDEXE) %s %s%s %s" % (objstr, out, nameexe, libstr)
    print
    self.make_targets["executables"].append(nameexe)
    self.make_targets["clean"].append(nameexe)
    self.update_depend(objects)

  def make_boost_python_module(self, name, objects_and_libs):
    objects = objects_and_libs[0]
    libs = objects_and_libs[1] + ("boost_python",)
    objstr = self.format_objects(objects)
    if   (self.platform == "mingw32"):
      self.mingw32_pyd(name, objstr, libs)
    elif (self.platform == "vc60"):
      self.vc60_pyd(name, objstr, libs)
    elif (self.platform == "win32_mwcc"):
      self.win32_mwcc_pyd(name, objstr, libs)
    else:
      self.unix_so(name, objstr, libs)
    print
    self.update_depend(objects)

  def unix_so(self, name, objstr, libs):
    libstr = self.format_libs(libs, ("$(PYLIB)", "$(LDMATH)"))
    nameso = name + ".so"
    print "%s: %s" % (nameso, objstr)
    print "\t$(LD) $(LDDLL) -o %s %s %s" \
          % (nameso, objstr, libstr)
    print "\tcp %s $(%sLIB_PYTHONDIR_UNIX)" % (nameso,upper(self.package.name))
    print
    self.make_targets["boost_python_modules"].append(nameso)

  def mingw32_pyd(self, name, objstr, libs):
    libstr = self.format_libs(libs, ("$(PYLIB)",))
    namepyd = name + ".pyd"
    namedef = name + ".def"
    print "%s: %s %s" % (namepyd, namedef, objstr)
    print (  "\tdllwrap -s --driver-name g++ --entry _DllMainCRTStartup@12"
           + " --target=i386-mingw32 --dllname %s --def %s"
           + " %s %s") % (namepyd, namedef, objstr, libstr)
    print "\tcopy %s $(%sLIB_PYTHONDIR_WIN)" % (namepyd,upper(self.package.name))
    print
    print "%s:" % (namedef,)
    print "\techo EXPORTS > %s" % (namedef,)
    print "\techo \tinit%s >> %s" % (name, namedef)
    print
    self.make_targets["boost_python_modules"].append(namepyd)

  def vc60_pyd(self, name, objstr, libs):
    libstr = self.format_libs(libs, ("$(PYLIB)",))
    namepyd = name + ".pyd"
    print "%s: %s" % (namepyd, objstr)
    print (  "\t$(LD) $(LDDLL) /out:%s /export:init%s %s"
           + " %s") % (namepyd, name, objstr, libstr)
    print "\tcopy %s $(%sLIB_PYTHONDIR_WIN)" % (namepyd,upper(self.package.name))
    self.make_targets["boost_python_modules"].append(namepyd)

  def win32_mwcc_pyd(self, name, objstr, libs):
    libstr = self.format_libs(libs, ("$(PYLIB)",))
    namepyd = name + ".pyd"
    print "%s: %s" % (namepyd, objstr)
    print (  "\t$(LD) $(LDDLL) -o %s %s %s") % (
      namepyd, objstr, libstr)
    print "\tcopy %s $(%sLIB_PYTHONDIR_WIN)" % (namepyd,upper(self.package.name))
    self.make_targets["boost_python_modules"].append(namepyd)

  def make_clean(self):
    print "clean_unix:"
    for f in self.make_targets["clean"]:
      print "\trm -f " + f
    print "\trm -f *.o *.a *.so *.pyc"
    print "\trm -f *.obj *.lib *.exp *.idb *.exe *.def *.pyd"
    print "\trm -rf cxx_repository so_locations ii_files"
    print
    print "clean_win:"
    for f in self.make_targets["clean"]:
      print "\t-del " + f
    for ext in ("o", "a", "so", "pyc",
                "obj", "lib", "exp", "idb", "exe", "def", "pyd"):
      print "\t-del *." + ext
    print
    if (self.platform in ("mingw32", "vc60", "win32_mwcc")):
      print "clean: clean_win"
    else:
      print "clean: clean_unix"
    print

  def write(self, file):
    old_sys_stdout = sys.stdout
    sys.stdout = file
    try:
      self.head()
      if (hasattr(self, "libraries")):
        for name in self.libraries.keys():
          self.make_library(name, self.libraries[name])
      if (hasattr(self, "executables")):
        for name in self.executables.keys():
          self.make_executable(name, self.executables[name])
      if (hasattr(self, "boost_python_modules")):
        for name in self.boost_python_modules.keys():
          self.make_boost_python_module(name, self.boost_python_modules[name])
      self.file_management()
      self.tail()
    finally:
      sys.stdout = old_sys_stdout
