# $Id$

import sys, os
from string import split, strip, maketrans, translate, join
transl_table_slash_backslash = maketrans("/", "\\")

class write_makefiles:

  def __init__(self, configuration):
    self.platform = strip(configuration[0])
    if (not (self.platform in ("tru64_cxx", "linux_gcc", "irix_CC",
                               "mingw32", "vc60"))):
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
    # remove empty lines at beginning
    while (len(strip(self.macros[0])) == 0): del self.macros[0]
    self.dependencies()

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
    if (self.platform == "vc60"):
      doto = ".obj"
    else:
      doto = ".o"
    s = ""
    for obj in objects:
      s = s + " " + obj + doto
    return s[1:]

  def tail(self):

    all = []
    for t in ("libraries", "executables", "boost_python_modules"):
      s = join(self.make_targets[t])
      if (len(s)):
        print "%s: %s" % (t, s)
        all.append(t)
    if (len(all)):
      print "compile:", join(all)
    print

    if (hasattr(self, "make_test")):
      self.make_test()

    doto = self.format_objects(("",))
    print "CPPOPTS=$(STDFIXINC) $(STDOPTS) $(WARNOPTS) $(OPTOPTS) \\"
    print "        $(CCTBXINC) $(BOOSTINC) $(PYINC)"
    print
    print ".SUFFIXES: %s .cpp" % (doto,)
    print
    print ".cpp%s:" % (doto,)
    print "\t$(CPP) $(CPPOPTS) -c $*.cpp"
    print

    self.make_clean()

    if (self.platform != "vc60"):
      print "depend:"
      if (self.platform == "mingw32"):
        print "\t@type Makefile.nodepend"
      else:
        print "\t@cat Makefile.nodepend"
      for src in self.make_targets["depend"]:
        print "\t@$(CPP) $(CPPOPTS) $(MAKEDEP) %s.cpp" % (src,)
      print

  def file_management(self):
    print "softlinks:"
    for srcf in self.files:
      print "\t-ln -s $(CCTBX_UNIX)/" + srcf + " ."
    print
    print "cp:"
    for srcf in self.files:
      print "\t-cp $(CCTBX_UNIX)/" + srcf + " ."
    print
    print "unlink:"
    for srcf in self.files:
      f = split(srcf, "/")[-1]
      print "\t-test -L %s && rm %s" % (f, f)
    print
    print "rm:"
    for srcf in self.files:
      print "\t-rm " + split(srcf, "/")[-1]
    print
    if (self.platform in ("mingw32", "vc60")):
      print "copy:"
      for srcf in self.files:
        f = translate(srcf, transl_table_slash_backslash)
        print "\t-copy $(CCTBX_WIN)\\" + f
      print
      print "del:"
      for srcf in self.files:
        print "\t-del " + split(srcf, "/")[-1]
      print

  def update_depend(self, objects):
    for obj in objects:
      if (not obj in self.make_targets["depend"]):
        self.make_targets["depend"].append(obj)

  def make_library(self, name, objects):
    objstr = self.format_objects(objects)
    if (self.platform != "vc60"):
      lib = "lib" + name + ".a"
      print "%s: %s" % (lib, objstr)
      if (self.platform == "mingw32"):
        print "\t-del %s" % (lib,)
      else:
        print "\trm -f %s" % (lib,)
      if   (self.platform == "tru64_cxx"):
        print "\tar r %s %s cxx_repository/*.o" % (lib, objstr)
      elif (self.platform == "irix_CC"):
        print "\t$(CPP) -ar -o %s %s" % (lib, objstr)
      else:
        print "\tar r %s %s" % (lib, objstr)
    else:
      lib = "lib" + name + ".lib"
      print "%s: %s" % (lib, objstr)
      print "\t-del %s" % (lib,)
      print "\t$(LD) -lib /nologo /out:%s %s" % (lib, objstr)
    print
    self.make_targets["libraries"].append(lib)
    self.update_depend(objects)

  def make_executable(self, name, objects, libs = "$(LDMATH)"):
    objstr = self.format_objects(objects)
    if (not self.platform in ("mingw32", "vc60")):
      nameexe = name
    else:
      nameexe = name + ".exe"
    if (self.platform != "vc60"):
      out = "-o "
    else:
      out = "/out:"
    print "%s: %s" % (nameexe, objstr)
    print "\t$(LD) $(LDEXE) %s %s%s %s" % (objstr, out, nameexe, libs)
    print
    self.make_targets["executables"].append(nameexe)
    self.make_targets["clean"].append(nameexe)
    self.update_depend(objects)

  def make_boost_python_module(self, name, objects):
    objstr = self.format_objects(objects)
    if   (self.platform == "mingw32"):
      self.mingw32_pyd(name, objstr)
    elif (self.platform == "vc60"):
      self.vc60_pyd(name, objstr)
    else:
      self.unix_so(name, objstr)
    print
    self.update_depend(objects)

  def unix_so(self, name, objstr):
    nameso = name + ".so"
    print "%s: %s" % (nameso, objstr)
    print "\t$(LD) $(LDDLL) -o %s %s $(BOOST_PYTHONLIB) $(PYLIB) $(LDMATH)" \
          % (nameso, objstr)
    print
    self.make_targets["boost_python_modules"].append(nameso)

  def mingw32_pyd(self, name, objstr):
    namepyd = name + ".pyd"
    namedef = name + ".def"
    print "%s: %s %s" % (namepyd, namedef, objstr)
    print (  "\tdllwrap -s --driver-name g++ --entry _DllMainCRTStartup@12"
           + " --target=i386-mingw32 --dllname %s --def %s"
           + " %s $(BOOST_PYTHONLIB) $(PYLIB)") % (namepyd, namedef, objstr)
    print
    print "%s:" % (namedef,)
    print "\techo EXPORTS > %s" % (namedef,)
    print "\techo \tinit%s >> %s" % (name, namedef)
    print
    self.make_targets["boost_python_modules"].append(namepyd)

  def vc60_pyd(self, name, objstr):
    namepyd = name + ".pyd"
    print "%s: %s" % (namepyd, objstr)
    print (  "\t$(LD) $(LDDLL) /out:%s /export:init%s %s"
           + " $(BOOST_PYTHONLIB) $(PYLIB)") % (namepyd, name, objstr)
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
    if (self.platform in ("mingw32", "vc60")):
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
      if (hasattr(self, "examples")):
        for name in self.examples.keys():
          self.make_executable(name, self.examples[name],
                               "$(CCTBXLIB) $(LDMATH)")
      if (hasattr(self, "boost_python_modules")):
        for name in self.boost_python_modules.keys():
          self.make_boost_python_module(name, self.boost_python_modules[name])
      self.file_management()
      self.tail()
    finally:
      sys.stdout = old_sys_stdout
