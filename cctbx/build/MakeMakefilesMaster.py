# $Id$

import sys, os
from string import split, maketrans, translate
transl_table_slash_backslash = maketrans("/", "\\")

class write_makefiles:

  def head(self, target):
    self.all = []
    self.depend = []
    self.clean = []
    print r"""# Automatically generated file.
# To customize, edit the build/MakeMakefilesMaster.py and MakeMakefiles.py.
# Run "python MakeMakefiles.py" to create the customized Makefiles
# for all platforms.

# Usage:
#
#   Create a new empty directory anywhere (preferably not in the cctbx tree).
#   Copy this Makefile to that new directory and rename it to "Makefile"
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

    if (not target in ("mingw32", "vc60")):
      print r"""
CCTBX_UNIX=/net/cci/rwgk/cctbx
CCTBXINC=-I$(CCTBX_UNIX)
"""
    else:
      print r"""
CCTBX_UNIX=/net/cci/rwgk/cctbx
CCTBX_WIN=L:\cctbx
CCTBXINC=-I$(CCTBX_WIN)
"""

    if (target == "vc60"):
      print r"""
CCTBXLIB=..\..\lib\lib*.lib
"""
    else:
      print r"""
CCTBXLIB=-L../../lib -leltbx -lsgtbx -luctbx
"""

    if (not target in ("mingw32", "vc60")):
      print r"""
BOOSTINC=-I/net/cci/rwgk/boost
"""
    else:
      print r"""
BOOSTINC=-I"L:\boost"
"""

    if   (target == "tru64_cxx"):
      print r"""
BOOST_PYTHONLIB=/net/anaconda/scratch1/rwgk/boost/libboost_python.a
"""
    elif (target == "linux_gcc"):
      print r"""
BOOST_PYTHONLIB=/net/legless/scratch1/rwgk/boost/libboost_python.a
"""
    elif (target == "irix_CC"):
      print r"""
BOOST_PYTHONLIB=/net/rattler/scratch1/rwgk/boost/libboost_python.a
"""
    elif (target == "mingw32"):
      print r"""
BOOST_PYTHONLIB="L:\mingw32\boost\libboost_python.a"
"""
    elif (target == "vc60"):
      print r"""
BOOST_PYTHONLIB="L:\vc60\boost\boost_python.lib"
"""

    if   (target == "mingw32"):
      print r"""
PYEXE="C:\Program files\Python\python.exe"
PYINC=-I"C:\usr\include\python1.5"
PYLIB="C:\usr\lib\libpython15.a"
"""
    elif (target == "vc60"):
      print r"""
PYEXE="C:\Program files\Python\python.exe"
PYINC=-I"C:\Program files\Python\include"
PYLIB="C:\Program files\Python\libs\python15.lib"
"""
    else:
      print r"""
PYEXE=/usr/local/Python-1.5.2/bin/python
PYINC=-I/usr/local/Python-1.5.2/include/python1.5
#PYEXE=/usr/local/Python-2.0/bin/python
#PYINC=-I/usr/local/Python-2.0/include/python2.0
PYLIB=
"""

    if (target in ("tru64_cxx", "irix_CC")):
      print r"""
STDFIXINC=-I/net/cci/xp/C++_C_headers
"""

    if   (target == "tru64_cxx"):
      print r"""
CPP=cxx
LD=cxx
STDOPTS=-std strict_ansi
WARNOPTS=-msg_display_number -msg_disable 450
OPTOPTS=-g
MAKEDEP=-Em
LDEXE=
LDDLL=-shared -expect_unresolved 'Py*' -expect_unresolved '_Py*'
LDMATH=-lm
"""
    elif (target == "linux_gcc"):
      print r"""
CPP=g++
LD=g++
STDOPTS=-ftemplate-depth-21
WARNOPTS=
OPTOPTS=-g
MAKEDEP=-M
LDEXE=
LDDLL=-shared
LDMATH=-lm
"""
    elif (target == "irix_CC"):
      print r"""
CPP=CC -LANG:std -n32 -mips4
LD=CC -LANG:std -n32 -mips4
STDOPTS=
WARNOPTS=-woff 1001,1234,1682
OPTOPTS=-g
MAKEDEP=-M
LDEXE=-LD_MSG:off=15,84
LDDLL=-shared -LD_MSG:off=15,84
LDMATH=-lm
"""
    elif (target == "mingw32"):
      print r"""
CPP=g++
LD=g++
STDOPTS=-ftemplate-depth-21
WARNOPTS=
OPTOPTS=-g
MAKEDEP=-M
LDEXE=
LDDLL=
LDMATH=-lm
"""
    elif (target == "vc60"):
      print r"""
CPP=cl.exe
LD=link.exe
STDOPTS=/nologo /MD /GR /GX /Zm200
WARNOPTS=
OPTOPTS=
LDEXE=/nologo /incremental:no
LDDLL=/nologo /incremental:no /dll
LDMATH=
"""

    print "all: all_"
    print

  def tail(self, toolbox, target):

    if (len(self.all) != 0):
      s = "all_:"
      for t in self.all: s = s + " " + t
      print s
      print

    if (hasattr(self, "make_test")):
      self.make_test()

    doto = format_objects(target, ("",))
    print r"""
CPPOPTS=$(STDFIXINC) $(STDOPTS) $(WARNOPTS) $(OPTOPTS) \
        $(BOOSTINC) $(PYINC) $(CCTBXINC)

.SUFFIXES: %s .cpp

.cpp%s:
	$(CPP) $(CPPOPTS) -c $*.cpp
""" % (doto, doto)

    if (target in ("mingw32", "vc60")):
      print r"""
copy_makefile:
	copy $(CCTBX_WIN)\%s\Makefile.%s Makefile
""" % (translate(toolbox, transl_table_slash_backslash), target)

    self.make_clean(target)

    if (target != "vc60"):
      print "depend:"
      if (target == "mingw32"):
        print "\t@type Makefile.nodepend"
      else:
        print "\t@cat Makefile.nodepend"
      for src in self.depend:
        print "\t@$(CPP) $(CPPOPTS) $(MAKEDEP) %s.cpp" % (src,)
      print

  def file_management(self, target):
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
    if (target in ("mingw32", "vc60")):
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
      if (not obj in self.depend):
        self.depend.append(obj)

  def make_library(self, target, name, objects):
    objstr = format_objects(target, objects)
    if (target != "vc60"):
      lib = "lib" + name + ".a"
      print "%s: %s" % (lib, objstr)
      if (target == "mingw32"):
        print "\t-del %s" % (lib,)
      else:
        print "\trm -f %s" % (lib,)
      if   (target == "tru64_cxx"):
        print "\tar r %s %s cxx_repository/*.o" % (lib, objstr)
      elif (target == "irix_CC"):
        print "\t$(CPP) -ar -o %s %s" % (lib, objstr)
      else:
        print "\tar r %s %s" % (lib, objstr)
    else:
      lib = "lib" + name + ".lib"
      print "%s: %s" % (lib, objstr)
      print "\t-del %s" % (lib,)
      print "\t$(LD) -lib /nologo /out:%s %s" % (lib, objstr)
    print
    self.all.append(lib)
    self.update_depend(objects)

  def make_executable(self, target, name, objects, libs = "$(LDMATH)"):
    objstr = format_objects(target, objects)
    if (not target in ("mingw32", "vc60")):
      nameexe = name
    else:
      nameexe = name + ".exe"
    if (target != "vc60"):
      out = "-o "
    else:
      out = "/out:"
    print "%s: %s" % (nameexe, objstr)
    print "\t$(LD) $(LDEXE) %s %s%s %s" % (objstr, out, nameexe, libs)
    print
    self.all.append(nameexe)
    self.update_depend(objects)
    self.clean.append(nameexe)

  def make_boost_python_module(self, target, name, objects):
    objstr = format_objects(target, objects)
    if   (target == "mingw32"):
      self.mingw32_pyd(name, objstr)
    elif (target == "vc60"):
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
    self.all.append(nameso)

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
    self.all.append(namepyd)

  def vc60_pyd(self, name, objstr):
    namepyd = name + ".pyd"
    print "%s: %s" % (namepyd, objstr)
    print (  "\t$(LD) $(LDDLL) /out:%s /export:init%s %s"
           + " $(BOOST_PYTHONLIB) $(PYLIB)") % (namepyd, name, objstr)
    self.all.append(namepyd)

  def make_clean(self, target):
    print "clean_unix:"
    for f in self.clean:
      print "\trm -f " + f
    print "\trm -f *.o *.a *.so *.pyc"
    print "\trm -f *.obj *.lib *.exp *.idb *.exe *.def *.pyd"
    print "\trm -rf cxx_repository so_locations ii_files"
    print
    print "clean_win:"
    for f in self.clean:
      print "\t-del " + f
    for ext in ("o", "a", "so", "pyc",
                "obj", "lib", "exp", "idb", "exe", "def", "pyd"):
      print "\t-del *." + ext
    print
    if (target in ("mingw32", "vc60")):
      print "clean: clean_win"
    else:
      print "clean: clean_unix"
    print

  def write(self, toolbox, target):
    self.head(target)
    if (hasattr(self, "libraries")):
      for name in self.libraries.keys():
        self.make_library(target, name,
                          self.libraries[name])
    if (hasattr(self, "executables")):
      for name in self.executables.keys():
        self.make_executable(target, name,
                             self.executables[name])
    if (hasattr(self, "examples")):
      for name in self.examples.keys():
        self.make_executable(target, name,
                             self.examples[name],
                             "$(CCTBXLIB) $(LDMATH)")
    if (hasattr(self, "boost_python_modules")):
      for name in self.boost_python_modules.keys():
        self.make_boost_python_module(target, name,
                                      self.boost_python_modules[name])
    self.file_management(target)
    self.tail(toolbox, target)

  def all_targets(self, toolbox):
    for target in ("tru64_cxx", "linux_gcc", "irix_CC",
                   "mingw32", "vc60"):
      makefile = "Makefile." + target
      print "Creating", toolbox, makefile
      f = open(makefile, "w")
      sys.stdout = f
      self.write(toolbox, target)
      sys.stdout = sys.__stdout__
      f.close()

def format_objects(target, objects):
  if (target == "vc60"):
    doto = ".obj"
  else:
    doto = ".o"
  s = ""
  for obj in objects:
    s = s + " " + obj + doto
  return s[1:]
