#!/usr/bin/python

from __future__ import absolute_import, division, print_function

import os.path
import sys
from optparse import OptionParser

def run(args, prologue=None, epilogue=None, out=sys.stdout):
  parser = OptionParser(
    description="Generate the dispatcher include file for using "+
      "locally installed GUI (and related) libraries.  (Used in Phenix "+
      "installer and binary cctbx_plus builds.)")
  parser.add_option("--suffix", dest="suffix", action="store",
    help="Suffix for dispatcher_include file", default="gui")
  parser.add_option("--build_dir", dest="build_dir", action="store",
    help="Build directory (and dispatcher destination)", default=os.getcwd())
  parser.add_option("--prologue", dest="prologue", action="store",
    help="Additional dispatcher include to prepend")
  parser.add_option("--epilogue", dest="epilogue", action="store",
    help="Additional dispatcher include to append")
  parser.add_option("--base_dir", dest="base_dir", action="store",
    help="Directory for base packages (Python etc.)", default=None)
  parser.add_option("--ignore_missing_dirs", dest="ignore_missing_dirs",
    action="store_true", help="Don't raise error if GTK paths don't exist")
  parser.add_option("--quiet", action="store_true")
  parser.add_option("--use_conda", dest="use_conda", action="store_true",
    help="Use conda dependencies")
  # this is detached in case we need to upgrade GTK - setting GTK_PATH is
  # essential for the themes to be used correctly
  parser.add_option("--gtk_version", dest="gtk_version", action="store",
    help="Version number (major.minor.rev) for GTK+", default="2.10.0")
  options, args = parser.parse_args(args)
  build_path = options.build_dir
  if (not os.path.isdir(build_path)):
    raise OSError("The specified build directory (%s) does not exist." %
      build_path)
  build_path = os.path.abspath(build_path)
  base_path = options.base_dir
  if (base_path is None):
    base_path = os.path.join(build_path, "base")
  if (not os.path.isdir(base_path)):
    raise OSError("%s does not exist." % base_path)
  base_path = os.path.abspath(base_path)
  ld_library_paths = [ os.path.join(base_path, "lib"), ]
  if (sys.platform.startswith('linux')):
    lib64_path = os.path.join(base_path, "lib64")
    if (os.path.isdir(lib64_path)):
      ld_library_paths.append(lib64_path)
  dyld_library_paths = ld_library_paths
  if not options.use_conda and sys.platform == "darwin":
    dyld_library_paths += [
    os.path.join(base_path,"Python.framework","Versions","Current","lib"), ]
  check_libs = ld_library_paths
  if (sys.platform == "darwin") : check_libs = dyld_library_paths
  if sys.platform != "win32":
    for lib_dir in check_libs :
      if (not os.path.isdir(lib_dir)):
        raise OSError("%s does not exist." % lib_dir)
  gtk_path = os.path.join(base_path, "lib", "gtk-2.0", options.gtk_version)
  if ((sys.platform.startswith("linux")) and
      (not os.path.isdir(gtk_path)) and
      (not options.ignore_missing_dirs)):
    raise OSError("The path for the specified version of GTK+ does not "+
      "exist (%s)." % gtk_path)
  dispatcher = os.path.join(build_path, "dispatcher_include_%s.sh" %
    options.suffix)
  if sys.platform == "win32":
    dispatcher = os.path.join(build_path, "dispatcher_include_%s.bat" %
      options.suffix)
  f = open(dispatcher, "w")
  if (prologue is not None):
    f.write(prologue + "\n")
  if (options.prologue is not None):
    if not os.path.exists(options.prologue) and options.prologue.find("\n")>-1:
      f.write("%s\n" % options.prologue)
    else:
      f.write(open(options.prologue).read() + "\n")
  if sys.platform != "win32" and not options.use_conda:
    print("""\
# include at start
if [ "$LIBTBX_DISPATCHER_NAME" != "libtbx.scons" ] && \
   [ -z "$PHENIX_TRUST_OTHER_ENV" ]; then
  # work around broken library environments
  LD_LIBRARY_PATH=""
  PYTHONPATH=""
fi
# include before command
#
# GUI dependencies
CCTBX_BUILD_BASE="%s"
LIBTBX_OS_NAME=`uname -s`
if [ "$PHENIX_GUI_ENVIRONMENT" = "1" ]; then
  if [ -z "$DISABLE_PHENIX_GUI" ]; then
    export BOOST_ADAPTBX_FPE_DEFAULT=1
    export BOOST_ADAPTBX_SIGNALS_DEFAULT=1
  fi
  # echo $LIBTBX_OS_NAME
  if [ "$LIBTBX_OS_NAME" = "Linux" ]; then
    export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
    export OLD_XDG_DATA_DIRS=$XDG_DATA_DIRS
    if [ -z "$LD_LIBRARY_PATH" ]; then
      LD_LIBRARY_PATH=%s
    else
      LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH
    fi
    export LD_LIBRARY_PATH
#
    unset GTK_MODULES
    GTK_PATH=$CCTBX_BUILD_BASE/lib/gtk-2.0/%s
    PANGO_RC_FILE=$CCTBX_BUILD_BASE/etc/pango/pangorc
    GTK_IM_MODULE_FILE=$CCTBX_BUILD_BASE/etc/gtk-2.0/gtk.immodules
    GDK_PIXBUF_MODULE_FILE=$CCTBX_BUILD_BASE/etc/gtk-2.0/gdk-pixbuf.loaders
    GTK2_RC_FILES=$CCTBX_BUILD_BASE/share/themes/Clearlooks/gtk-2.0/gtkrc
    FONTCONFIG_PATH=$CCTBX_BUILD_BASE/etc/fonts
    FONTCONFIG_FILE=$CCTBX_BUILD_BASE/etc/fonts/fonts.conf
    if [ -z "$XDG_DATA_DIRS" ]; then
      XDG_DATA_DIRS=$CCTBX_BUILD_BASE/share:/usr/share
    else
      XDG_DATA_DIRS=$CCTBX_BUILD_BASE/share:$XDG_DATA_DIRS
    fi
    GNOME_DISABLE_CRASH_DIALOG=1
    export PANGO_RC_FILE
    export GTK_IM_MODULE_FILE
    export GDK_PIXBUF_MODULE_FILE
    export FONTCONFIG_PATH
    export FONTCONFIG_FILE
    export XDG_DATA_DIRS
    export GTK_PATH
    export GTK2_RC_FILES
    export GNOME_DISABLE_CRASH_DIALOG
    PATH=$CCTBX_BUILD_BASE/bin:$PATH
    export PATH
  fi
fi
""" % (base_path, ":".join(ld_library_paths), ":".join(ld_library_paths), options.gtk_version), file=f)
  # restore some variables for conda
  if sys.platform != "win32" and options.use_conda:
    print("""
# include at start
if [ "$LIBTBX_DISPATCHER_NAME" != "libtbx.scons" ] && \
   [ -z "$PHENIX_TRUST_OTHER_ENV" ]; then
  # work around broken library environments
  LD_LIBRARY_PATH=""
  DYLD_LIBRARY_PATH=""
  DYLD_FALLBACK_LIBRARY_PATH=""
  PYTHONPATH=""
fi
# include before command
if [ "$PHENIX_GUI_ENVIRONMENT" = "1" ]; then
  if [ -z "$DISABLE_PHENIX_GUI" ]; then
    export BOOST_ADAPTBX_FPE_DEFAULT=1
    export BOOST_ADAPTBX_SIGNALS_DEFAULT=1
  fi
fi
""", file=f)
  # QBio DivCon paths
  if sys.platform != "win32":
    print("""
if [ ! -z "$QB_PYTHONPATH" ]; then
  export PYTHONPATH=$PYTHONPATH:$QB_PYTHONPATH
fi
if [ ! -z "$QB_LD_LIBRARY_PATH" ]; then
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$QB_LD_LIBRARY_PATH
fi
if [ ! -z "$QB_DYLD_LIBRARY_PATH" ]; then
  export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$QB_DYLD_LIBRARY_PATH
fi
""", file=f)
  if (epilogue is not None):
    f.write(epilogue + "\n")
  if (options.epilogue is not None):
    f.write(open(options.epilogue).read() + "\n")
  f.close()
  if (not options.quiet):
    print("Wrote %s" % dispatcher, file=out)
    print("You should now run libtbx.refresh to regenerate dispatchers.", file=out)

# obsolete???
pymol_paths = """
 if [ "$PHENIX_MTYPE" != "mac-ppc-osx" ] && \
    [ "$PHENIX_MTYPE" != "mac-intel-osx" ] && \
    [ "$PHENIX_MTYPE" != "mac-intel-osx-x86_64" ]; then
   export PYMOL_PATH=$PHENIX/pymol
 fi
"""

if (__name__ == "__main__"):
  run(sys.argv[1:])
