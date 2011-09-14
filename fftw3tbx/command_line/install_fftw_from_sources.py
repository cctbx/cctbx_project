import fftw3tbx
from libtbx import easy_run
from libtbx.utils import Sorry, Usage
from libtbx.str_utils import show_string
import libtbx.introspection
import libtbx.load_env
import sys, os

def install(tgz, precision, n_proc):
  assert precision in ["double", "float"]
  #
  print "cd %s" % show_string(libtbx.env.build_path)
  os.chdir(libtbx.env.build_path)
  print
  #
  command = "gunzip -c %s | tar xvf -" % show_string(tgz)
  print command
  lines = easy_run.fully_buffered(
    command=command).raise_if_errors().stdout_lines
  print
  if (len(lines) > 10):
    lines = lines[:3] \
          + ["... (%d lines not shown)" % (len(lines)-6)] \
          + lines[-3:]
  print " ", "\n  ".join(lines)
  print
  if (len(lines) == 0):
    raise Sorry(
      "Source code file appears to be empty: %s" % show_string(tgz))
  src = os.path.dirname(lines[0])
  if (not os.path.isdir(src)):
    raise Sorry(
      "Expected source directory does not exist: %s" % show_string(src))
  #
  print "cd %s" % show_string(src)
  os.chdir(src)
  print
  #
  if (precision == "float"):
    s = " --enable-single"
    libfftw3 = fftw3tbx.libfftw3f
  else:
    s = ""
    libfftw3 = fftw3tbx.libfftw3
  command = "./configure --prefix=%s --enable-shared%s" % (show_string(
    libtbx.env.under_build("base")), s)
  print command
  easy_run.call(command=command)
  print
  if (not os.path.isfile("Makefile")):
    raise Sorry(
      "./configure command failed: Makefile does not exist.")
  #
  command = "make"
  if (n_proc is not None): command += " -j%d" % n_proc
  print command
  easy_run.call(command=command)
  print
  f = ".libs/" + libfftw3
  if (not os.path.isfile(f)):
    raise Sorry(
      "make command failed: %s does not exist." % show_string(f))
  #
  command = "make install"
  print command
  easy_run.call(command=command)
  print
  for f in ["base/include/"+fftw3tbx.fftw3_h,
            "base/lib/"+libfftw3]:
    f = libtbx.env.under_build(f)
    if (not os.path.isfile(f)):
      raise Sorry(
        "make install command failed: %s does not exist." % show_string(f))
  #
  print "cd %s" % show_string(libtbx.env.build_path)
  os.chdir(libtbx.env.build_path)
  print
  #
  command = "rm -rf %s" % show_string(src)
  print command
  easy_run.call(command=command)
  print

def run(args):
  print
  if (len(args) != 1):
    raise Usage("%s fftw-*.tar.gz" % libtbx.env.dispatcher_name)
  tgz = args[0]
  if (not os.path.isfile(tgz)):
    raise Sorry(
      "Not a file: %s" % show_string(tgz))
  tgz = os.path.abspath(tgz)
  n_proc = libtbx.introspection.number_of_processors()
  install(tgz=tgz, precision="double", n_proc=n_proc)
  install(tgz=tgz, precision="float", n_proc=n_proc)
  #
  command = "libtbx.refresh"
  print command
  easy_run.call(command=command)
  print
  #
  command = "libtbx.scons"
  if (n_proc is not None): command += " -j%d" % n_proc
  print command
  easy_run.call(command=command)
  print

if (__name__ == "__main__"):
  run(sys.argv[1:])
