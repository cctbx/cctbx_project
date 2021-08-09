#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

import sys
from   libtbx import easy_run
from   os     import listdir



def check_libcpp():

  #______________________________________________________________________________
  # get the libstd version that the python executable links against
  #
  python_libstdcxx_so = None
  ldd_exec = easy_run.fully_buffered(command="ldd " + str(sys.executable))
  for line in ldd_exec.stdout_lines:
    if line.strip().startswith("libstdc++.so"):
      python_libstdcxx_so = line.split()[0]
      break

  if python_libstdcxx_so is None:
    return False
  #
  #------------------------------------------------------------------------------

  #______________________________________________________________________________
  # Go through all installed *.so's and ensure that they are linking to the same
  # libstdc++ as python executable
  #
  lib_dir = join(
      libtbx.env.lib_path._anchor._path, libtbx.env.lib_path.relocatable
  )
  libs = [join(lib_dir, name) for name in listdir(lib_dir) if name.endswith(".so")]
  for lib in libs:
    ldd_lib = easy_run.fully_buffered(command="ldd " + str(lib))
    for line in ldd_lib.stdout_lines:
      if line.strip().startswith("libstdc++.so"):
        lib_libstdcxx_so = line.split()[0]
        if lib_libstdcxx_so != python_libstdcxx_so:
          raise SystemError(
            "FATAL: libstdc++.so mismatch\n"
            "sys.executable: " + str(sys.executable) + "\n"
            "lib_file:" + str(lib)
          )
        break
  #
  #------------------------------------------------------------------------------

  return True



if __name__=="__main__":

  if not sys.platform.startswith("linux"):
    print("This test only works on Linux at the moment")
  else:
    if check_libcpp():
      print("Success!")
    else:
      print(
        "This test was not run because your python executable was not linked "
        "against libstdc++ => this test is not needed."
      )
