from __future__ import absolute_import, division, print_function
import sys

def create_script(
      bundle,
      top_modules,
      test_py="`libtbx.show_dist_paths boost_adaptbx`/tst_rational.py",
      minimum_python_version="2.3"):
  return """\
#! /bin/csh -f

set install_root="$cwd"
set bundle="%(bundle)s"
set sources="$cwd/${bundle}_sources"
set build="$cwd/${bundle}_build"

set build_mode=release

set minimum_python_version=%(minimum_python_version)s

unsetenv PYTHONHOME
unsetenv PYTHONPATH
unsetenv PYTHONSTARTUP
unsetenv PYTHONDEBUG
unsetenv PYTHONINSPECT
unsetenv PYTHONSUPPRESS
unsetenv PYTHONUNBUFFERED
unsetenv PYTHONVERBOSE
unsetenv PYTHONCASEOK
set noglob
unalias cat
unalias cut
unalias cd
unalias grep
unalias head
unalias ls
unalias mkdir
unalias sort
unalias tr
unalias uname

if ("`uname`" == "Darwin") then
  limit stacksize 8192k
endif

if (-f "$sources/TAG") then
  echo "Build tag:"
  cat "$sources/TAG"
endif

if (-d "$sources/boost") then
  set have_sources=1
else
  set have_sources=0
endif

set required_python_version=None
if (-e "$build/lib/PYTHON_VERSION_MAJOR_MINOR") then
  set required_python_version=`grep -v '^#' "$build/lib/PYTHON_VERSION_MAJOR_MINOR"`
endif

set aborted="Installation aborted."

set number_of_cpus=None
set python_exe=None

if ($#argv <= 2) goto process_args

show_usage:
  echo "usage: $0 [/path/to/python/executable] [number_of_cpus]"
  exit 1

def_try_python_exe:
  set python_version_full=None
  set python_version=None
  set ok_minimum_python_version=0
  set q=`("$trial_python_exe" -c 'import sys; from string import upper, split, find, join; print upper("ok"); v=split(sys.version)[0]; print v; ac,ic=split(v, ".")[:2]; print ac+"."+ic; am,im=split("'"$minimum_python_version"'", ".")[:2]; print int(int(ac) < int(am)); print int(int(ac) == int(am) and int(ic) < int(im))' |& cat)`
  if ($#q == 5) then
    if (X"$q[1]" == X"OK") then
      set python_version_full="$q[2]"
      set python_version="$q[3]"
      if (X"$q[4]" == X"0" && X"$q[5]" == X"0") then
        set ok_minimum_python_version=1
      endif
    endif
  endif
  unset q
  goto "$try_python_rtnpt"

def_cmdln_number_of_cpus:
  if (X"$arg" == X) goto show_usage
  set n=(`echo "$arg" | tr -s '[0-9]' '[ *]'`)
  if ($#n == 0) then
    set n=(`echo "$arg"`)
    if ($#n == 1) then
      set number_of_cpus="$n[1]"
      echo "Number of available CPUs: $number_of_cpus"
      if ("$number_of_cpus" == "0") exit 0
    endif
  endif
  goto "$cmdln_number_of_cpus_rtnpt"

def_cmdln_python:
  set trial_python_exe="$arg"
  set try_python_rtnpt=try_cmdln_python_return
  goto def_try_python_exe
  try_cmdln_python_return:
  if ("$python_version" == "None") then
    echo "$0"": command line argument is not a working Python executable: $arg"
    exit 1
  endif
  if (! $ok_minimum_python_version) then
    echo "$0"": command line argument is not a sufficiently recent Python: $arg"
    exit 1
  endif
  set python_exe="$trial_python_exe"
  if (   "$required_python_version" != "None" \
      && "$required_python_version" != "$python_version") then
    echo "$0"": command line argument $python_exe (Python version" \
         "$python_version_full) is not the required version" \
         "($required_python_version)."
    exit 1
  endif
  goto have_python_exe

def_prompt_number_of_cpus:
  if ("$number_of_cpus" == "None") then
    set counter=0
    show_prompt:
    @ counter++
    echo -n "Please enter the number of available CPUs [1]: "
    set arg=(`echo "$<"`)
    if ($#arg == 0) set arg="1"
    set cmdln_number_of_cpus_rtnpt=prompt_n_return
    goto def_cmdln_number_of_cpus
    prompt_n_return:
    if ("$number_of_cpus" == "None") then
      echo "Not a number!"
      if ($counter > 2) then
        echo "Giving up."
        echo "$aborted"
        exit 1
      endif
      echo "Please try again."
      goto show_prompt
    endif
  endif
  goto "$prompt_number_of_cpus_rtnpt"

process_args:

if ($#argv == 1) then
  set arg="$1"
  set cmdln_number_of_cpus_rtnpt=argv1_return
  goto def_cmdln_number_of_cpus
  argv1_return:
  if ("$number_of_cpus" == "None") then
    goto def_cmdln_python
  endif
else if ($#argv == 2) then
  set arg="$1"
  set cmdln_number_of_cpus_rtnpt=argv2_return1
  goto def_cmdln_number_of_cpus
  argv2_return1:
  if ("$number_of_cpus" == "None") then
    set arg="$2"
    set cmdln_number_of_cpus_rtnpt=argv2_return2
    goto def_cmdln_number_of_cpus
    argv2_return2:
    if ("$number_of_cpus" == "None") goto show_usage
    set arg="$1"
    goto def_cmdln_python
  else
    set arg="$2"
    goto def_cmdln_python
  endif
endif

set python_sources=None
if ("$required_python_version" == "None") then
  set python_sources=(`ls | grep '^Python-'`)
  if ($#python_sources == 0) then
    set python_sources=None
  else if ($#python_sources == 1) then
    set python_sources="$python_sources[1]"
  else
    echo "ERROR: Multiple Python source code directories:"
    set noglob
    foreach d ($python_sources)
      echo "         $d"
    end
    echo "       Move or remove all but one directory."
    exit 1
  endif
endif

if ("$python_sources" != "None") then
  set prompt_number_of_cpus_rtnpt=inst_py_src_return
  goto def_prompt_number_of_cpus
  inst_py_src_return:
  #
  echo "Installing $python_sources from sources"
  mkdir -p "$build"
  cd "$build"
  cd ..
  cd "$python_sources"
  set py_install_log="../py_install_log"
  echo "Configuring Python"
  if ("X`uname`" == "XHP-UX") then
    # tested with HP aC++/ANSI C B3910B A.06.06 [Nov 7 2005]
    env CC=cc CXX="aCC -mt" BASECFLAGS="+DD64 -mt" LDFLAGS="+DD64 -lxnet" ./configure --without-gcc --prefix="$build/base" >& "$py_install_log"
    if (-f pyconfig.h) then
      grep -v _POSIX_C_SOURCE pyconfig.h > zz; mv zz pyconfig.h
    endif
  else
    ./configure --prefix="$build/base" >& "$py_install_log"
  endif
  echo "Compiling Python. This may take a while."
  if ("$number_of_cpus" == "1" || X"`uname`" != X"Linux") then
    echo "Command: make" >>& "$py_install_log"
    make >>& "$py_install_log"
  else
    echo "Command: make -j $number_of_cpus" >>& "$py_install_log"
    make -j "$number_of_cpus" >>& "$py_install_log"
  endif
  echo "Installing Python"
  make install >>& "$py_install_log"
  echo "Done installing Python."
  cd "$install_root"
  set trial_python_exe="$build/base/bin/python"
  set try_python_rtnpt=try_base_bin_python_return
  goto def_try_python_exe
  try_base_bin_python_return:
  if ("$python_version" == "None") then
    echo "ERROR: Python installation failed."
    echo "       Please check the log file for errors:"
    echo "         $py_install_log"
    exit 1
  endif
  set python_exe="$trial_python_exe"
  goto have_python_exe
endif

set trial_python_exe="$build/base/bin/python"
if (-e "$trial_python_exe") then
  set try_python_rtnpt=try_python_bb_return
  goto def_try_python_exe
  try_python_bb_return:
  if ("$python_version" == "None") then
    echo "Sorry: The installer is not suitable for this platform."
    echo "$aborted"
    exit 1
  endif
  set python_exe="$trial_python_exe"
  goto have_python_exe
endif
#
set trial_python_exe="python"
set try_python_rtnpt=try_python_return
goto def_try_python_exe
try_python_return:
if ("$python_version" != "None") then
  if (   "$required_python_version" == "None" \
      || "$required_python_version" == "$python_version") then
    if ($ok_minimum_python_version) then
      set python_exe="$trial_python_exe"
      goto have_python_exe
    endif
    echo "INFO: python on PATH (version $python_version_full) is not a" \
         "recent enough version (minimum $minimum_python_version)."
  else
    echo "INFO: python on PATH (version $python_version_full) is not the" \
         "required version ($required_python_version)."
  endif
endif
#
if ("`uname`" == "Darwin") then

  set lib_fw_vers="/Library/Frameworks/Python.framework/Versions"
  set sys_lib_fw_vers="/System/Library/Frameworks/Python.framework/Versions"

  set pyv="$required_python_version"
  if ("$pyv" == "None") then
    set pyv="`ls "$lib_fw_vers" "$sys_lib_fw_vers" | grep '^[1-9]' | sort -r | head -1`"
    if ("X$pyv" == "X") set pyv=None
  endif
  if (-d "$lib_fw_vers/$pyv") then
    set python_exe="$lib_fw_vers/$pyv/bin/python"
    if (-x "$python_exe") goto have_python_exe
    if (-f "$python_exe") then
      echo "INFO: $python_exe is not executable."
    else
      echo "INFO: $python_exe does not exist."
    endif
  endif
  if (-d "$sys_lib_fw_vers/$pyv") then
    set python_exe="$sys_lib_fw_vers/$pyv/bin/python"
    if (-x "$python_exe") goto have_python_exe
    if (-f "$python_exe") then
      echo "INFO: $python_exe is not executable."
    else
      echo "INFO: $python_exe does not exist."
    endif
  endif
#
endif
#
if ("$required_python_version" != "None") then
  set trial_python_exe="/usr/local/bin/python$required_python_version"
  set try_python_rtnpt=try_ulb_python_return
  goto def_try_python_exe
  try_ulb_python_return:
  if ("$python_version" != "None") then
    set python_exe="$trial_python_exe"
    goto have_python_exe
  endif
endif
#
if ("$required_python_version" == "None") then
  set trial_python_exe="/usr/bin/python"
  set try_python_rtnpt=try_ub_python_return
  goto def_try_python_exe
  try_ub_python_return:
  if ("$python_version" != "None" && $ok_minimum_python_version) then
    set python_exe="$trial_python_exe"
    goto have_python_exe
  endif
endif
#
if ("$required_python_version" != "None") then
  set trial_python_exe="/usr/bin/python$required_python_version"
  set try_python_rtnpt=try_ub_pythonv_return
  goto def_try_python_exe
  try_ub_pythonv_return:
  if ("$python_version" != "None") then
    set python_exe="$trial_python_exe"
    goto have_python_exe
  endif
endif

if ("$required_python_version" == "None") then
  set pyv=""
else
  set pyv=" (version $required_python_version)"
endif
echo ""
echo "Sorry: Cannot find a suitable Python interpreter$pyv."
echo ""
echo "  Please download an installer with Python included or"
echo "  run this command with the full path to a suitable Python"
echo "  executable as an argument, e.g.:"
echo ""
if ("$number_of_cpus" == "None") then
  set n=""
else
  set n=" $number_of_cpus"
endif
echo "    $0 some/path/bin/python$n"
echo ""
echo "$aborted"
exit 1

have_python_exe:

if ($have_sources) then
  set prompt_number_of_cpus_rtnpt=have_py_exe_return
  goto def_prompt_number_of_cpus
  have_py_exe_return:
endif

"$python_exe" -V

echo ""
echo "Precompiling all .py files. This may take a minute or two."
"$python_exe" "$sources/libtbx/command_line/py_compile_all.py"

echo ""
if (! -d "$build") mkdir -p "$build"
cd "$build"
echo "Configuring $bundle build directory"
"$python_exe" "$sources/libtbx/configure.py" --current_working_directory="$build" --build="$build_mode" %(top_modules)s
source setpaths.csh

if ($have_sources) then
  echo ""
  echo "Installing $bundle modules. This may take a while."
  libtbx.scons -j "$number_of_cpus" .
endif

set test_py="%(test_py)s"
if (-f "$test_py") then
  echo ""
  echo "Running a selected test"
  set cmd='libtbx.python "'"$test_py"'"'
  echo "$cmd"
  eval $cmd
endif

cat << EOT

***
*** csh and tcsh users:
***     To use this installation in a new shell or process run the command:
***
***         source "$build/setpaths.csh"
***
***     You may want to add this line to your .cshrc file.
***
*** sh and bash users:
***     To use this installation in a new shell or process run the command:
***
***         . "$build/setpaths.sh"
***
***     You may want to add this line to your .profile or .bashrc file.
***
EOT
""" % vars()

if (__name__ == "__main__"):
  assert len(sys.argv) == 3
  sys.stdout.write(create_script(sys.argv[1], sys.argv[2]))
