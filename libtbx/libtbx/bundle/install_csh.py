import sys

def create_script(bundle, top_modules):
  return """\
#! /bin/csh -f

set install_root=$cwd
set bundle=%(bundle)s
set sources=$cwd/${bundle}_sources
set build=$cwd/${bundle}_build

unalias cat
unalias cd
unalias grep
unalias ls
unalias mkdir

set pythonhome_was_defined=0
if ($?PYTHONHOME) then
  set pythonhome_was_defined=1
  unsetenv PYTHONHOME
endif

if (-f "$sources/TAG") then
  echo "Build tag:"
  cat "$sources/TAG"
endif

if (-d $sources/boost) then
  set have_sources=1
else
  set have_sources=0
endif

set python_exe=None
set build_mode=release

if ("`uname`" == "Darwin") then
  set python_exe=/Library/Frameworks/Python.framework/Versions/2.3/bin/python
  if (! -x $python_exe) then
    set python_exe=/System$python_exe
  endif
  $python_exe -V
  if ($status != 0) then
    echo "Under Mac OS 10 Python 2.3 must be pre-installed."
    echo "Please refer to the following web page for more information:"
    echo "http://cci.lbl.gov/cctbx_build/mac_os_x_notes.html"
    exit 1
  endif
endif

if ($have_sources == 0) then

  if (-d $build/python) then
    set python_exe=$build/python/bin/python
  else if ($python_exe == None) then
    set python_exe=python
  endif

else

  if ($#argv == 0) then
    echo -n "Please enter the number of available CPU's [1]: "
    set n_cpu_s=(`echo "$<"`)
    if ($#n_cpu_s > 1) then
      echo "Not a number! Please try again."
      exit 1
    else if ($#n_cpu_s == 0) then
      set n_cpu_s=1
    else
      set n_cpu_s="$n_cpu_s[1]"
    endif
    if ($n_cpu_s == "0") exit 0
  else if ($#argv == 1) then
    set n_cpu_s=$1
  else
    echo "usage: install_from_scratch.csh number_of_cpu_s"
    exit 1
  endif
  echo "Number of available CPU's: $n_cpu_s"

  if ($python_exe == None) then

    set python_sources=(`ls | grep Python-`)
    if ($#python_sources == 0) then
      set python_sources=None
    else if ($#python_sources == 1) then
      set python_sources=$python_sources[1]
    else
      echo "ERROR: Multiple Python source code directories."
      echo "       Move or remove all but one directory."
      exit 1
    endif

    if ($python_sources == None) then

      echo "Trying to find a pre-installed Python:"
      set python_exe=None
      if (-x $build/python/bin/python) then
        $build/python/bin/python -V |& head -1
        if ($status == 0) then
          set python_exe=$build/python/bin/python
        endif
      endif
      if ($python_exe == None) then
        python -V |& head -1
        if ($status == 0) then
          set python_exe=python
        endif
      endif
      if ($python_exe == None) then
        python2 -V |& head -1
        if ($status == 0) then
          set python_exe=python2
        endif
      endif
      if ($python_exe == None) then
        echo "Cannot find a working Python."
      else
        set python_version=(`$python_exe -V |& tr "." " "`)
        if ($python_version[2]) then
          set minor=`echo $python_version[3] | cut -c-1`
          if ($minor < 2 || ($minor == 2 && $#python_version == 3)) then
            echo "A more recent Python version is required (2.2.1 or higher)."
            set python_exe=None
          endif
        endif
      endif
      if ($python_exe == None) then
        echo "Please use the installer that includes the Python sources."
        exit 1
      endif

    else

      echo "Installing $python_sources from sources"
      mkdir -p $build
      cd $build
      cd ..
      cd $python_sources
      set py_install_log=../py_install_log
      echo "Configuring Python"
      ./configure --prefix="$build/python" >& $py_install_log
      echo "Compiling Python. This may take a while."
      make >>& $py_install_log
      echo "Installing Python"
      make install >>& $py_install_log
      echo "Done installing Python."
      cd $install_root
      set python_exe=$build/python/bin/python
      $python_exe -V
      if ($status != 0) then
        echo "ERROR: Python installation failed."
        echo "Please check the log file for errors:"
        echo "  $py_install_log"
        exit 1
      endif

    endif

  endif

  mkdir -p $build

endif

echo ""
echo "Precompiling all .py files. This may take a minute or two."
$python_exe $sources/libtbx/libtbx/command_line/py_compile_all.py

echo ""
cd $build
echo "Configuring $bundle build directory"
$python_exe $sources/libtbx/configure.py --build=$build_mode %(top_modules)s
source setpaths.csh

if ($have_sources != 0) then
  echo ""
  echo "Installing $bundle modules. This may take a while."
  libtbx.scons -j $n_cpu_s .
endif

set test_py="$SCITBX_DIST/lbfgs/boost_python/tst_lbfgs.py"
if (-f "$test_py") then
  echo ""
  echo "Running a selected test"
  set cmd="python $test_py"
  echo "$cmd"
  eval $cmd
endif

if ($pythonhome_was_defined != 0) then
  echo ""
  echo "####################"
  echo "     \!\!\!\!\!\!\!\!\!       "
  echo "     ATTENTION        "
  echo "     \!\!\!\!\!\!\!\!\!       "
  echo "####################"
  echo "PYTHONHOME conflict\!"
  echo "####################"
  echo "This $bundle installation will work only if the PYTHONHOME"
  echo "environment variable is *not* defined. Therefore the"
  echo "PYTHONHOME variable will be undefined each time this"
  echo "$bundle installation is initialized. Software packages"
  echo "that require PYTHONHOME (such as the MGL Tools) cannot be"
  echo "used in the same shell. Please open a separate shell for"
  echo "each application."
endif

cat << EOT

***
*** To use this installation in a new shell or process run the command:
***
***     source $LIBTBX_BUILD/setpaths.csh
***
*** You may want to add this line to your .cshrc file.
***
EOT
""" % vars()

if (__name__ == "__main__"):
  assert len(sys.argv) == 3
  sys.stdout.write(create_script(sys.argv[1], sys.argv[2]))
