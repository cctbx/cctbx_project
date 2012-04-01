import sys

def create_script(bundle, top_modules, single_dir=False):
  py_major, py_minor = sys.version_info[:2]
  script = r"""@echo off

set PYTHONHOME=
set PYTHONSTARTUP=
set PYTHONDEBUG=
set PYTHONINSPECT=
set PYTHONSUPPRESS=
set PYTHONUNBUFFERED=
set PYTHONVERBOSE=
set PYTHONCASEOK=1
"""
  if (single_dir) :
    script += r"""
cd %(bundle)s
"""
  script += r"""
if not exist %(bundle)s_sources\TAG goto find_python
echo.
echo Build tag:
type %(bundle)s_sources\TAG

:find_python
echo.
echo Trying to find Python:
cd %(bundle)s_build
if not exist base\python\python.exe goto try_plain_python
set python=base\python\python
call "%%python%%" -V
if %%errorlevel%% == 0 goto have_python
:try_plain_python
set python=python
call "%%python%%" -V
if %%errorlevel%% == 0 goto have_python
set python=C:\python%(py_major)d%(py_minor)d\python
call "%%python%%" -V
if %%errorlevel%% == 0 goto have_python
cd ..
echo.
echo Cannot find a Python interpreter.
echo.
echo Please download an installer with Python included
echo or add a matching Python to the PATH environment variable.
echo.
echo Installation aborted.
echo.
goto end

:have_python
echo.
echo Configuring %(bundle)s build directory
call "%%python%%" ..\%(bundle)s_sources\libtbx\configure.py %(top_modules)s
set el=%%errorlevel%%
cd ..
if not %%el%% == 0 goto end
call %(bundle)s_build\setpaths_all.bat

if not exist %(bundle)s_install_finalize.bat goto run_tests
echo Running final setup script %(bundle)s_install_finalize.bat
call %(bundle)s_install_finalize.bat

:run_tests
if not exist "%%BOOST_ADAPTBX_DIST%%\tst_rational.py" goto skip_test
echo.
echo Running a selected test
echo libtbx.python "%%BOOST_ADAPTBX_DIST%%\tst_rational.py"
call libtbx.python "%%BOOST_ADAPTBX_DIST%%\tst_rational.py"
if not %%errorlevel%% == 0 goto end
:skip_test

echo.
echo ***
echo *** To use this installation in a new shell or process run the command:
echo ***
echo ***     %%LIBTBX_BUILD%%\setpaths.bat
echo ***
echo *** You may want to add this line to some startup file.
echo ***

:end
if not defined LIBTBX_BATCH_INSTALL goto end_prompt
if not %%LIBTBX_BATCH_INSTALL%% == 0 goto final_exit
:end_prompt
pause
:final_exit
"""
  return script % vars()

if (__name__ == "__main__"):
  assert len(sys.argv) == 3
  sys.stdout.write(create_script(sys.argv[1], sys.argv[2]))
