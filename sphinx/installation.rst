
.. _installation:
.. contents:: Table of Contents

---------------------
Installation overview
---------------------

Periodically, the most recent cctbx source code including all its
dependencies is automatically exported from the source code repositories
(git) and bundled into compressed binary installers that are published
at:

  - http://cci.lbl.gov/cctbx_build/

This web page provides compressed binary distributions for a
variety of platforms. It is most convenient to use these binary bundles
if possible. Installation is very simple and fast. After uncompressing the
bundle, run the ``install --prefix=<installation location>`` command to
install on macOS and Linux. On Windows, uncompressing the zip file will provide
a working copy of cctbx. There is not installation script.

----------------------------------------------------
Manually building from sources under Linux and macOS
----------------------------------------------------

Please note: **The following instructions are for developers!**

Building from sources requires Python 2.7, 3.6 or newer and a C++
compiler. If you like to use the most recent Python, it can be
installed in the following way::

  tar -xf Python-2.7.13.tar.xz
  cd Python-2.7.13
  ./configure --prefix=/your/choice
  make
  make install

It may be convenient (but is not required) to add the directory
``/your/choice/bin`` to the command-line search ``PATH``, e.g. using
``csh``::

  set path=(/your/choice/bin $path)

Recent cctbx sources are available from the cctbx GitHub_ repository. To
download the repository and build cctbx, use bootstrap.py::

  mkdir <installation directory>
  cd <installation directory>
  wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py
  python ./bootstrap.py

On macOS, since ``wget`` is not available by default, use ``curl`` instead::

  curl https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py > bootstrap.py

After some time, this creates the subdirectories ``base``, ``base_tmp``,
``build``, and ``modules``. The ``base`` directory contains dependencies for
cctbx, ``base_tmp`` is a temporary directory for compiling dependencies (can be
deleted), ``build`` contains the compiled cctbx code, and ``modules`` contains
the source code for cctbx. To keep ``bootstrap.py`` in your
``<installation directory>`` up-to-date with the version in ``modules``, you
can create a symbolic link to that version::

  cd <installation directory>
  rm bootstrap.py
  ln -s ./modules/cctbx_project/libtbx/auto_build/bootstrap.py

Within the ``build`` directory, cctbx creates a
file ``setpaths.csh`` (among others). This file must be used to
initialize a new shell or process with the cctbx settings::

  source setpaths.csh

There is also a ``setpaths.sh`` for ``bash`` users.

To update your cctbx installation to the latest version, you can just run::

  python ./bootstrap.py

again in the ``<installation directory>``. This will update the source code in
the ``modules`` directory and recompile the changes (if necessary) in
``build``. Occasionally, dependencies in ``base`` are updated. When this
happens, just delete the ``base`` and ``base_tmp`` directories and rerun
the ``bootstrap.py`` command.

To compile any local changes to the source code, enter the ``build`` directory
and run::

  make

This will actually run the ``libtbx.scons`` command using all
available CPUs. You can also manually specify the number of CPUs to
use, for example::

  libtbx.scons -j 4

Note that ``libtbx.scons`` is just a thin wrapper around SCons_. The
`SCons documentation`_ applies without modification.

To run scripts with cctbx imports use the command::

  libtbx.python your_script.py

(You can also use ``scitbx.python``, ``cctbx.python``, ``iotbx.python``, etc.;
all these commands are equivalent.)

For example, to run some regression tests after the compilation is
finished enter::

  source build/setpaths_all.csh
  libtbx.python $SCITBX_DIST/run_tests.py
  libtbx.python $CCTBX_DIST/run_tests.py --Quick

The output should show many OK. A Python Traceback is an indicator
for problems.

-----------------------------------------------------------
Manually building from sources under Windows 2000 or higher
-----------------------------------------------------------

The cctbx installation requires Visual C++ 8.0 (Visual Studio .NET
2005) or higher.

To install Python under Windows it is best to use a binary
installer from the `Python download page <http://www.python.org/download/>`_.
The default choices presented by the installation wizard are usually fine.

Recent self-contained cctbx sources are available in the
self-extracting file
`cctbx_bundle.exe <http://cci.lbl.gov/cctbx_build/results/last_published/cctbx_bundle.exe>`_
published at the cctbx build page. To unpack this file in a new, empty
directory enter::

  cctbx_bundle.exe

This creates a subdirectory ``cctbx_sources``. The installation
procedure should be executed in another directory, e.g.::

  mkdir cctbx_build
  cd cctbx_build
  C:\python27\python.exe ..\cctbx_sources\libtbx\configure.py mmtbx

The last command initializes the ``cctbx_build`` directory and creates
a file ``setpaths.bat`` (among others). This file must be used to
initialize a new shell or process with the cctbx settings::

  setpaths.bat

To compile enter::

  libtbx.scons

On a machine with multiple CPUs enter::

  libtbx.scons -j N

where N is the number of CPUs available.

Note that ``libtbx.scons`` is just a thin wrapper around SCons_. The
`SCons documentation`_ applies without modification.

To run scripts with cctbx imports use the command::

  libtbx.python your_script.py

(You can also use ``scitbx.python``, ``cctbx.python``, ``iotbx.python``, etc.;
all these commands are equivalent.)

For example, to run some regression tests after the compilation is
finished enter::

  setpaths_all.bat
  libtbx.python %SCITBX_DIST%\run_tests.py
  libtbx.python %CCTBX_DIST%\run_tests.py --Quick

The output should show many OK. A Python Traceback is an indicator
for problems.

Back_

.. _Back: introduction.html

.. _SCons: http://www.scons.org/
.. _`SCons documentation`: http://www.scons.org/doc/HTML/scons-man.html
.. _Boost: http://www.boost.org/
.. _`boost SVN tree`: http://svn.boost.org/trac/boost/wiki/BoostSubversion
.. _CCP4: http://www.ccp4.ac.uk/
.. _GitHub: https://github.com/cctbx/cctbx_project/
