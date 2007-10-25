++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                              gltbx - GL Toolbox
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Python bindings for OpenGL
--------------------------

The *GL Toolbox* (gltbx) provides Python bindings for OpenGL ``GL_*``
and ``GLU_*`` defines, and ``gl*()`` and ``glu*()`` functions. The
original C interfaces are preserved as much as possible. Except
for the functions that may appear between calls of ``glBegin()``
and ``glEnd()``, OpenGL errors are converted to Python exceptions
(``RuntimeError``). All functions support Python keyword arguments.

The gltbx is part of the cctbx project hosts at:
http://cctbx.sourceforge.net/

See also:

  - http://www.opengl.org/
  - http://pyopengl.sourceforge.net/

Requirements
------------

In addition to the gltbx sources, the following are required:

  - OpenGL headers and runtime libraries in standard places
    (e.g. ``/usr/include`` and ``/usr/lib``)
  - Python 2.2 or higher (http://www.python.org/)
  - Boost source code, version 1.32.0 or higher (http://www.boost.org/)

Note that only the Boost sources are required. It is *not* necessary
to run the Boost.Build system.

Platforms
---------

The gltbx extensions are known to compile on these platforms:

  - Windows 2000/XP, Visual C++ 7.1 (Intel Xeon)
  - Redhat Linux 7.3, 8.0, 9.0, Enterprise Workstation 3 (Intel Xeon)
  - Fedora Core 3 (AMD Opteron)
  - Mac OS X 10.3, Xcode 1.5 (Apple G5)
  - Tru64 Unix, cxx 6.5-042 (HP Alpha)
  - IRIX 6.5, MIPSpro 7.3.1.3 (SGI Mips)

License
-------

BSD-style open source. For details see LICENSE.txt and COPYRIGHT.txt.

Implementation
--------------

The primary OpenGL man pages were copied via wget from
http://www.opengl.org/ (see ``wget_opengl_specs.csh``). The Python
script ``extract_opengl_specs.py`` was used to automatically generate
the file ``opengl_specs.txt``. This derived file is in the CVS and in
the distribution.

Another set of Python scripts is used to generate Boost.Python wrappers
based on the information in ``opengl_specs.txt``; this step takes
only fractions of a second. The Boost.Python wrappers are compiled
to Python extension modules using a C++ compiler. Depending on the
platform and the number of available CPUs this takes between 1 and
10 minutes and results in two extension modules, ``gltbx_gl_ext.so``
and ``gltbx_glu_ext.so``. These modules are imported via two small
Python scripts, ``gl.py`` and ``glu.py``.
For example::

  from gltbx.gl import *
  from gltbx.glu import *
  error = glGetError()
  print gluErrorString(error=error)

OpenGL functions with pointer arguments are handled like this (see
``wx_scene.py``)::

    model = [0]*16
    glGetDoublev(GL_MODELVIEW_MATRIX, model)

Some functions expect strings instead of Python lists, e.g.::

      glTexImage2D(
        GL_TEXTURE_2D, 0, GL_RGBA, self.stripeImageWidth,1,0,
        GL_RGBA, GL_UNSIGNED_BYTE, [self.stripeImage])

Non-const pointers require the Python string to be passed via a
list with one element, as above. On return the original string is
replaced with a new instance. Const pointers also accept Python
strings directly.

Currently the handling of improper pointer arguments is not safe
for all function. See the TODO list below. However, if the correct
arguments are passed everything should work as expected.

Examples
--------

Two simple examples of using the gltbx with wxPython
(http://www.wxpython.org/) are included: ``wx_scene.py`` and ``wx_cube.py``.
These are known to work under Windows 2000/XP using
``wxPython2.6-win32-unicode-2.6.1.0-py24.exe``.

The gltbx extensions are also expected to work with wxPython under
Unix but this has not been tested. It should also be possible to use
the gltbx extensions in combination with window systems other than
wxPython, but this has also not been tested.

gltbx source code
-----------------

The gltbx source code is available via CVS or as complete bundles:

  - CVS: http://cctbx.cvs.sourceforge.net/cctbx/gltbx/

  - Bundles: http://cci.lbl.gov/gltbx/

Installation
------------

Unpack the most recent ``gltbx_####_##_##_####.tar.gz`` in an empty
directory (``####_##_##_####`` is a time stamp). This creates the
subdirectories::

  scons
  libtbx
  boost_adaptbx
  gltbx

For information about SCons see: http://www.scons.org/ .

Create a build directory::

  mkdir build
  cd build

Initialize the build directory::

  python ../libtbx/configure.py boost=../boost_1_32_0 gltbx

Replace ``../boost_1_32_0`` with the proper path name. If the Boost
sources are not installed already, visit the Boost Download page
and unpack the sources in the same directory as the gltbx sources.

The ``libtbx/configure.py`` command above creates the subdirectories::

  build/bin # sh dispatchers for commands to be used
  build/lib # target directory for shared libraries and Python extensions
  build/gltbx # automatically generated Boost.Python wrappers for OpenGL

The ``libtbx/configure.py`` command also creates several ``setpaths*``
files. Source the correct file for your shell, e.g.::

  source setpaths.csh

All this does is add ``build/bin`` to ``$PATH``, thereby enabling
the command for compilation::

  libtbx.scons

This works similar to the good old ``make`` command, but uses SCons
instead. SCons supports parallel compilation. E.g. if 4 CPUs are
available::

  libtbx.scons -j 4

If all goes well, the compilation finishes with lines like these::

  scons: done building targets.
  usr+sys time: 1.13 seconds
  wall clock time: 2 minutes 39.38 seconds (159.38 seconds total)

At this point the installation is complete. If wxPython is installed
try::

  python ../gltbx/wx_scene.py

To use the gltbx installation in a new shell, source the appropriate
``setpaths`` script.

The Windows build procedure is very similar. Simply replace slashes by
backslashes and enter ``setpaths`` instead of ``source setpaths.csh``.

See also: http://cctbx.sourceforge.net/current/installation.html

Pick and choose
---------------

The gltbx extension modules in ``build/lib`` are relocatable, i.e. can
be copied to other places outside the gltbx installation. However, the
Boost.Python runtime library (e.g. ``build/lib/libboost_python.so``)
has to be copied along with the extensions.

The gltbx source code including the auto-generated files in
``build/gltbx`` can also be copied to places outside the gltbx. All
C++ files are platform independent and could be used in other build
systems. It would also be easy to make the source code generators
independent of the ``libtbx`` module if there is an interest.

Prioritized TODO list (help appreciated!)
-----------------------------------------

- All ``gl*()`` and ``glu*()`` functions are wrapped, but some
  functions throw "not implemented" exceptions (``cd build/gltbx;
  grep 'not implemented' *.cpp``). Special code has to be added to
  replace these exceptions. The special code should be included
  in the Python script generating the function wrappers,
  ``generate_functions_bpl.py`` (look for ``special_wrappers``
  and ``special_wrapper_bodies`` to see some examples).

- Many functions taking pointers as arguments do not verify
  that the passed data arrays have the correct size. If data
  of the wrong size are passed, this can result in segmentation
  faults or bus errors. Special code to check expected sizes
  of pointees has to be added, throwing exceptions that can
  be handled in Python.

- Currently three ``gl*()`` functions are disabled because they
  are not available on all platforms listed above. Proper
  defines have to be added to enable these functions on
  platforms where they are available.

Motivation
----------

The main motivation for us was to simplify our build process, and to
minimize additional dependencies. Since we distribute our software to
non-programmers as source code it is important to us that the build
procedure is fully automatic and self-contained. We are already heavily
using Boost.Python and wxPython, and our cross-platform ``libtbx``
build system enables us to add extensions without intruding into Python
installations. In contrast to PyOpenGL, with the gltbx we don't have to
worry about Numeric and GLUT installations and compatibilities.

Contact
-------

Please send suggestions and problem reports to: cctbx@cci.lbl.gov
