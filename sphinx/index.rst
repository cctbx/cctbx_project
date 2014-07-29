
+++++++++++++++++++++++++++++++++++++
Computational Crystallography Toolbox
+++++++++++++++++++++++++++++++++++++

.. _introduction:
.. toctree::
   :maxdepth: 1

   installation
   tour
   build_system
   history
   libtbx/index
   boost_adaptbx/index
   iotbx/index
   scitbx/index
   cctbx/cctbx
   mmtbx/index


Welcome to CCTBX's documentation!
=================================

The Computational Crystallography Toolbox (cctbx) is being developed as the open
source component of the PHENIX_ system. The goal of the PHENIX project is to
advance automation of macromolecular structure determination. PHENIX depends on
the cctbx, but not vice versa. This hierarchical approach enforces a clean
design as a reusable library. The cctbx is therefore also useful for
small-molecule crystallography and even general scientific applications.

To maximize reusability and, maybe even more importantly, to give individual
developers a notion of privacy, the cctbx is organized as a set of smaller
modules. This is very much like a village (the cctbx project) with individual
houses (modules) for each family (groups of developers, of any size including
one).

The cctbx code base is available without restrictions and free of charge to all
interested developers, both academic and commercial. The entire community is
invited to actively participate in the development of the code base. A
sophisticated technical infrastructure that enables community based software
development is provided by SourceForge_. This service is also free of charge and
open to the entire world.

The cctbx is designed with an open and flexible architecture to promote
extendability and easy incorporation into other software environments. The
package is organized as a set of ISO C++ classes with Python_ bindings. This
organization combines the computational efficiency of a strongly typed compiled
language with the convenience and flexibility of a dynamically typed scripting
language in a strikingly uniform and very maintainable way.

Use of the Python interfaces is highly recommended, but optional. The cctbx can
also be used purely as a C++ class library.

.. contents:: Table of Contents

High level organization
=======================

The SourceForge cctbx project currently contains these modules:

libtbx
------

The build system common to all other modules. This includes a very thin wrapper
around the SCons_ software construction tool.

:ref:`API Documentation for libtbx <libtbx>`

boost_adaptbx
-------------

A very small adaptor toolbox with platform-independent instructions for building
the Boost.Python_ library.

:ref:`API Documentation for boost_adaptbx <boost_adaptbx>`

scitbx
------

Libraries for general scientific computing (i.e. libraries that are not specific
to crystallographic applications): a family of high-level C++ array types, a
fast Fourier transform library, and a C++ port of the popular L-BFGS
quasi-Newton minimizer, all including Python_ bindings. These libraries are
separated from the crystallographic code base to make them easily accessible for
non-crystallographic application developers.

:ref:`API Documentation for scitbx <scitbx>`

cctbx
-----

Libraries for general crystallographic applications, useful for both
small-molecule and macro-molecular crystallography. The libraries in the cctbx
module cover everything from algorithms for the handling of unit cells to
high-level building blocks for refinement algorithms. Note the distinction
between the cctbx *project* and the cctbx *module*. In retrospect we should have
chosen a different name for the project, but the current naming reflects how the
modules have evolved and it would be too disruptive to start a grand renaming.

:doc:`API Documentation for cctbx <cctbx/cctbx>`

iotbx
-----

Libraries for reading and writing established file formats.


:ref:`API Documentation for iotbx <iotbx>`

mmtbx
-----

Libraries for macromolecular crystallography.

:ref:`API Documentation for mmtbx <mmtbx>`

Tour
====

Tour of the :ref:`cctbx <tour>`.

Installation
============

:ref:`Installation instructions <installation>` for both binary installation and
installation from sources.

The cctbx :ref:`build system <build_system>` is based on SCons_.

Reference Documentation
=======================

`cctbx C++ interfaces <http://cctbx.sourceforge.net/current/c_plus_plus/namespaces.html>`_

Most documented C++ interfaces are also available at the Python
layer. Unfortunately the documentation tools available are not capable of
merging the documentations. Therefore Python users need to also consult the C++
documention.

Links
=====

* cctbx - Automatic multi-platform builds
* cctbx - Public SVN repository
* cctbx - Project page at SourceForge
* cctbx - Development history

Acknowledgments
===============

We would like to thank David Abrahams for creating the amazing Boost.Python_
library and for patiently supporting the entire open source community. We would
like to thank Airlie McCoy for allowing us to adapt some parts of the Phaser
package (FFT structure factor calculation). Kevin Cowtan has contributed
algorithms for the handling of reciprocal space asymmetric units. We are also
grateful for his development of the `Clipper <http://www.ysbl.york.ac.uk/~cowtan/clipper/clipper.html>`_ library from which we have adapted
some source code fragments. Our work was funded in part by the US Department of
Energy under Contract No. DE-AC03-76SF00098. We gratefully acknowledge the
financial support of NIH/NIGMS.

The cctbx SVN development tree is hosted by SourceForge

Contact
=======

`cctbx@cci.lbl.gov <mailto:cctbx@cci.lbl.gov>`_


.. _Boost: http://www.boost.org/

.. _Boost.Python: http://www.boost.org/libs/python/doc/

.. _Clipper: http://www.ysbl.york.ac.uk/~cowtan/clipper/clipper.html

.. _PHENIX: http://phenix-online.org/

.. _Python: http://www.python.org/

.. _SCons: http://www.scons.org/

.. _SourceForge: http://sourceforge.net/

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
