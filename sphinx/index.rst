
+++++++++++++++++++++++++++++++++++++
Computational Crystallography Toolbox
+++++++++++++++++++++++++++++++++++++

.. _introduction:
.. toctree::
   :maxdepth: 1

   installation
   tour
   history
   libtbx/index
   boost_adaptbx/index
   iotbx/index
   scitbx/index
   cctbx/cctbx
   mmtbx/mmtbx


Welcome to CCTBX's documentation!
=================================

.. contents:: Table of Contents

High level organization
=======================

The SourceForge cctbx project currently contains these modules.  The core
libraries required for most other applications are libtbx, boost_adaptbx,
scitbx, cctbx, and usually iotbx.  Functionality specific to macromolecules
and small molecules lives in mmtbx and smtbx, respectively.

libtbx
------

The build system common to all other modules. This includes a very thin wrapper
around the SCons_ software construction tool.  It also contains many useful
frameworks and utilities to simplify application development, including tools
for regression testing, parallelization across multiprocessor systems and
managed clusters, and a flexible, modular configuration syntax called PHIL
(Python Hierarchial Interface Language) used throughout the CCTBX.

:ref:`API Documentation for libtbx <libtbx>`

boost_adaptbx
-------------

A very small adaptor toolbox with platform-independent instructions for building
the Boost.Python_ library.

:ref:`API Documentation for boost_adaptbx <boost_adaptbx>`

scitbx
------

Libraries for general scientific computing (i.e. libraries that are not specific
to crystallographic applications).  This includes a family of high-level C++
array types, a fast Fourier transform library, and a C++ port of the popular
L-BFGS quasi-Newton minimizer, and many mathematical utilities, all including
Python_ bindings. These libraries are
separated from the crystallographic code base to make them easily accessible for
non-crystallographic application developers.

:ref:`API Documentation for scitbx <scitbx>`

cctbx
-----

Libraries for general crystallographic applications, useful for both
small-molecule and macro-molecular crystallography. The libraries in the cctbx
module include algorithms and data structures for the handling
of crystal symmetry, basic geometry restraints, reflection data, atomic
displacement parameters, X-ray scattering, and high-level building blocks for
refinement algorithms.
Note the distinction between the CCTBX *project* and the cctbx *module*.

:doc:`API Documentation for cctbx <cctbx/cctbx>`

iotbx
-----

Libraries for reading and writing common file formats, including PDB, CIF,
many reflection formats, electron density maps, and sequences.

:ref:`API Documentation for iotbx <iotbx>`

mmtbx
-----

Functionality specific to macromolecular crystallography.  This includes
all of the machinery required for setup of geometry restraints, bulk solvent
correction and scaling, analysis of macromolecular diffraction data,
calculation of weighted map coefficients, and most of the methods implemented
in phenix.refine.  The majority of infrastructure for the MolProbity
validation server (and Phenix equivalent) is also located here.

:doc:`API Documentation for mmtbx <mmtbx/mmtbx>`

xfel
-----

Software for processing serial data collected using an X-ray free electron laser.
Includes spotfinding, integration, data clustering/filering and merging tools.

:doc:`API Documentation for xfel <xfel/xfel>`

smtbx
-----

Functionality specific to small-molecule crystallography, including a complete
refinement program (smtbx.refine).

:doc:`API Documentation for smtbx <smtbx>`

dxtbx
-----

The Diffraction Image Toolbox, a library for handling X-ray detector data
of arbitrary complexity from a variety of standard formats.  (Also used by
routines in iotbx.)

:doc:`API documentation for dxtbx <dxtbx>`

Tour
====

* Tour of the :ref:`cctbx <tour>`.
* Tour of the `scitbx (by Michael Hohn) <http://cci.lbl.gov/~hohn/scitbx-tour.html>`_
* Tour of the `array_family (by Michael Hohn) <http://cci.lbl.gov/~hohn/array-family-tour.html>`_

Tutorials
=========

* `Free-electron laser data processing (cctbx.xfel) <http://cci.lbl.gov/xfel>`_
* `IUCr 2008 (Software Fayre) (rigid_body_refinement_core.py)
  <http://cctbx.sourceforge.net/iucr2008/rigid_body_refinement_core.html>`_
* `SBGrid 2008 (Quo Vadis) (iotbx.pdb)
  <http://cctbx.sourceforge.net/sbgrid2008/tutorial.html>`_
* `Siena 2005 IUCr Crystallographic Computing School
  <http://cctbx.sourceforge.net/siena2005/>`_
* `scitbx/rigid_body/essence subset
  <http://cctbx.sourceforge.net/scitbx_rigid_body_essence/>`_

Installation
============

:ref:`Installation instructions <installation>` for both binary installation and
installation from sources.

The cctbx build system is based on SCons_.

Reference Documentation
=======================

`cctbx C++ interfaces <http://cctbx.sourceforge.net/current/c_plus_plus/namespaces.html>`_

Most documented C++ interfaces are also available at the Python
layer. Unfortunately the documentation tools available are not capable of
merging the documentations. Therefore Python users need to also consult the C++
documention.

Links
=====

* `cctbx - Automatic multi-platform builds <http://cci.lbl.gov/cctbx_build/>`_
* `cctbx - Public git repository <https://github.com/cctbx/cctbx_project>`_

Contact
=======

`cctbx@cci.lbl.gov <mailto:cctbx@cci.lbl.gov>`_

There is also a mailing list (cctbxbb) where developers can answer questions from
the community. To join/manage your subscription, please visit http://www.phenix-online.org/mailman/listinfo/cctbxbb.


.. _Boost: http://www.boost.org/

.. _Boost.Python: http://www.boost.org/libs/python/doc/

.. _Clipper: http://www.ysbl.york.ac.uk/~cowtan/clipper/clipper.html

.. _PHENIX: http://phenix-online.org/

.. _Python: http://www.python.org/

.. _SCons: http://www.scons.org/

.. _GitHub: https://github.com

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
