
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
   mmtbx/mmtbx


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
development is provided by GitHub_. This service is also free of charge and
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

Other
-----

Many additional libraries have more specialized functionality, including:

* spotfinder - fast detection of Bragg peaks in diffraction images
* ucif - the core CIF I/O library (used by iotbx)
* rstbx - Reciprocal Space Toolbox, used for data processing
* gltbx - OpenGL bindings, including a wxPython-based viewer framework
* crys3d - Modules for the display of molecules, electron density, and reciprocal space data
* fable - a program (and compatibility library) for porting Fortran77 to C++
* wxtbx - wxPython controls used in the Phenix GUI and various utilities
* cma_es - a library of derivative-free optimization methods

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


Newsletter articles and examples
================================

* `IUCr Computing Commission No. 1, 2003/01 <http://cci.lbl.gov/publications/download/iucrcompcomm_jan2003.pdf>`_
   - State of the Toolbox: an overview of the Computational Crystallography Toolbox (CCTBX)
* `IUCr Computing Commission No. 2, 2003/07 <http://cci.lbl.gov/publications/download/iucrcompcomm_jul2003.pdf>`_
   - Fast triplet generator for direct methods
   - Gallery of direct-space asymmetric units
* `IUCr Computing Commission No. 3, 2004/01 <http://cci.lbl.gov/publications/download/iucrcompcomm_jan2004.pdf>`_
   - Reduced cell computations
   - Determination of lattice symmetry
   - N-Gaussian approximations to scattering factors
   - Fast structure-factor gradients
   - Universal reflection file reader
* `IUCr Computing Commission No. 4, 2004/08 <http://cci.lbl.gov/publications/download/iucrcompcomm_aug2004.pdf>`_
   - Geometry restraints
   - Bulk solvent correction and scaling
* `IUCr Computing Commission No. 5, 2005/01 <http://cci.lbl.gov/publications/download/iucrcompcomm_jan2005.pdf>`_
   - Phil and friends (:doc:`latest Phil documentation <libtbx/libtbx.phil>`)
   - Refinement tools
   - Reflection statistics
   - Double coset decomposition
   - iotbx.mtz
* `IUCr Computing Commission No. 6, 2005/09 <http://cci.lbl.gov/publications/download/iucrcompcomm_sep2005.pdf>`_
   - See also `Sienna 2005 tutorials <http://cctbx.sourceforge.net/siena2005/>`_
* `IUCr Computing Commission No. 7, 2006/11 <http://cci.lbl.gov/publications/download/iucrcompcomm_nov2006.pdf>`_
   - iotbx.pdb
   - mmtbx.alignment
   - scitbx.math.superpose
   - mmtbx.super
* `IUCr Computing Commission No. 8, 2007/11 <http://cci.lbl.gov/publications/download/iucrcompcomm_nov2007.pdf>`_
   - Refinement tools for small-molecule crystallographers
   - Phil developments (:doc:`latest Phil documentation <libtbx/libtbx.phil>`)
* `IUCr Computing Commission No. 8, 2008/10 <http://cci.lbl.gov/publications/download/iucrcompcomm_oct2008.pdf>`_
   - Symmetry in crystallographic applications
* `IUCr Computing Commission No. 8, 2009/11 <http://cci.lbl.gov/publications/download/iucrcompcomm_nov2009.pdf>`_
   - Experience converting a large Fortran-77 program to C++
* `CCP4 newsletter No. 42, Summer 2005: The Phenix refinement framework <http://www.ccp4.ac.uk/newsletters/newsletter42/articles/Afonine_GrosseKunstleve_Adams_18JUL2005.doc>`_
* `CCP4 newsletter No. 42, Summer 2005: Characterization of X-ray data sets <http://www.ccp4.ac.uk/newsletters/newsletter42/articles/CCP4_2005_PHZ_RWGK_PDA.doc>`_
* `CCP4 newsletter No. 43, Winter 2005: Xtriage and Fest: automatic assessment of X-ray data and substructure structure factor estimation <http://www.ccp4.ac.uk/newsletters/newsletter43/articles/PHZ_RWGK_PDA.pdf>`_
* `CCP4 newsletter No. 44, Summer 2006: Exploring Metric Symmetry <http://www.ccp4.ac.uk/newsletters/newsletter44/articles/explore_metric_symmetry.html>`_

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

* `cctbx - Automatic multi-platform builds <http://cci.lbl.gov/cctbx_build/>`_
* `cctbx - Public git repository <https://github.com/cctbx/cctbx_project>`_

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

The cctbx git development tree is hosted by GitHub.

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
