==============================================================
scitbx - general-purpose scientific programming infrastructure
==============================================================

.. _scitbx:

The :py:mod:`scitbx` module contains many of the core library routines
required for any computational science project
(i.e. libraries that are not specific to crystallographic applications):
a family of high-level C++ array types, a fast Fourier transform library,
gradient-driven optimization algorithms (:py:mod:`scitbx.lbfgs` and
:py:mod:`scitbx.lstbx`), matrix manipulation (:py:mod:`scitbx.matrix`,
among others), and a variety of general-purpose
mathematical functions (primarily in :py:mod:`scitbx.math`).
These libraries are separated from the
crystallographic code base to make them easily accessible for
non-crystallographic application developers.

Submodules
==========

.. toctree::
    :maxdepth: 1

    scitbx.array_family
    scitbx.matrix
