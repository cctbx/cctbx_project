
.. _iotbx:


iotbx - file readers and writers
================================

The iotbx module contains most tools for reading and writing the standard
formats used by both macromolecular and small-molecule crystallographers,
including PDB, CIF, MTZ, and various other file types.  For some formats the
resulting data will be encapsulated in objects defined in :py:mod:`cctbx`
and/or :py:mod:`scitbx`; others have custom classes, in particular PDB/mmCIF
files which have their own complex internal representation independent of the
X-ray scattering properties.

Subpackages
-----------

.. toctree::
    iotbx.managers
    iotbx.file_reader
    iotbx.cif
    iotbx.pdb
    iotbx.shelx
    iotbx.xds
    iotbx.reflection_file_reader
