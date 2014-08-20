mmtbx - macromolecular crystallography
======================================

The :py:mod:`mmtbx` module has two purposes: firstly, it introduces additional
crystallographic algorithms that are traditionally specific to macromolecular
structures, such as bulk-solvent correction, absolute scaling, maximum
likelihood targets, non-crystallographic symmetry, weighted maps, and
real-space refinement; and secondly, it provides many
higher-level classes for querying and manipulating molecular structures and
writing full-featured macromolecular crystallography applications.  Most of the
core functions of the program `phenix.refine
<http://www.phenix-online.org/documentation/reference/refinement.html>`_ are
part of mmtbx, and the module also includes the complete code for the
programs ``phenix.xtriage``, ``phenix.maps``, ``phenix.ensemble_refinement``,
``phenix.model_vs_data``, and ``phenix.pdbtools``, plus the backend logic for
the `MolProbity server <http://molprobity.biochem.duke.edu>`_.

For new developers, the module :py:mod:`mmtbx.command_line` provides a
simplified frontend for running routine setup tasks (such as bulk solvent
correction or geometry restraints interpretation) and populating many of the
standard data structures.  (The source code also provides a suitable example
should you wish to instantiate these objects yourself.)  However, since the
resulting objects tend to be especially complex, some knowledge of their
internals and purpose is helpful.

The most important modules in :py:mod:`mmtbx` can be summarized as follows:

- :py:mod:`mmtbx.f_model`: encapsulates most functionality for calculating
  F(model) (i.e. F(calc)) and comparing it against experimental data.  The
  calculation of weighted electron density maps usually starts from this
  point.
- :py:mod:`mmtbx.monomer_library.pdb_interpretation`: processing of models
  (in PDB or mmCIF format) and interpretation with respect to a library of
  geometry restraints (provided as part of Phenix, but the CCP4 Monomer
  Library may also be used).
- :py:mod:`mmtbx.scaling`: absolute scaling routines plus analyses of
  experimental data to assess quality and detect possible pathologies such as
  twinning.  This module is slightly misnamed, since it is primarily used as
  the backend for ``Xtriage``.
- :py:mod:`mmtbx.validation`: tools for assessing model quality, based either
  on geometric criteria or fit to data.  Also see :py:mod:`mmtbx.rotamer`, the
  lower-level code for analyzing conformation.
- :py:mod:`mmtbx.refinement`: routines for model optimization against either
  reciprocal space or real-space (i.e. electron density map) data.
- :py:mod:`mmtbx.pdbtools`: a small library (as well as a standalone program)
  for various common model manipulations.

Subpackages
-----------

.. toctree::
    mmtbx.alignment
    mmtbx.building
    mmtbx.bulk_solvent
    mmtbx.cablam
    mmtbx.chemical_components
    mmtbx.command_line
    mmtbx.conformation_dependent_library
    mmtbx.den
    mmtbx.density_modification
    mmtbx.dynamics
    mmtbx.f_model
    mmtbx.f_model_info
    mmtbx.find_peaks
    mmtbx.geometry
    mmtbx.geometry_restraints
    mmtbx.grow_density
    mmtbx.hydrogens
    mmtbx.ias
    mmtbx.invariant_domain
    mmtbx.ions
    mmtbx.kinemage
    mmtbx.lattice
    mmtbx.map_tools
    mmtbx.maps
    mmtbx.masks
    mmtbx.max_lik
    mmtbx.model
    mmtbx.model_statistics
    mmtbx.model_vs_data
    mmtbx.monomer_library
    mmtbx.msa
    mmtbx.ncs
    mmtbx.pdbtools
    mmtbx.polygon
    mmtbx.real_space
    mmtbx.real_space_correlation
    mmtbx.refinement
    mmtbx.restraints
    mmtbx.rotamer
    mmtbx.rsr
    mmtbx.scaling
    mmtbx.secondary_structure
    mmtbx.solvent
    mmtbx.tls
    mmtbx.ligands
    mmtbx.torsion_restraints
    mmtbx.twinning
    mmtbx.utils
    mmtbx.validation
    mmtbx.wwpdb
    mmtbx.disorder

Submodules
----------

.. toctree::

    mmtbx.resolve_resources
    mmtbx.xmanip_tasks
    mmtbx.arrays
    mmtbx.pdb_distances
    mmtbx.pdb_symmetry
    mmtbx.superpose
    mmtbx.xmanip

Module contents
---------------

.. automodule:: mmtbx
    :members:
    :undoc-members:
    :show-inheritance:
