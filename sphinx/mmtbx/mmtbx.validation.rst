========================
mmtbx.validation package
========================

The ``validation`` module combines all aspects of model validation, both with
respect to geometry and against experimental data.  Many of these were adapted
from the `MolProbity web server <http://molprobity.biochem.duke.edu>`_, and
continue to be used for that purpose.  They are also available in the Phenix
GUI as a standalone program and as an accessory to phenix.refine.  However, the
individual analyses may also be run separately and used to guide various
decision-making during model-building and refinement.

Base classes
------------

With few exceptions, all analyses in the validation framework inherit from a
common set of base classes (or use them internally).  This provides a unified
API for accessing similar information.

.. automodule:: mmtbx.validation
    :members:
    :undoc-members:
    :show-inheritance:

-----------
Subpackages
-----------

.. toctree::

    mmtbx.validation.molprobity

----------
Submodules
----------

Nearly a dozen different analyses may be performed.  Note that many of these
require additional programs and/or data not distributed with CCTBX except in
the context of `Phenix <http://www.phenix-online.org>`_:

- **Reduce** and **Probe**: standalone C++ programs for adding hydrogens and
  analyzing atomic contacts, respectively.  These are available from the
  `Richardson lab <http://kinemage.biochem.duke.edu`>_.
- **suitename**: standalone C program for analyzing RNA geometry, also from the
  Richardson lab.
- **Geometry restraints**: this can be substituted by the standard CCP4 monomer
  library, but Phenix includes its own set of restraints (partially derived
  from CCP4's).
- **Ramachandran and rotamer distributions**: these contain frequencies for
  each conformation based on a library of high-quality X-ray structures.

These are all essentially freely available; contact the developers if you
require specific files.  If you are using a Phenix installation to perform
CCTBX development, you will already have access to everything necessary.

Rotalyze - protein sidechain rotamer analysis
---------------------------------------------

Also available as a standalone program, ``phenix.rotalyze``.

.. automodule:: mmtbx.validation.rotalyze
    :members:
    :undoc-members:
    :show-inheritance:

Ramalyze - Ramachandran plot analysis
-------------------------------------

Also available as a standalone program, ``phenix.ramalyze``.

.. automodule:: mmtbx.validation.ramalyze
    :members:
    :undoc-members:
    :show-inheritance:

Clashscore - all-atom contacts using Reduce and Probe
-----------------------------------------------------

Also available as a standalone program, ``phenix.clashscore``.

.. automodule:: mmtbx.validation.clashscore
    :members:
    :undoc-members:
    :show-inheritance:

C-beta deviations
-----------------

Also available as a standalone program, ``phenix.cbetadev``.

.. automodule:: mmtbx.validation.cbetadev
    :members:
    :undoc-members:
    :show-inheritance:

mmtbx.validation.restraints module
----------------------------------

.. automodule:: mmtbx.validation.restraints
    :members:
    :undoc-members:
    :show-inheritance:

mmtbx.validation.rna_validate module
------------------------------------

.. automodule:: mmtbx.validation.rna_validate
    :members:
    :undoc-members:
    :show-inheritance:

mmtbx.validation.analyze_peptides module
----------------------------------------

.. automodule:: mmtbx.validation.analyze_peptides
    :members:
    :undoc-members:
    :show-inheritance:

mmtbx.validation.experimental module
------------------------------------

.. automodule:: mmtbx.validation.experimental
    :members:
    :undoc-members:
    :show-inheritance:

mmtbx.validation.graphics module
--------------------------------

.. automodule:: mmtbx.validation.graphics
    :members:
    :undoc-members:
    :show-inheritance:

mmtbx.validation.ligands module
-------------------------------

.. automodule:: mmtbx.validation.ligands
    :members:
    :undoc-members:
    :show-inheritance:

mmtbx.validation.model_properties module
----------------------------------------

.. automodule:: mmtbx.validation.model_properties
    :members:
    :undoc-members:
    :show-inheritance:

mmtbx.validation.sequence module
--------------------------------

.. automodule:: mmtbx.validation.sequence
    :members:
    :undoc-members:
    :show-inheritance:

mmtbx.validation.utils module
-----------------------------

.. automodule:: mmtbx.validation.utils
    :members:
    :undoc-members:
    :show-inheritance:

mmtbx.validation.waters module
------------------------------

.. automodule:: mmtbx.validation.waters
    :members:
    :undoc-members:
    :show-inheritance:
