"""endo_exo: QM-region builder engine.

Grows a QM region around seed sites by breadth-first traversal of the
covalent graph and caps dangling bonds with hydrogens.  The
:class:`~mmtbx.geometry_restraints.endo_exo.builder.QMRegionBuilder`
orchestrates the per-seed pipeline; the individual engine classes live in
the sibling modules (seeds, graph, cutting, grow, capping).

This package re-exports nothing; submodules are imported directly, as in
``from mmtbx.geometry_restraints.endo_exo.builder import QMRegionBuilder``,
so that a single engine piece can be used without importing the whole
pipeline.
"""

from __future__ import absolute_import, division, print_function
