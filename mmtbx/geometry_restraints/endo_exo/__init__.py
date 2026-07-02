"""endo_exo -- QM-region builder engine.

Grows a QM region around seed sites by breadth-first traversal of the
covalent graph, caps dangling bonds with hydrogens, and estimates the net
charge of the region.  The :class:`~mmtbx.geometry_restraints.endo_exo.builder.QMRegionBuilder`
orchestrates the per-seed pipeline; the individual engine classes live in
the sibling modules (seeds, graph, cutting, grow, capping, charge).

Submodules are imported directly (e.g. ``from
mmtbx.geometry_restraints.endo_exo.builder import QMRegionBuilder``) so that
the lightweight engine pieces can be used without pulling in the full
pipeline.
"""

from __future__ import absolute_import, division, print_function
