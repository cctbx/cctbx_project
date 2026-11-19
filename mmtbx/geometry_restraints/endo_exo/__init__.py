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

Used as a library, it takes a model and returns the regions, without files
or file names::

    from mmtbx.geometry_restraints.endo_exo.builder import (
      QMRegionBuilder, default_params)

    params = default_params()
    params.buffer.radius = 4.0
    for region in QMRegionBuilder(params, logger=log).run(model):
      region.label, region.model, region.cap_iseqs  # see QMRegionBuilder.run

``mmtbx.programs.endo_exo`` is the command-line program on top: it names
the regions and writes them.
"""
