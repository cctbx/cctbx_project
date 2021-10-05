# Probe2

This is a CCTBX re-implementation of the Richardson Labs' Probe program from https://github.com/rlabduke/probe
It is named probe2 so that it can be built and used alongside the original probe (which can still be built as
a module) during regression testing.  It is intended that future development and maintenance will happen
on Probe2 and not on the original code base.

## C++ classes

The C++ classes, wrapped for use in Python, make use of CCTBX and Boost structures and define all of
their classes and functions in the molprobity::probe namespace.

* **Common.h:** Definition of the type of floating-point number to be used as a coordinate (Coord) and of
a three-dimensional **Point** representation to be used internally.

* **SpatialQuery.*:** Defines a **SpatialQuery** class that can be used to quickly determine how many atoms
are within a specified minimum and maximum distance from a point in space.  Atoms can be dynamically added
and removed from the structure to support rapid queries in a system where some atoms are moving (this is
needed by Reduce2).

* **DotSpheres.*:** Defines a **DotSphere** class that will evenly place dots on the surface of a sphere with
specified radius and a **DotSphereCache** class that constructs DotSpheres with specified radii, only making
a single sphere for any radius.

* **Scoring.*:** Defines a set of classes and helper functions that can be used to compute interactions between
single dots or sets of dots and nearby atoms.  This is used by both Probe2 and Reduce2 to perform dot scoring
and is a combination of the codes from the original Probe and original Reduce codes.  The **DotScorer** is the
workhorse, with the other classes providing structured input and output for this class.
