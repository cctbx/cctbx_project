# Probe2

This is a CCTBX re-implementation of the Richardson Labs' Probe program from https://github.com/rlabduke/probe
It is named probe2 so that it can be built and used alongside the original probe (which can still be built as
a module) during regression testing.  It is intended that future development and maintenance will happen
on Probe2 and not on the original code base.

## C++ Classes

The C++ classes, wrapped for use in Python, make use of CCTBX and Boost structures and define all of
their classes and functions in the molprobity::probe C++ namespace.  These are wrapped in Python and
can be imported as **mmtbx_probe_ext**.

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
and is a combination of the codes from the original Probe and original Reduce codes.  The **DotScorer** class is the
workhorse, with the other classes providing structured input and output for this class.

* **boost_python/probe_pbl.cpp:** As is usual for CCTBX projects, this file contains the wrapper code for the
above functions that uses Boost Python to wrap them to be called from Python.  It places 

* **tst_probe.cpp** Compiles into a program that tests all of the C++ methods.  It is compiled by CMakeLists.txt
but not by SConscript, so it is not normally used and is not deployed within CCTBX.

## Python Modules

As mentioned above, the C++ classes can be imported using the **mmtbx_probe_ext** module.

There are also a set of Python classes that provide support for applications that want to use the C++ classes.

* **AtomTypes.py:** This contains classes and functions that help determine the characteristics of atoms, which is needed
to know how to handle them during probe calculations.  Some of these tables contain information from the
original Probe and Reduce code that has been superceded by the same information from various CCTBX tables,
but one of them contains information not found elsewhere in CCTBX: the **IsAromatic()** function tells whether
a specified atom from a specified residue is part of an aromatic ring in a standard residue.

* **Helpers.py.** This contains helper functions needed by both Probe2 and Reduce2.  See the file itself for a complete
list of functions and parameters.  **getBondedNeighborLists()** converts bond proxy information into a dictionary of
lists, providing a list of bonded atoms looked up by atom.  **compatibleConformations()** tells whether two atoms
are in compatible conformations.  **getAtomsWithinNBonds()** returns the list of atoms from compatible conformations
that are bonded within N hops from a specified atom.  **getExtraAtomInfo()** looks up extra information needed to
determine interactions and returns it encapsulated in a C++ structure reference.
**getPhantomHydrogensFor()** produces a list of potential hydrogens for a water oxygen that point
towards nearby acceptors (or all nearby atoms).  **fixupExplicitDonors()** adjusts the extra atom information for
a set of atoms once hydrogens have been explicitly added to the model, adjusting the donor status.  **rvec3()**
and **lvec3()** return **scitbx.matrix.rec** elements for left-side and right-side multiplication, enabling the use
of these built-in C++ types' methods from within Python.

## Testing

The file **ener_lib_molprobity.cif** is a drop-in replacement for the ener_lib.cif file that is found in
the chem_data/mon_lib directory under modules in a CCTBX build.  It has atomic radii consistent with
the original Probe code to enable regression tests between Probe2 and Probe.  It should not be used in
production once the CCTBX radii are adjusted (this is ongoing as of 10/5/2021).  Note that replacing this
file will impact other CCTBX programs that use atom radii.

Regression tests were performed between Probe2 and the original Probe code in October 2021 before the release
of Probe2:

