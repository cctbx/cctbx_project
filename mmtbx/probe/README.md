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
to know how to handle them during Probe2 calculations.  The **IsAromatic()** function tells whether
a specified atom from a specified residue is part of an aromatic ring in a standard residue.

* **Helpers.py.** This contains helper functions needed by both Probe2 and Reduce2.  See the file itself for a complete
list of functions and parameters.  **getBondedNeighborLists()** converts bond proxy information into a dictionary of
lists, providing a list of bonded atoms looked up by atom.  **compatibleConformations()** tells whether two atoms
are in compatible conformations.  **getAtomsWithinNBonds()** returns the list of atoms from compatible conformations
that are bonded within N hops from a specified atom.  **getExtraAtomInfo()** looks up extra information needed to
determine interactions and returns it encapsulated in a C++ structure reference.  **dihedralChoicesForRotatableHydrogens()**
determines which of the Hydrogens passed in should be paired with which of the potential dihedral-determining
atoms passed in to compute the dihedral angle for the Hydrogen; this is used for placement and angle description.
**getPhantomHydrogensFor()** produces a list of potential hydrogens for a water oxygen that point
towards nearby acceptors (or all nearby atoms).  **fixupExplicitDonors()** adjusts the extra atom information for
a set of atoms once hydrogens have been explicitly added to the model, adjusting the donor status.  **rvec3()**
and **lvec3()** return **scitbx.matrix.rec** elements for left-side and right-side multiplication, enabling the use
of these built-in C++ types' methods from within Python.

## Testing

Unit tests are implemented in each C++ class and Python module to verify that they perform as intended.
The tests performed include the following:
* **AtomTypes.py:** When it is run as a script, this file will raise an assertion on failure or exit
normally on success.  The **IsAromatic** function is tested by calling it with a small number of atoms
within residues to ensure that it markes aromatic ones and does not mark others.
* **Helpers.py:** When it is run as a script, this file will raise an assertion on failure or exit
normally on success.  It normally operates on a generated model, but can be run with the name of
a PDB or CIF file on the command line to specify a different model to use.  The script tests the
following internal functions:
    * **getExtraAtomInfo()** is called with useNeutronDistances and useImplicitHydrogenDistances
    both set to False, and with each set to True.  The radius, acceptor and donor status are compared
    against expected results for several atom types.  The abiliy to set and test the dummy-hydrogen
    status is tested for all atoms.  It is later tested against either a generated snippet or a PDB
    or CIF file specified on the command line to ensure it can work on a loaded model.
    * **getPhantomHydrogensFor()** is tested against a hand-constructed geometry that does not match
    any standard residue but has an acceptor with three nearby atoms, two of which are acceptors.
    Placement is run with different occupancy thresholds and acceptorOnly (placedHydrogenRadius is
    not varied from the default).  It is verified that we get as many phantom Hydrogens as expected
    for each case and that all of them point towards non-Oxygen atoms.
    * **getBondedNeighborLists()** uses a specific snippet from 1xso for which we know the correct
    bonding relationship.  We verify that the number of neighbors of each atom matches what is expected.
    * **compatibleConformations()** is tested by making atoms in alternates "", "A", and "B"
    and verifying that the compatible ones are marked as such and the incompatible ones are not.
    * **isPolarHydrogen**() is tested using a PDB snipped from hydrogenated 4fen that we know the
    answer for and verifying that all hydrogens are correctly marked as polar or non-polar.
    * **getAtomsWithinNBonds()** is tested by running it on a PDB snippet, both with a limited
    number of steps for non-Hydrogen cases and without, both starting on a Hydrogen and not.  The
    counts are verified against hand-counted expected values.
    * **dihedralChoicesForRotatableHydrogens()** is tested to make sure it selects the correct atoms
    from a PDB snippet, and to verify that it raises an exception when it does not have the information
    needed to determine this.
    * **rvec3()** and **lvec3()** are tested to ensure that they produce vectors that can be subtracted
    to form proper-length results.  These verify that the underlying scitbx.matrix classes are working
    as expected.
* **SpatialQuery.cpp:** exposes a **SpatialQuery_test()** method that returns an empty string on success
and an error message on failure.  It runs the following tests on the SpatialQuery structure that Probe2
uses to determine which atoms are near a point in space:
    * Constructed number of bins in X, Y, and Z when creating a particular region in space with a
    specified bin size matches expectations.
    * Inverted upper and lower grid points swaps the order of the edges; the region is smaller than
    the bin size and it is tested that a single grid point is constructed.
    * Atom addition and removal is tested on an empty grid, ensuring that an atom can be inserted
    exactly once and can be removed exactly once.
    * Atoms are inserted in a multi-element grid and it is verified that only the atom near a specified
    grid point is returned.  Also, the use of a minimum distance it tested to ensure that no atoms
    are returned when the atom is too close to the query location.
    * A synthetic model consisting of many hydrogens on a grid is used to construct a query structure
    from a model.   It is verified that the expected number of grid points are created based on the
    size of the structure and the default bin size.
    * In a query that covers much more than the size of the entire grid, it is verified that all
    atoms are returned for queries either inside the grid or outside the grid on all sides.
    * A test is run for the return of the correct number of neighbors when a query is made that includes
    multiple bins, each of which includes an atom.


A **tst_probe.py** script is located in the mmtbx/regression folder within this project.  this
script runs all of the unit tests described above and then computes a probe score on a generated
small molecule (or on a PDB or CIF file is one is specified on the command line).  This program raises
an assertion if it finds a problem and exits normally when everything is okay.  This verifies that
the Python/C++ linkage is working correctly and that the code works on a standard model.

The file **ener_lib_molprobity.cif** is a drop-in replacement for the ener_lib.cif file that is found in
the chem_data/mon_lib directory under modules in a CCTBX build.  It has atomic radii consistent with
the original Probe code to enable regression tests between Probe2 and Probe.  It should not be used in
production once the CCTBX radii are adjusted (this is ongoing as of 10/5/2021).  Note that replacing this
file will impact other CCTBX programs that use atom radii.

Regression tests were performed between Probe2 and the original Probe code in October 2021 before the release
of Probe2:

**@todo**