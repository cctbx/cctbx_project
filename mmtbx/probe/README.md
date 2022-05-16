# Probe2

This is a CCTBX re-implementation of the Richardson Labs' Probe program from https://github.com/rlabduke/probe
It is named probe2 so that it can be built and used alongside the original probe (which can still be built as
a module) during regression testing.  It is intended that future development and maintenance will happen
on Probe2 and not on the original code base.

**Notes:**
* Probe2 requires the chem_data modules from Phenix.

## C++ Classes

The C++ classes, wrapped for use in Python, make use of CCTBX and Boost structures and define all of
their classes and functions in the molprobity::probe C++ namespace.  These are wrapped in Python and
can be imported as **mmtbx_probe_ext**.  As described below, the use of the Helpers.create*() functions
is preferred to constructing these objects directly because it handles all Probe Phil parameters, including
any that are added in the future.

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
above functions that uses Boost Python to wrap them to be called from Python.

* **tst_probe.cpp** Compiles into a program that tests all of the C++ methods.  It is compiled by CMakeLists.txt
but not by SConscript, so it is not normally used and is not deployed within CCTBX.

## Python Modules

As mentioned above, the C++ classes can be imported using the **mmtbx_probe_ext** module.

There are also a set of Python classes that provide support for applications that want to use the C++ classes.

* **AtomTypes.py:** This contains classes and functions that help determine the characteristics of atoms, which is needed
to know how to handle them during Probe2 calculations.  The **IsAromatic()** function tells whether
a specified atom from a specified residue is part of an aromatic ring in a standard residue.

* **Helpers.py.** This contains helper functions needed by both Probe2 and Reduce2.  See the file itself for a complete
list of functions and parameters.  Notable ones include:
    * **probe_phil_parameters**: These are a description of CCTBX Phil parameters that control the
behavior of Probe2 library functions.  They can be used by a client program to automatically update its
command-line arguments whenever the library is updated.  This can be done by adding
`from mmtbx.probe import Helpers; master_phil_str += Helpers.probe_phil_parameters` in the Program Template.
    * **create** functions for SpatialQuery, DotSphereCache, and DotScorer C++ objects take in the probe
Phil parameters as an argument and apply them to the object construction as needed.  These functions will be
updated as new Phil parameters (and corresponding object parameters) are added, so that a client program can
make use of new default and settable options without having to change their code.  For this reason, the
use of these functions is preferred to constructing the objects directly.
    * **getBondedNeighborLists()** converts bond proxy information into a dictionary of
lists, providing a list of bonded atoms looked up by atom.
    * **compatibleConformations()** tells whether two atoms
are in compatible conformations.
    * **getAtomsWithinNBonds()** returns the list of atoms from compatible conformations
that are bonded within N hops from a specified atom.
    * **getExtraAtomInfo()** looks up extra information needed to
determine interactions and returns it encapsulated in a C++ structure reference.
    * **dihedralChoicesForRotatableHydrogens()**
determines which of the Hydrogens passed in should be paired with which of the potential dihedral-determining
atoms passed in to compute the dihedral angle for the Hydrogen; this is used for placement and angle description.
    * **getPhantomHydrogensFor()** produces a list of potential hydrogens for a water oxygen that point
towards nearby acceptors (or all nearby atoms).
    * **fixupExplicitDonors()** adjusts the extra atom information for
a set of atoms once hydrogens have been explicitly added to the model, adjusting the donor status.
    * **rvec3()** and **lvec3()** return **scitbx.matrix.rec** elements for left-side and right-side multiplication, enabling the use
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
* **SpatialQuery.cpp:** exposes a **SpatialQuery_test()** function that returns an empty string on success
and an error message on failure.  It runs the following tests on the SpatialQuery class that Probe2
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
* **DotSpheres.cpp:** exposes a **DotSpheres_test()** function that returns an empty string on success
and an error message on failure.  It calls the test() methods on the DotSphere and DotSphereCache classes
that Probe2 uses to construct probe spheres around atoms, which test the following:
    * Sphere creation should construct zero radius and no dots with negative density and radius parameters.
    * Sphere with very small density should have a single dot.
    * The number of dots produced is close to the number requested based no density and radius.  It scales
    with the density.
    * Sphere with reasonable density produces dots within all 8 octants around the origin.
    * Cache creates a single sphere with the correct radius the first time it is called.
    * Cache called a second time with the same radius returns the same sphere.
    * Cache called a third time with a different radius gives a different sphere with correct radius.
* **Scoring.cpp:** exposes a **Scoring_test()** function that returns an empty string on success
and an error message on failure.  It tests the following:
    * The **atom_charge()** function returns the correct polarity and
    sign for the following values: "--", "-", "", "+", "++", "+2", "-1", "0".
    * The **closest_contact()** function produces a distance above surface that properly scales
    as the dot moves further from the atom center, negative inside and positive outside.
    * The **closest_contact()** function applied to a dot and atom at the same non-origin location
    with radius 1 produces a -1 distance and a projection that is 1 away from the atom center.
    * The **closest_contact()** function produces the same projected contact distance for all points in the
    cardinal and diagonal directions away an atom at the origin.
    * **DotScorer.check_dot()** can detect clashes, determine the expected cause, for synthetic probe
    locations finds only the expected type for various bumps: WorseOverlap, BadOverlap, Clash, SmallOverlap,
    NoOverlap, CloseContact, WideContact.
    * **DotScorer.check_dot()** reports Ignore overlap type and Invalid interaction type when the probe is
    not near any atom.
    * **DotScorer.score_dots()** Construct test cases with all combinations of charges and extra information,
    holding the radii of the neighbor atom and probe atom constant.  Do this in combination with adding or
    not adding an excluded atom that completely covers the neighbor atom.  Tests against all of these
    cases to ensure that the behavior is as expected in each case.
    * **DotScorer.score_dots()** Test behavior of weak hydrogen bonds and their interaction with dummy Hydrogens
    for all combination of weakHBonds, source and target being a dummy hydrogen.
    * **DotScorer.score_dots()** Sweep an atom from just touching to far away and make sure the attract
    curve is monotonically decreasing to 0.
    * **DotScorer.score_dots()** Test the setting of weights for the various subscores: the ratio of scores matches
    the ratio of weights for some entries and is linear for others.
    * **DotScorer.score_dots()** Test the setting of bond-gap distances.  Testing for bad bumps present
    exactly when expected.  Test with non-donor Hydrogens, uncharged donor Hydrogen, charged Hydrogen donor.
    * **DotScorer.score_dots()** Test the control of occupancy level.
    * **DotScorer.score_dots()** Test behavior when given invalid parameters.
    * **DotScorer.count_surface_dots()** behaves as expected for non-overlapping, fully-overlapping,
    and partially-overlapping atoms (partial overlap only tested to be between the others).

A **tst_probe.py** script is located in the mmtbx/regression folder within this project.  this
script runs all of the unit tests described above and then computes a probe score on a generated
small molecule (or on a PDB or CIF file is one is specified on the command line).  This program raises
an assertion if it finds a problem and exits normally when everything is okay.  This verifies that
the Python/C++ linkage is working correctly and that the code works on a standard model.

The  **mmtbx/programs/probe2.py** script uses the probe library routines above to construct kinemages
of interactions and surfaces and to generate summary scores for interactions within a model.  It is
a replacement for the **mmtbx.probe** program, producing similar outputs but taking different
command-line options.  The following tests are provided for this program:
* There is a Test() function defined within the module that will test all of its non-class functions.
To run it, the module can be imported and then probe2.Test() called.
It will fail with an assertion failure if there is a problem with the tests:
    * **_condense():** Verify that this method works when sorting and when sorting and condensing.
    * **_totalInteractionCount():** Verify that this sums the counts of all interaction types.
    * **_color_for_gap():** Verify that this returns the correct color for hydrogen bonds and for
    a sample of values.
    * **_color_for_atom_class():** Verify that this returns the correct color for a sample of values.
* @todo

**@todo Not yet tested:**
* **annularDots()** and related functions.
* **overlapScale** parameter to DotScorer::check_dot(); this impacts the visual appearance of spikes but no scores.

## Regression testing

Regression tests were performed between Probe2 and the original Probe code in March-April 2022 before the release
of Probe2.  These were facilitated by the https://github.com/ReliaSolve/new_probe_regression repository.
The following sets of command-line arguments were used to compare the two when run with different comparisons.
These are listed for the file 1bti_reduced but similar for others:
- Comparing dot scores and extra atom information:
    - `probe -quiet -kin -mc -self all -count -sepworse -DUMPATOMS outputs/1bti_reduced.orig.dump 1bti_reduced.pdb > outputs/1bti_reduced.old.out`
    - `mmtbx.probe2 source_selection=all approach=self count_dots=True output.separate_worse_clashes=True output.file_name=outputs/1bti_reduced.new.out output.dump_file_name=outputs/1bti_reduced.new.dump 1bti_reduced.pdb`
    - The two atom-dump files were compared to ensure that all atoms were present in both files,
    that they were within 0.05 Angstroms of each other, and that they matched in radius and in acceptor, donor,
    and metallic status.
    - The 'tot' lines from the old and new outputs were compared to ensure that they were identical.
    (**Note:** Probe2 does not separately report Car contacts; they are included in the C reports.)
- Comparing dot distributions for -self case:
    - `probe -quiet -DUMPH2O -kin -mc -self all -sepworse 1bti_reduced.pdb > outputs/1bti_reduced.old.out`
    - `mmtbx.probe2 source_selection=all output.add_kinemage_keyword=True record_added_hydrogens=True approach=self output.separate_worse_clashes=True output.file_name=outputs/1bti_reduced.new.out 1bti_reduced.pdb`
    - Both the original and new output files were opened in King, with each set of dots turned on and off to see
    if there were visible missing or moved dots between the two cases.  This includes looking at the
    Phantom Hydrogens added in each case.
- Comparing dot distributions for -surface case:
    - `probe -quiet -kin -mc -outside all 1bti_reduced.pdb > outputs/1bti_reduced.old.surface.kin`
    - `mmtbx.probe2 source_selection=all output.add_kinemage_keyword=True approach=surface output.file_name=outputs/1bti_reduced.new.surface.kin 1bti_reduced.pdb`
    - Both the original and new output files were opened in King, with each set of dots turned on and off to see
    if there were visible missing or moved dots between the two cases.

The following files were tested and passed the above tests.  Each input file was first run through Reduce2 to produce
the corresponding _reduced.pdb file.
- 1bti (tests alternates): Matches.
- 1xso (tests ions): Probe2 properly identifies a Phantom Hydrogen as bonded to its Oxygen in one instance
where Probe did not because it was too close to be considered bonded, resulting in a slight difference in the
surface representation.  Other outputs match.
- 3wrp (tests too-close waters): Probe2 properly identifies clashes between Water Oxygens and heavy atoms,
possibly because Probe considers them to be bonded.
- **@todo:** Other molecules and snippets.

- **@todo:** Test -once and -both on specific portions of a specific structure.
- **@todo:** Test "not water" selection on 3wrp to ensure that selections are working.

**@todo**

