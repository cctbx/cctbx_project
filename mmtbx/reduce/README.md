# Reduce2

This is a CCTBX re-implementation of the Richardson Labs' Reduce program from https://github.com/rlabduke/reduce
It is named reduce2 so that it can be built and used alongside the original reduce (which can still be built as
a module) during regression testing.  It is intended that future development and maintenance will happen
on Reduce2 and not on the original code base. The **reduce2** program is a command-line program that
uses the modules in this directory, along with those in mmtbx/hydrogens and mmtbx/probe to place
and optimize Hydrogens on a model.

**Notes:**
* Reduce2 uses the Probe2 module, which requires the chem_data modules from Phenix.

## Python Modules

There are a set of Python modules that implement the core optimization functions.

* **Movers.py:** This contains a set of classes that implement various "Movers", which are
responsible for describing the potential states of sets of atoms, along with atom deletions.
Some handle rotations of hydrogens around a bonded parent and others handle flips of rings,
along with fixups to keep bond lengths correct. They do not provide the motions or deletions,
but merely describe them.

* **InteractionGraph.py:** This contains functions to determine which Movers may interact,
providing checks for all pairs of atoms in a pair of Movers for all of their possible locations.
It is used to determine a graph of Movers that may interact, which is used to determine "Cliques"
that must be jointly optimized and "Singletons" which can be independenly optimized.

* **Optimizers.py:** This contains classes to optimize a set of Movers, including adding them
to a model and determining their optimal states. User code should use the **FastOptimizer** class,
which is derived from other classes that provide more basic, but slower, approaches to optimization
so that the different accelerations can be tested independently.

## Testing

Unit tests are implemented in each Python module to verify that they perform as intended.
When it is run as a script, each file will raise an assertion on failure or exit normally on success.
The tests performed include the following:
* **Movers.py:** Tests the following:
    * **probe.ExtraAtomInfo** is tested to ensure that information in copies can be independently changed.
    * **_MoverNull** base class is tested to verify that it behaves as expected.
    * **_MoverRotator** private class is tested to verify that it behaves as expected.
    * **MoverSingleHydrogenRotator** class is tested to verify that it behaves as expected.
    * **MoverNH3Rotator** class is tested to verify that it behaves as expected.
    * **MoverAromaticMethylRotator** class is tested to verify that it behaves as expected.
    * **MoverTetrahedralMethylRotator** class is tested to verify that it behaves as expected.
    * **MoverAmideFlip** class is tested to verify that it behaves as expected.
    * **MoverHisFlip** class is tested to verify that it behaves as expected.
    * **FixUp()** method is tested to verify that it behaves as expected; it is used by Flips.
* **InteractionGraph.py:** Tests the following:
    * **InteractionGraphAllPairs** function is tested to verify that it behaves as expected.
* **Optimizers.py:** When run with the name of a PDB file, will attempt to optimize it. When run
    with no arguments, tests the following:
    * **AlternatesInModel()** function is tested to verify that it behaves as expected.
    * **GetAtomsForConformer()** function is tested to verify that it behaves as expected.
    * **MoverSingleHydrogenRotator** optimization is tested to verify that it behaves as expected.
    * **MoverNH3Rotator** optimization is tested to verify that it behaves as expected.
    * **MoverHisFlip** optimization is tested to verify that it behaves as expected,
    including locking down the flip when there is a bond with a nearby ion.
    * Optimization of a multi-element clique is tested to ensure it behaves as expected with all
    types of optimizers.
    * The occupancy and B-factor cut-offs for water Oxygens is tested to ensure they are handled.
    * An example model (if run with no command-line arguments) or the specified model file (if
    given on the command line) is optimized to ensure that a complete run can be done. If the
    --dumpAtoms command-line flag was given, an **atomDump.pdb** file will be written with the
    description of the extra atom info for each atom in the model, along with a **deleteme.pdb**
    file that contains the optimized model.

A **tst_reduce.py** script is located in the mmtbx/regression folder within this project.  this
script runs all of the unit tests described above and then several specific test cases are run
to ensure that various corner cases are properly handled.

Version 1.2.2 is the code that passed all of the above tests.

## Regression testing

Regression tests were performed between Reduce2 and the original Reduce code in June-July 2023 before
the release of Reduce2. These were facilitated by the https://github.com/ReliaSolve/new_reduce_regression
repository. The tests include adding flip movers in both cases and then running Probe2 to score the results from
the reduce2 and original reduce outputs and then compares the total scores to see which is better.

The following files were tested, resulting in the same or better scores for most. The ones with worse
scores point to a set of outstanding issues with hydrogen placement that are still being investigated.
- 1bti: Reduce2 better.
- 1crn: Reduce2 better.
- 1ehz: Reduce2 better.
- 1otf: Reduce2 better.
- 1pq7_fragment: Same.
- 1ubq-nh3_test: **Reduce2 worse**.
- 1xso: Reduce2 better.
- 1yk4_snip: Reduce2 better.
- 2mbw_fragment: **Reduce2 worse**.
- 3cp5_fragment: Reduce2 better.
- 4fen: Reduce2 better.
- 5dka_fragment: Reduce2 better.
- 6xhv: Reduce2 better.
- 6t5h_fragment: Reduce2 better.
- 6tte_fragment: **Reduce2 worse**.
- 6vw1_fragment: **Reduce2 worse**.
- 7c31: Reduce2 better.

Additional tests were performed using the **comparison_file** command-line option to compare the
behavior of reduce2 with original reduce on a Mover-by-Mover basis. 3VYK was studied along
with 1DFU. These tests turned up other issues that are being investigated. The 1XSO file was
also studied, which has 164 Movers; the overall score was 75 higher for Reduce2, issues were
created (and then resolved) for individual residues that were significantly lower.

As of 7/21/2023, the issues under investigation are described in the status and issues document at
https://docs.google.com/document/d/1ogzS6QlBnPJRGU2CNzHvw0KhV8E5fXdwO1TwUEM7GQc/edit?usp=sharing 
which is being used to track ongoing issues. The major issues are:
- Some hydrogens are being placed improperly (sometimes too many sometimes too few).
- The restraints information for many files in the PDB is not available, so many files cannot be run.
- The mostly-Python Reduce2 is significantly slower than the original C++ reduce.
