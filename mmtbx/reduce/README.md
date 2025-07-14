# Reduce2

This is a CCTBX re-implementation of the Richardson Labs' Reduce program from https://github.com/rlabduke/reduce
It is named reduce2 so that it can be built and used alongside the original reduce (which can still be built as
a module) during regression testing.  It is intended that future development and maintenance will happen
on Reduce2 and not on the original code base. The **reduce2** program is a command-line program that
uses the modules in this directory, along with those in mmtbx/hydrogens and mmtbx/probe to place
and optimize Hydrogens on a model.

**Notes:**
* Reduce2 uses the probe2 module, which requires the chem_data module.
* As of 7/14/2025, reduce2 is switching to new default parameters for the relative weighting of
  external contacts (which remains the same) and both hydrogen bonds and collisions (whic are
  increasing by 10x). This is to work as expected with a radius larger than 0 (it is being switched
  to 0.25, which matches the probe2 default). The original Reduce had switched to a radius of 0.0
  and was not considering external contacts; this fixes that without causing hydrogen bonds to be
  broken. See https://github.com/cctbx/cctbx_project/issues/1072 for details. The defaults for
  probe2 are not being changed, so its default behavior will be different from reduce2 until the
  issue can be fully resolved and both set to the same defaults. The probe radius and relative
  weights can be set in both probe2 and reduce2 using the probe2 Phil parameters if different
  behavior is desired.

**Installing and running:** *Reduce2* is part of the mmtbx module, which is part of the CCTBX
distribution. It can be built and run as part of the CCTBX build process. The CCTBX install and
build processes are described in the README.md file in the root of the CCTBX distribution at
https://github.com/cctbx/cctbx_project/blob/master/README.md. Reduce2 does require the Monomer
library described at https://github.com/cctbx/cctbx_project?tab=readme-ov-file#monomer-library
to be installed. You can download a .conda package from the linked-to releases page and then
after activating your conda environment, run `conda install -c cctbx cctbx_monomer_library` to install it.

Once the system has been installed and the environment configured,
the Reduce2 program can be run from the command line as `mmtbx.reduce2`, and the --help or
--show-defaults options used to see the available options.  Older versions of the program
require you to set output.description_file_name to specify where descriptive text output
will be written, current versions set this by default to add a .txt extension rather than
a .cif or .pdb extension to the output file name.

# C++ Classes

The C++ classes, wrapped for use in Python, make use of CCTBX and Boost structures and define
all of their classes and functions in the molprobity::reduce C++ namespace. These are wrapped
in Python and can be imported as mmtbx_reduce_ext.

* **PositionReturn.*:** Structure to hold atom-behavior information returned from Mover methods.

* **InteractionGraph.*:** Defines the **PairsOverlap** function to determine whether two Movers
have any atoms that interact in any of their possible conformations.

* **Optimizers.*:** Defines the **OptimizerC** class, which determines the best conformations
for a set of Movers in a clique. It caches the scores for atoms to avoid recalculation when all
Movers they interact with are in the same configuration. It splits the clique into smaller
subcliques and optimizes them independently, then combines the results.

* **boost_python/reduce_bpl.cpp:** As is usual for CCTBX projects, this file contains the wrapper
code for the above functions that uses Boost Python to wrap them to be called from Python. It also
defines a wrapper for the **RotatePointDegreesAroundAxisDir** function, which wraps
**scitbx::math::rotate_point_around_axis()** to quickly rotate a point around an axis.

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
to a model and determining their optimal states. User code should use the **Optimizer** class,
which wraps the C++ OptimizerC class.

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
* **PositionReturn.cpp:** Its **PositionReturn_test** function verifies that the C++ structure
is properly wrapped and can be accessed.
* **InteractionGraph.cpp:** Its tests are run from Python, so its **InteractionGraph_test**
test function currently always returns success.
* **Optimizers.cpp:** its **Optimizers_test** function tests the static methods in the class.
The other methods and the class itself is tested from Python.
    * **generateAllStates()** function is tested to verify that it behaves as expected.
    * **nChooseM()** function is tested to verify that it behaves as expected.
    * **subsetGraph()** function is tested to verify that it behaves as expected.
    * **findVertexCut()** function is tested to verify that it behaves as expected.

A **tst_reduce.py** script is located in the mmtbx/regression folder within this project.  this
script runs all of the unit tests described above and then several specific test cases are run
to ensure that various corner cases are properly handled.

A **tst_mmtbx_reduce_ext.py:** Script is located in the mmtbx/reduce folder that runs the test
functions in the C++ classes via their python wrappers. When run, prints nothing on success and
raises an assertion on failure.

## Regression testing

Regression tests were performed between Reduce2 and the original Reduce code in June-July 2023 before
the release of Reduce2. These were facilitated by the https://github.com/ReliaSolve/new_reduce_regression
repository. The tests include adding flip movers in both cases and then running Probe2 to score the results from
the reduce2 and original reduce outputs and then compares the total scores to see which is better.

The following files were tested, resulting in the same or better scores for most. The ones with worse
scores appear to be due to hydrogen placement difference for non-movable hydrogens where Reduce2 is
using more-correct angles that sometimes have lower scores in Probe2 but are nonetheless more accurate.

Additional tests were performed using the **comparison_file** command-line option to compare the
behavior of reduce2 with original reduce on a Mover-by-Mover basis. 3VYK was studied along
with 1DFU. These tests turned up other issues that are being investigated. The 1XSO file was
also studied, which has 164 Movers; the overall score was 75 higher for Reduce2, issues were
created (and then resolved) for individual residues that were significantly lower.

As of 2/29/2024, the issues under investigation are described in the status and issues document at
https://docs.google.com/document/d/1ogzS6QlBnPJRGU2CNzHvw0KhV8E5fXdwO1TwUEM7GQc/edit?usp=sharing 
which is being used to track ongoing issues. The major issues are:
- Ionic bonds are not being handled properly, causing flip Movers to be adjustable near ions.
- The restraints information for many files in the PDB is not available, so many files cannot be run.
- The mostly-Python Reduce2 is significantly slower than the original C++ reduce.
