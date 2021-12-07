# Contributing to the Computational Crystallography Toolbox (cctbx)

## cctbx Developer Guidance

##### Table of Contents  
- [Short version](#shortversion) 
- [Long version](#longversion) 
- [Guidance for Developing Tests](#guidancetests) 
  - [Location](#testlocation)
  - [Name](#testname)
  - [Run Time](#testruntime)
  - [3 steps for adding steps](#teststeps)
  - [Guidelines for writing tests](#testguide)
  - [Some useful cctbx tools for tests](#testtools)
  - [How to use models and data](#testmodelsdata)
  - [Using non-cctbx modules](#testnoncctbxmodules)
  - [How to handle errors](#testerrors)
  - [Miscellaneous](#testmixed)
  - [Example for a test](#testexample)
  - [Good etiquette](#testetiquette)

<a name="shortversion"/>

## Short version

As discussed at the Gordon Research Conference for Diffraction Methods in Structural Biology, 07/30/18, the following guidelines are agreed upon for developers wishing to contribute code to cctbx:

* cctbx is open source.
* We welcome third party contributions, but note that there is a cost for managing third party contributions.
* cctbx is not a free-for-all.
* Contributors should be encouraged to review their plans for contributions before doing a large amount of work and submitting a pull request, but that shouldn't preclude users from submitting requests.
* Interface changes and other breaking changes should be communicated to the bulletin boards (cctbxbb@phenix-online.org for example).

**Expectations for contributions**
* Include documentation: docstrings for functions, in-line comments, and help strings for command line programs.
* Include tests.  Testing includes unit tests (individual functions) and integration tests (full program), and possibly multiple of each with different inputs.
* Run the cctbx tests.  If appropriate, run tests from other projects affected by the new code.
* Do not add external dependencies without discussion.

**Expectations for developers**
* We will respond to contributions, at least with "I don't have time to review this right now".
* We can deny pull requests, but if we do so, we should do so with constructive feedback, following consistent rules.
* The build meeting team will review tickets and pull requests to ensure regular feedback.

<a name="longversion"/>

## Long version


The Computational Crystallography Toolbox (cctbx) is a large code base under active development by several groups worldwide. There are more than 1 million lines of code, 500 commits per month, and 20 active developers. It is therefore very important that all contributors follow some basic guidelines. Keep in mind that the intention is for the cctbx to provide a fully featured code base for crystallographic calculations while also remaining lightweight and straightforward to compile and install.

### 1. No new dependencies unless absolutely necessary

  * It is important to keep the cctbx easy to compile and install. The inclusion of third-party packages that have their own dependencies is therefore strongly discouraged. This applies especially if a dependency on a new compiler language is introduced.

  * Check if the desired functionality already exists in the cctbx.

  * Ask other developers if the functionality exists already in the code base.

  * If a new dependency is needed, **discuss with the other developers** before adding it.

### 2. Avoid code duplication

  * It is often the case that the desired functionality already exists somewhere in the cctbx. Therefore:
    
    * Review the code base prior to implementing new functionality.
    
    * Ask other cctbx developers if they know that this functionality exists already
    
  * Add new context independent functions into appropriate modules, not into specialized code. Example: a place for the function that calculates distance between two points is scitbx; however, it may be very tempting to inline this function into a specialized code as needed thus creating code duplication.

  * Use constants from a central place. Example: use `math.pi` instead of defining
`pi = 3.14` every time it is needed. Note: there are plenty of other constants available, such as `scitbx::constants::two_pi_sq`; add more as needed. (Make sure not to use an OS-dependent constant)

### 3. Use the cctbx coding style

  * The Python Style guide is generally useful
Note however that where cctbx deviates from PEP8, follow cctbx (for example using 2 spaces to indent instead of 4).

  * Constructs should be used to make subsequent use and testing as easy to debug as possible (i.e. standard out and standard error output)

  * Printing output

    * Use `print(bla, file=log)` and not `print(bla)`.

    * Use `show` function to print a summary or result of code execution instead of in-lining print statements directly into the code.

    * Avoid unconditional printing.

  * Consistency. If editing an existing file, follow the code style of that file.

  * Code clarity. Ideally, clearly written code does not need documentation (however, even clearly written code may benefit from documentation!). Having this in mind:

    * Use self-explicable function/class/argument/file names. A clear long name is better than a short cryptic name.

    * While valid exceptions exist, generally if a function does not fit a page then it is likely to benefit from being broken into smaller functions.

  * Common modern computation courses emphasize the importance of in-code documentation

    * Since what is clear to one person might not so clear to another (or to the same person after several months or years).

    * Auto documentation creation (like sphinx) use the in-code documentation strings.

    * Time saving: a developer can read the documentation to understand what a function does, what are the arguments’ types and formats and what it returns. This is much better than having to read and figure out what a function needs and what it does.

    * Longevity: Since developers and collaborators change over time and since the cctbx and related systems such as Phenix and DIALS are quite complex software, uncommented code might be a barrier for continuous development.

    * Update documentation to reflect code changes

  * Developers are encouraged to correct severe deviations from coding styles found anywhere in codebase.

  * Developers are encouraged to modify code comments when unclear, outdated or in a format that does not render well in the SPHINX automated documentation.

### 4. Run find clutter before commits

  * `libtbx.find_clutter` primarily checks for a few common errors:

    * No `from __future__ import division`. Use `libtbx.add_from_future_import_division` to fix. Important: always use the `//` operator to indicate integer division; the `/` operator is exclusively for floating point division.

    * Tabs or trailing whitespace: use `libtbx.clean_clutter` to fix.

    * Unused imports: use `libtbx.find_unused_imports_crude` to find them and remove them.

    * The use of bare except statements is prohibited, since it causes the `try...except` to catch `KeyboardInterrupt`. At a minimum use `except Exception`.

  * Before submitting your code, be sure to test it again after fixing problems by running `libtbx.find_clutter`.

### 5. Help and guidance

  * Use the cctbx mailing list to ask for guidance from other developers, and locate specific features in the current code base.

  * The other cctbx developers are an invaluable resource that can be used to help with getting started in cctbx development.

  * Link to the mailing list: http://www.phenix-online.org/mailman/listinfo/cctbxbb

### 6. Send an email to the mailing list stating the intent to submit a new code tree within the cctbx.

  * There are several sub projects within the cctbx, often embodied in the form of code trees within the main cctbx project. New code trees should not be introduced without good reason, and if abandoned the code tree should be removed. This reduces the accumulation of distracting code in the cctbx over time.

  * Periodically, unused code trees may be removed from the cctbx to minimize clutter.

### 7. Subscribe to git check-in alerts sent via email (code changes alerts) and review the diffs.

  * This will minimize surprises later when someone changes someone’s code.

  * Any inefficiencies or bugs spotted in diffs should be pointed out to a respective contributor or fixed by anyone who found them first.

### 8. Python do’s and don’ts:

  * Use inheritance to specialize classes whenever possible to avoid the duplication of code.

  * Put most imports inside the method whenever possible, thus avoiding imports within the global space of a module. Nested imports create huge runtime overhead, particularly as the code base has grown so large over time.

  * Never use `import *`, thus it should always be clear within a module where a name comes from, i.e., `from math import sin, cos, pi`. The only exception is within the `__init__.py` module of a package, where it is permissible to import all items from a C++ extension module.

  * Never use `isinstance()`: a method should not be forced to inquire what type an argument is in order to know how to perform. Instead, the method is entitled to expect arguments to conform to an interface specification; for example if the method prints an object, it should be expected (or documented) that the object should have a `__str__()` method. More details of this discussion may be found at http://www.canonical.org/~kragen/isinstance

### 9. Tests (see more below)

  * Any newly added functionality requires a unit test.

  * Developers are encouraged to add unit tests to any existing functionality found to have no unit test.

  * Any bug fix requires a regression test.

### 10. Miscellaneous:
  
  * Some IDEs have Code inspection tools and Style formatting tools, that can help maintain the style and avoid other code pitfalls (for example: pycharm)

<a name="guidancetests"/>

## Guidance for Developing Tests

There are several reasons for writing good tests:
- **Tests ensure that the code does what it was intended to do:** During development, tests give us confidence that the code still works whenever new functionality is added. For established library code, tests can show us if code modifications -  in lower-level libraries or during refactoring -  yield errors. It may also happen that code intended to fix a bug actually produces conflicts.
- **Tests exemplify how to use your code:** Given the limited number of developers, we don’t have time to write extensive documentation. Therefore, unit tests can be considered as minimal examples that illustrate how to use the code. 
- **Tests may help debugging:** It can be helpful to ask a user to run certain tests (or all tests of a module) on their system.

*Note*: not all tests available in the code base are good examples to follow. When adding a new test please follow the guidelines below. Ask questions on the cctbx mailing list.

<a name="testlocation"/>

### 1. Location:
  * Each module has a directory called regression. This is where tests should be placed. For example: `mmtbx/regression/`

  * Each module has a file called `run_tests.py`, for example `mmtbx/run_tests.py`. This file runs all tests listed in it. If you forget to add your test in that file, it won't be executed.

  * *Note*: Historically, many test files were added next to the actual implementation (in the same directory). Those should eventually be moved to the regression directory.

<a name="testname"/>

### 2. Name:

  * The general name template for a test file is `tst_xxx.py`, where "xxx" may be the name of functionality or file being tested. Example: `tst_miller.py`.

<a name="testruntime"/>

### 3. Run time (per test file `tst_xxx.py`):

  * The faster, the better. Generally execution time should be well under 30 seconds, in exceptional cases 60 seconds should be the absolute max.

  * Slow tests may be marked as such by adding them to a list named `tst_list_slow` in the `run_tests.py` file.  These tests will not be executed by default, but can be run by adding the flag slow_tests=True to libtbx.run_tests_parallel. Slow tests should only be added with good reason, such as if 1) they are designed to test large problems that will always run slowly, such as constructing and inverting large matrices or 2) they are designed to reproduce results from a publication.

<a name="teststeps"/>

### 4. Adding a test involves three steps:

  * Create a file `tst_xxx.py`

  * Place it into the `module/regression` directory

  * Edit `module/run_tests.py` file (otherwise the test will not be executed).

<a name="testguide"/>

### 5. Guidelines for writing tests:

  * Tests should be focused, clear and exercise one functionality at a time. A general template is shown in the inset below.
  
  * Think about tests **during development**. If you start writing tests once the tool is finished, then you most likely forget important corner cases.

  * One file may contain several tests.

  * Ideally, the test code - "exercise" in the example below - should not exceed a page in length, so it can be quickly read and understood. This makes it also easier to fix.

  * A docstring should state the test objective, means and expected result. Docstrings should be put in triple quotes.

  * If a test fails, the failure should be obvious by showing the full traceback (no swallowing tracebacks with printing "TEST FAILED").
 
   * Output actual values in a failed assert statement:

```assert a>b, "%f > %f assertion failed" % (a,b)```

<a name="testtools"/>

### 6. Some useful cctbx tools for tests:

  * `libtbx.test_utils` has many useful tools for tests. Use them as much as possible and add more as needed. Example: use `approx_equal` to assert that a number computed in the test is close to the expected result.

  * Auxiliary functionalities are available in specialized locations such as `cctbx/development`. An example is `cctbx/development/random_structure.py`, which generates a random atomic model. Don't hesitate to add more abstracted functionalities if they are used repeatedly across multiple tests and deemed useful for future tests.

  * Use `libtbx.easy_run.call`. Print the command that is about to run to standard output. Avoid  `libtbx.easy_run.fully_buffered` because it hides the output and the traceback. 

<a name="testmodelsdata"/>

### 7. How to use models and data

  * Avoid storing input files with models and data as much as possible as they take up too much space. Instead, use inputs that can be generated at run time.
  
  * Re-use models and data that are already available. For example, to use the model 1jyp, you can do it in two ways:
      * code example from `modules/phenix_regression/validation/tst_phi_psi_2.py`:
        ```
        from mmtbx.regression import model_1yjp
        ...
        def main():
          dm = DataManager()
          dm.process_model_str('testing', model_1yjp)
          model = dm.get_model()
        ```
      * code example from `modules/phenix_regression/command_line/test_rosetta_refine.py`:
        ```
        def exercise() :
          pdb_file = libtbx.env.find_in_repositories(
            relative_path="phenix_regression/pdb/1yjp_h.pdb",
            test=os.path.isfile)
          mtz_file = libtbx.env.find_in_repositories(
            relative_path="phenix_regression/reflection_files/1yjp.mtz",
            test=os.path.isfile)
          assert (not None in [pdb_file, mtz_file])
        ```
  * Many tests don't require a full length protein. Consider using fragments of a model that exemplify the property being tested. Save these fragment(s) as PDB records or mmCIF string in the body of the test file `tst_xxx.py`. The strings can be read as pdb_input object like this: `pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)`

  * If you put a short description of the property being tested in the REMARK section of a PDB record, then it immediatly clear what the fragment is used for.
  
  * You can use `random_structure` to generate a random model.
  
  * Diffraction data can be calculated from an atomic model, for example with phenix.fmodel()

<a name="testnoncctbxmodules"/>

### 8. Using non-cctbx modules

We test the cctbx independently of the other trees to make sure that it will function properly for non-Phenix developers. In practice, we may use imports of proprietary modules (such as phenix, elbow, phaser, solve_resolve) but these are always made "inline", i.e. imported as needed rather than at the top of the file, and they're usually optional extensions to modules that have independent functionality.

How to check if a module is available:

```
import libtbx.load_env
# Run test only if phenix module is available
if (not libtbx.env.has_module("phenix")) :
  print("phenix tree missing, skipping test")
else :
  # run test
  tst_000()
```
<a name="testerrors"/>

### 9. How to handle errors

  * Don't use standard error, as it won't be shown
  
  * It is OK to raise exceptions
  
  * Don't rely too much on parsing output files. Access results programatically.

<a name="testmixed"/>

### 10. Miscellaneous:

  * It is best to keep the structure and style of tests as similar as possible, so that anyone (and not only the test author) can fix a broken test if necessary. Remember, fixing broken tests is not a pleasant exercise and often is time consuming. Therefore, any test design that can make this task easier is greatly appreciated; one is keeping tests similar in structure and style.

  * Tests should be robust w.r.t. platform, compiler and rounding errors.

  * Broken tests stop others from committing their code. **Fixing a failed test is the highest priority.**

<a name="testexample"/>

### 11. Example for a test:

```
from __future__ import absolute_import, division, print_function
from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal

def exercise():
  """ Make sure 2*2 is 4. """
  x=2.
  result=x*x
  assert approx_equal(result, 4.)

if(__name__ == "__main__"):
  exercise()
  print(format_cpu_times())
  print("OK")
```

<a name="testetiquette"/>

### 12. Good etiquette:

  * Everyone should run `t96` after changing code and before committing.
  
  * Run only one instance of `t96` on anaconda at once 

  * If your changes result in failed tests: Review your code and fix the tests. Even if they have been written originally by other developers. If you break code, it is your responsibility to fix it.

  * Fix tests in a timely manner.

  * Don't cheat: For example don't raise the standard deviation without understanding why the values differ more than before.

  * Don't comment out tests to prevent them from being run
