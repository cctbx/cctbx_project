# LatticeRepLib

`cctbx/uctbx/lrl` contains Python bindings for parts of the LatticeRepLib
package developed by Lawrence C. Andrews and Herbert J. Bernstein.

Currently the LatticeRepLib code is on a fork with a reduced repo size and
some modifications to permit building in a cctbx installation.

## Contact
Daniel Paley
dwpaley@lbl.gov

## Installation

```
$ cd $modules
$ git clone git@github.com:dwpaley/LatticeRepLib
$ pushd LatticeRepLib; git checkout cctbx_build; popd
$ cd ../build
$ source conda_setpaths.sh
$ libtbx.configure LatticeRepLib
$ make
```

## Testing

```
$ cd $modules/cctbx_project/cctbx/uctbx/lrl
$ libtbx.python test_match_lattices.py
```

## Features

### LatticeMatcher

Given a list of lattices, transform each to the setting that matches best to
a reference unit cell. The similarity is measured by distance in the six-
dimensional space of Selling scalars, S6.
