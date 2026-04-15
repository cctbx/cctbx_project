# Computational Crystallography Toolbox
[![Build Status](https://dev.azure.com/cctbx/cctbx_project/_apis/build/status/Updates/Update%20build%20cache?branchName=master)](https://dev.azure.com/cctbx/cctbx_project/_build/latest?definitionId=8&branchName=master) [![Conda Version](https://img.shields.io/conda/vn/conda-forge/cctbx-base.svg)](https://anaconda.org/conda-forge/cctbx-base) [![Nightly conda package tests](https://github.com/cctbx/cctbx/actions/workflows/nightly_tests.yml/badge.svg?event=schedule)](https://github.com/cctbx/cctbx/actions/workflows/nightly_tests.yml) [![Conda Platforms](https://anaconda.org/conda-forge/cctbx-base/badges/platforms.svg)](https://anaconda.org/conda-forge/cctbx-base) [![DOI](https://img.shields.io/badge/DOI-10.1107/S0021889801017824-blue.svg)](https://doi.org/10.1107/S0021889801017824)

##### Table of Contents
- [Introduction](#introduction)
- [Installation](#install)
- [Building a development version](#developmentversion)
- [Contributing to the cctbx](#contributing)
- [Nightly builds and tests](#nightlybuilds)

<a name="introduction"/>

## Introduction

The Computational Crystallography Toolbox (cctbx) is being developed as the open source component of the [Phenix project](https://phenix-online.org). The goal of the Phenix project is to advance automation of macromolecular structure determination. Phenix depends on the cctbx, but not vice versa. This hierarchical approach enforces a clean design as a reusable library. The cctbx is therefore also useful for small-molecule crystallography and even general scientific applications.

The cctbx also provides some of the key component of the Olex 2 software. Olex 2 is dedicated to the workflow of small molecule crystallographic studies. It features a powerful and flexible refinement engine, olex2.refine, which is developed as part of the cctbx,
in the smtbx top-module.

To maximize reusability and, maybe even more importantly, to give individual developers a notion of privacy, the cctbx is organized as a set of smaller modules. This is very much like a village (the cctbx project) with individual houses (modules) for each family (groups of developers, of any size including one).

The cctbx code base is available without restrictions and free of charge to all interested developers, both academic and commercial. The entire community is invited to actively participate in the development of the code base. A sophisticated technical infrastructure that enables community based software development is provided by GitHub. This service is also free of charge and open to the entire world.

The cctbx is designed with an open and flexible architecture to promote extendability and easy incorporation into other software environments. The package is organized as a set of ISO C++ classes with Python bindings. This organization combines the computational efficiency of a strongly typed compiled language with the convenience and flexibility of a dynamically typed scripting language in a strikingly uniform and very maintainable way.

Use of the Python interfaces is highly recommended, but optional. The cctbx can also be used purely as a C++ class library.


<a name="install"/>

## Installation

The easiest way to install cctbx is through the [Conda package manager](https://docs.conda.io/en/latest/). We recommend the [Miniforge installers](https://github.com/conda-forge/miniforge?tab=readme-ov-file#miniforge3) since they provide a minimal environment and default to the `conda-forge` channel.

There are two packages available, `cctbx` and `cctbx-base`. The `cctbx` package is `cctbx-base` with some additional GUI packages (e.g. `wxpython`, `pyside2`, `ipython`).

With the `conda` command available, a new `cctbx-base` environment named `my_env` can be created with
```
conda create -n my_env -c conda-forge cctbx-base
```
To choose a specific version of Python, add the `python` package with the specific version
```
conda create -n my_env -c conda-forge cctbx-base python=3.11
```
Then the environment can be activated with
```
conda activate my_env
```

To install `cctbx-base` into the currently active environment, use
```
conda install -c conda-forge cctbx-base
```
The `python` package with a specific version can be added to change the version of `python` that is already installed in the active environment.

There are also regular python packages for `cctbx-base` and `cctbx`. They can
be installed with
```
pip install cctbx-base
```

### Monomer library

Some programs in cctbx require information about geometric restraints for molecules. This information is available in the `chem_data` conda package from our [releases](https://github.com/cctbx/cctbx_project/releases). Download the `chem_data` conda package and install in your active environment with
```
conda install <chem_data package>
```
Alternatively, the latest version of the `chem_data` conda package is available
on the `chem_data` channel from https://anaconda.org. You can install it with
```
conda install -c chem_data chem_data
```
Due to space limitations on anaconda.org, we can only keep the latest version
available.
The `chem_data` package is built from the [`chem_data`](https://gitlab.com/phenix_project/chem_data) and [`geostd`](https://github.com/phenix-project/geostd) repositories.

<a name="developmentversion"/>

## Building a development version

1. Download https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py in the directory where the cctbx and its dependencies shall be installed
2. Run `python bootstrap.py` (you may want to run it with the `--help` option first to discover the available options).
  - For better compatibility with newer operating systems, `conda` packages can be used for dependencies. Add the `--use-conda` flag and the command becomes `python bootstrap.py --use-conda`. This will run the `miniconda` installer if `conda` cannot be found. The environment with the dependencies will be located in the `conda_base` directory. See the description of the `--use-conda` flag from the `--help` output for more details.

The installation will take a long while but the script will verbosely describe what it does.

<a name="contributing"/>

## Contributing to the cctbx

For a more detailed description on how to contribute to the cctbx please visit our [contribution guide](https://github.com/cctbx/cctbx_project/blob/master/CONTRIBUTING.md).

<a name="nightlybuilds"/>

## Nightly builds and tests
[![Build Status](https://dev.azure.com/cctbx-release/feedstock-builds/_apis/build/status/nightly-feedstock?branchName=main)](https://dev.azure.com/cctbx-release/feedstock-builds/_build/latest?definitionId=11&branchName=main) [![Conda Version](https://img.shields.io/conda/vn/cctbx-nightly/cctbx-base.svg)](https://anaconda.org/cctbx-nightly/cctbx-base) [![Conda Platforms](https://anaconda.org/cctbx-nightly/cctbx-base/badges/platforms.svg)](https://anaconda.org/cctbx-nightly/cctbx-base)

A nightly build of the `conda` packages are available on the [`cctbx-nightly` channel](https://anaconda.org/cctbx-nightly/repo). To use these packages, prepend `-c cctbx-nightly` as a channel to the commands above. For example, the command to create a new `my_env` environment would become
<pre><code>
conda create -n my_env <b>-c cctbx-nightly</b> -c conda-forge cctbx-base
</code></pre>
This will use the `cctbx-base` package from the `cctbx-nightly` channel, but pull the remaining dependencies from `conda-forge`.

There are also nightly builds for the regular python packages. They can be
installed with
```
pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple cctbx-base
```
The `--extra-index-url` flag may be needed if `pip` is not able to find the
other dependencies on TestPyPI.

Nightly builds are only updated if there are additional commits from the previous build.

A subset of tests is run on the current `cctbx-base` packages from the `conda-forge` and `cctbx-nightly` channels every night (10 pm Pacific) to test compatibility with the latest packages from `conda-forge`. Additional source files for `fable` and `antlr3` are needed for the tests. The nightly test details can be viewed by clicking the "Nightly conda package tests" badge near the beginning of this README.
