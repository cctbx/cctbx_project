# Computational Crystallography Toolbox

[![Code Issues](https://www.quantifiedcode.com/api/v1/project/ab81909bad9d437b842c37a4da4f3c53/badge.svg)](https://www.quantifiedcode.com/app/project/ab81909bad9d437b842c37a4da4f3c53)

## Introduction

The Computational Crystallography Toolbox (cctbx) is being developed as the open source component of the PHENIX system. The goal of the PHENIX project is to advance automation of macromolecular structure determination. PHENIX depends on the cctbx, but not vice versa. This hierarchical approach enforces a clean design as a reusable library. The cctbx is therefore also useful for small-molecule crystallography and even general scientific applications.

The cctbx also provides some of the key component of the Olex 2 software. Olex 2 is dedicated to the workflow of small molecule crystallographic studies. It features a powerful and flexible refinement engine, olex2.refine, which is developed as part of the cctbx,
in the smtbx top-module.

To maximize reusability and, maybe even more importantly, to give individual developers a notion of privacy, the cctbx is organized as a set of smaller modules. This is very much like a village (the cctbx project) with individual houses (modules) for each family (groups of developers, of any size including one).

The cctbx code base is available without restrictions and free of charge to all interested developers, both academic and commercial. The entire community is invited to actively participate in the development of the code base. A sophisticated technical infrastructure that enables community based software development is provided by GitHub. This service is also free of charge and open to the entire world.

The cctbx is designed with an open and flexible architecture to promote extendability and easy incorporation into other software environments. The package is organized as a set of ISO C++ classes with Python bindings. This organization combines the computational efficiency of a strongly typed compiled language with the convenience and flexibility of a dynamically typed scripting language in a strikingly uniform and very maintainable way.

Use of the Python interfaces is highly recommended, but optional. The cctbx can also be used purely as a C++ class library.

## Installation

The easiest way to set up a development environment from scratch is to:

1. Download https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py in the directory where the cctbx and its dependencies shall be installed
2. 
  - On Linux or Mac OS execute it: `python bootstrap.py` (you may want to run it with the `--help` option first to discover the available options).
  - On Windows follow the instructions detailed on https://github.com/cctbx/cctbx_project/wiki/How-to-build-CCTBX-on-Windows.

The installation will take a long while but the script will verbosely describe what it does.
