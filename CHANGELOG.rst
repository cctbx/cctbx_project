2021.1
======

* Improvements to cctbx.HKLviewer for displaying reciprocal space data

2020.12
=======

* BIOMT/MTRIX matrices in model reading

  * Added option to loosen handling of improper matrices in DataManager
  * Make behavior conistent between mmCIF and PDB formats

* Improvements to map_model_manager

  * Better handling of cases when information is missing
  * Calculate the RMSD of matching residues between models

2020.11
=======

* Updated API for fetching data from RCSB

2020.10
=======

* Added basic ``flex`` arrays for fixed width integer types (`#533 <https://github.com/cctbx/cctbx_project/pull/533>`_)

  * Signed types (``int8``, ``int16``, ``int32``, ``int64``)
  * Unsigned types (``uint8``, ``uint16``, ``uint32``, ``uint64``)
  * Additional functions may be wrapped in the future to support these types

* Improved building of downstream software with ``cctbx`` conda package

  * In some cases, the location of ``annlib`` is not found properly

2020.8
======

* First release on GitHub and conda-forge
