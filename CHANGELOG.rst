2021.10
=======

* Initial migration of MolProbity functionality (probe and reduce) to mmtbx
* Initial tool for likelihood-based map preparation for docking
* Improvements to ADP refinement for real-space

2021.9
======

* Improved structure factor calculation at ultra-low resolution
* Improved processing of prediced models
* Added diffBragg to simtbx for modeling pixels in stills to improve structure factors
* Added suitename to mmtbx for classifying RNA

2021.8
======

* Added tools for processing predicted models based on the error estimate
* Updated list of modified amino and nucleic acids
* Better handling of sequence files with empty sequences

2021.7
======

* Fix pickling error with anomalous_probability_plot
* Fix bug in reading data CIF file with paired data and sigma arrays of
  different sizes
* Added functions for retrieving H-bond types and Van der Waals radii to
  the model manager class
* Sequence validation will only use protein or nucleic acid residues for
  alignment

2021.6
======

* More improvements to bulk solvent masking for multiple regions
* Updates to ensemble refinement
* Enable conversion of some numpy types into flex types instead of
  requiring that the types match (e.g. int to float is now supported)

2021.5
======

* Improved bulk solvent masking with support for multiple regions

2021.4
======

* Improved parsing of reciprocal space data in CIF
* CCTBXParser can handle intermixed arguments for Python >= 3.7
* Consolidate management of conda depenencies with conda-devenv

2021.3
======

* Initial support for native compilation on Apple Silicon
* Real-space refinement of occupancies and isotropic ADP
* Improvements in map_model_manager

  * Split up map and model by NCS groups
  * Create new map_model_manager with resampled maps

2021.2
======

* Improved remediator code for converting PDB version 2 format to version 3
* Add compilation support for Boost 1.72 and 1.74

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
