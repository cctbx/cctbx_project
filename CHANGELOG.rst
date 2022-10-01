2022.9
======

* Improvements to cctbx.HKLviewer for viewing output from xtriage or xtricorder
* Improve stability of prepare_map_for_docking
* DataManager models can now have multiple types (e.g. x-ray, neutron, etc.)
* Fix any_reflection_file_reader when reading "hklf+ins/res" (#787)
* Fix binary and source installers when /usr/bin/python is not available (#788)

2022.8
======

* Update to handle new PAE format from AlphaFold
* CCTBXParser can parse DataManager PHIL parameters and provide more
  information about ambiguous parameters
* Initial support for density dependent restraints
* Bug fixes and updates for cctbx.HKLviewer

2022.7
======

* Clean up multiple mmtbx tests
* Updated to RCSB API v2 for accessing data

2022.6
======

* Added --dry-run flag to CCTBXParser to validate inputs
* Bug fixes for cctbx.HKLviewer
* Updates to reduce2 and probe2
* Updates to restraints based on quantum mechanics

2022.5
======

* Added function to return fmodel object in DataManager
* Added options for cubic box and soft masking to resolve_cryo_em
* Updates to tool for preparing maps for docking
* More improvements to cctbx.HKLviewer

2022.4
======

* Added option to keep unmerged data in DataManager
* Bug fixes for cctbx.HKLviewer
* Improvements to finding water in maps

2022.3
======

* Added mmtbx.process_predicted_model command
* Updates and bug fixes to cctbx.HKLviewer
* Added methods for improved handling of heavy hydrogens in model

2022.2
======

* Added --quiet flag to CCTBXParser to suppress output
* Updates to restraints based on quantum mechanics
* Updates and bug fixes to cctbx.HKLviewer

2022.1
======

* Added option to use legacy bulk solvent mask
* Added option to any_reflection_file_reader to control averaging of
  anomalous data columns in MTZ files
* Fix SHELXF formatting where integer values may be interpreted as decimal

2021.12
=======

* Initial support for Python 3.10
* Added support for outputting multi-model mmCIF files
* Adjusted mask gridding for bulk solvent

2021.11
=======

* Added right-handed nucleic acids for DNA
* Improved handling of different unit cells in MTZ file
* Avoid division by zero when rotating 0 degrees
* Added option to ignore secondary strucure annotations when reading
  models through the DataManager
* Initial support for restraints based on quantum mechanics
* Improved consistency of binning by d_star_sq

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
