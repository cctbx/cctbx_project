PRIME: Post-RefInement and MErging program for XFEL
***Quick guide***
      prime.run your_input_file.phil
  To obtain the input file run:
      prime.run > your_input_file.phil
  You can modify parameters inside the .phil file to match your experiments.
  Parameters most used include:
      - data = /path/to/your/integration/results (you can specify multiple lines of data)
      - d_min, d_max are the minimum and maximum resolutions used to post-refine and merge
        your dataset.
***Update 2016-07-25***
  NOTE MAJOR UPGRADE: For advance usage, the program is now divided into three seperate command lines.
  prime.genref
      Generates a reference set (if needed) and prepare (initial scaling) integration pickles
  prime.merge
      Merge for an mtz file. The mtz files are stored in run_no/mtz with new file marked with
      incremental number. This commands allow you to do merging without re-postrefine
      the data pickles. This way you can change resolution cut-off or I/sigI filter for
      merging only. See .phil merge parameters set.
  prime.postrefine
      From prime.genref and prime.merge results, prime.postrefine will perform one round of
      post-refinement. You'll need to run prime.merge to obtain the mtz file and merging statistics.
  All commands work on the same .phil format as used in prime.run.
***Tutorials***
  Visit Cctbx.xfel wiki page: http://viper.lbl.gov/cctbx.xfel/index.php/Cctbx.prime for detail of
  program inputs and running test cases.
***Citation***
  Uervirojnangkoorn et al. (2015). Enabling X-ray Free Electron Laser Crystallography for Challenging
  Biological Systems from a Limited Number of Crystals. https://elifesciences.org/content/4/e05421.
