class PackageStructure:
  name       = "cctbx"
  supporting = () # other structured packages that provide supporting code
  externals  = ("external/boost_python",)
  toolboxes  = ("misc", "uctbx", "sgtbx", "arraytbx",
               "adptbx", "eltbx", "sftbx", "fftbx", "lbfgs", "dmtbx")
  examples   = ("examples/cpp",)
  tbx_subpkg = ("arraytbx","eltbx")
