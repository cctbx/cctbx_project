class PackageStructure:
  name       = "cctbx"
  supporting = () # other structured packages that provide supporting code
  externals  = ("external/boost_python",)
  toolboxes  = ("misc", "uctbx", "sgtbx", "miller", "eltbx", "arraytbx",
               "adptbx", "sftbx", "fftbx", "lbfgs", "dmtbx", "mintbx")
  examples   = ("examples/cpp",)
  tbx_subpkg = ("arraytbx","eltbx")
