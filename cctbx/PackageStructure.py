class PackageStructure:
  name       = "cctbx"
  supporting = () # other structured packages that provide supporting code
  externals  = ("external/boost_python",)
  toolboxes  = ("uctbx", "sgtbx", "arraytbx",
               "adptbx", "eltbx", "sftbx", "fftbx", "lbfgs")
  examples   = ("examples/cpp",)
  tbx_subpkg = ("arraytbx","eltbx")
