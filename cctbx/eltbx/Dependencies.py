# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.files = (
      "global/error.cpp",
      "eltbx/basic.cpp",
      "eltbx/tinypse.cpp",
      "eltbx/tinypsemodule.cpp",
      "eltbx/icsd_radii.cpp",
      "eltbx/icsd_radiimodule.cpp",
      "eltbx/wavelengths.cpp",
      "eltbx/wavelengthsmodule.cpp",
      "eltbx/caasf_it1992.cpp",
      "eltbx/caasf_it1992module.cpp",
      "eltbx/caasf_wk1995.cpp",
      "eltbx/caasf_wk1995module.cpp",
      "eltbx/efpfdp.cpp",
      "eltbx/fpfdpmodule.cpp",
      "eltbx/henke.cpp",
      "eltbx/henkemodule.cpp",
      "eltbx/sasaki.cpp",
      "eltbx/sasakimodule.cpp",
      "eltbx/neutron.cpp",
      "eltbx/neutronmodule.cpp",
      "eltbx/tst_tinypse.py",
      "eltbx/tst_icsd_radii.py",
      "eltbx/tst_wavelengths.py",
      "eltbx/tst_caasf_it1992.py",
      "eltbx/tst_caasf_wk1995.py",
      "eltbx/tst_henke.py",
      "eltbx/tst_sasaki.py",
      "eltbx/tst_neutron.py",
    )

    lib = (
      "error",
      "basic",
      "tinypse",
      "icsd_radii",
      "wavelengths",
      "caasf_it1992",
      "caasf_wk1995",
      "efpfdp",
      "henke",
      "sasaki",
      "neutron",
    )

    self.libraries = {
      "eltbx": lib,
    }

    self.boost_python_modules = {
      "tinypse": ("tinypsemodule", "tinypse", "basic", "error"),
      "icsd_radii": ("icsd_radiimodule", "icsd_radii", "basic", "error"),
      "wavelengths": ("wavelengthsmodule", "wavelengths", "basic", "error"),
      "caasf_it1992": ("caasf_it1992module", "caasf_it1992", "basic", "error"),
      "caasf_wk1995": ("caasf_wk1995module", "caasf_wk1995", "basic", "error"),
      "fpfdp": ("fpfdpmodule", "efpfdp", "basic", "error"),
      "henke": ("henkemodule", "henke", "efpfdp", "basic", "error"),
      "sasaki": ("sasakimodule", "sasaki", "efpfdp", "basic", "error"),
      "neutron": ("neutronmodule", "neutron", "basic", "error"),
    }

  def make_test(self):
    print r"""
tst:
	python tst_tinypse.py
	python tst_icsd_radii.py
	python tst_wavelengths.py
	python tst_caasf_it1992.py
	python tst_caasf_wk1995.py
	python tst_henke.py
	python tst_sasaki.py
	python tst_neutron.py
"""
