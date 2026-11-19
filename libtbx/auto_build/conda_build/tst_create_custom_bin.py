"""The Windows copies of the dispatchers in <prefix>\\phenix_bin must still
compute LIBTBX_PREFIX as <prefix>\\Library, which the original derives from
its own location in Library\\bin."""
from libtbx.auto_build.conda_build import create_custom_bin

ORIGINAL = [
  "@setlocal",
  "@set LIBTBX_PREFIX=%~dp0",
  "@set LIBTBX_PREFIX=%LIBTBX_PREFIX:~0,-1%",
  '@for %%F in ("%LIBTBX_PREFIX%") do @set LIBTBX_PREFIX=%%~dpF',
  "@set LIBTBX_PREFIX=%LIBTBX_PREFIX:~0,-1%",
  "@set LIBTBX_BUILD=%LIBTBX_PREFIX%\\..\\Library\\share\\cctbx",
  "@set LIBTBX_DISPATCHER_NAME=%~nx0",
]

RELOCATED = [
  "@setlocal",
  '@for %%F in ("%~dp0..\\Library") do @set LIBTBX_PREFIX=%%~fF',
  "@set LIBTBX_BUILD=%LIBTBX_PREFIX%\\..\\Library\\share\\cctbx",
  "@set LIBTBX_DISPATCHER_NAME=%~nx0",
]

def exercise():
  # the four lines deriving LIBTBX_PREFIX from %~dp0 become one that
  # resolves ..\Library from the copy's location
  assert create_custom_bin.relocate_bat_dispatcher(ORIGINAL) == RELOCATED, \
    create_custom_bin.relocate_bat_dispatcher(ORIGINAL)
  # line endings of a file read on any platform do not matter
  assert create_custom_bin.relocate_bat_dispatcher(
    [l + "\r\n" for l in ORIGINAL]) == RELOCATED
  # applying it to a relocated copy changes nothing
  assert create_custom_bin.relocate_bat_dispatcher(RELOCATED) == RELOCATED
  # a script without the block is returned unchanged
  other = ["@echo off", "@set X=%~dp0"]
  assert create_custom_bin.relocate_bat_dispatcher(other) == other

def exercise_copy_bin_matches_the_package_name_exactly():
  # conda-meta holds <name>-<version>-<build>.json; a package whose name
  # merely starts with the requested one (phenix_html or phenix-restraints
  # for phenix, cctbx-base for cctbx) is not the package, and every
  # requested package's file list is used, not the last one
  import json
  import os
  import shutil
  import tempfile
  tmp = tempfile.mkdtemp()
  try:
    os.makedirs(os.path.join(tmp, "conda-meta"))
    os.makedirs(os.path.join(tmp, "bin"))
    for name, files in (("phenix-1.0-0", ["bin/phenix.refine"]),
                        ("phenix_html-1.0-0", ["bin/phenix_html.tool"]),
                        ("phenix-restraints-1.0-0", ["bin/phenix.restraints"]),
                        ("cctbx-1.0-0", ["bin/cctbx.only"]),
                        ("cctbx-base-1.0-0", ["bin/cctbx.python"])):
      with open(os.path.join(tmp, "conda-meta", name + ".json"), "w") as f:
        json.dump({"name": name.rsplit("-", 2)[0], "files": files}, f)
      for file in files:
        open(os.path.join(tmp, file), "w").close()
    create_custom_bin.copy_bin(tmp, "phenix_bin", ["phenix", "cctbx-base"])
    copied = sorted(os.listdir(os.path.join(tmp, "phenix_bin")))
    assert copied == ["cctbx.python", "phenix.refine"], copied
  finally:
    shutil.rmtree(tmp, ignore_errors=True)

if __name__ == "__main__":
  exercise()
  exercise_copy_bin_matches_the_package_name_exactly()
  print("OK")
