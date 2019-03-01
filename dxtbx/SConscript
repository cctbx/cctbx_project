import sys
import os
import libtbx.load_env
Import("env_etc")

env_etc.dxtbx_dist = libtbx.env.dist_path("dxtbx")
env_etc.dxtbx_include = os.path.dirname(env_etc.dxtbx_dist)
env_etc.dxtbx_includes = []
env_etc.dxtbx_common_includes = [env_etc.base_include,
                                 env_etc.libtbx_include,
                                 env_etc.scitbx_include,
                                 env_etc.boost_adaptbx_include,
                                 env_etc.boost_include,
                                 env_etc.dxtbx_include] + env_etc.cbflib_common_includes
env_etc.dxtbx_libs = ["tiff", "cbf"]
env_etc.dxtbx_hdf5_libs = ["hdf5"]
env_etc.dxtbx_lib_paths = [env_etc.base_lib, env_etc.libtbx_lib, os.path.join(sys.prefix, 'lib')]
env_etc.dxtbx_hdf5_lib_paths = []
if (sys.platform == "win32" and env_etc.compiler == "win32_cl"):
  env_etc.dxtbx_libs = ["libtiff", "cbf", "boost_python"]
  env_etc.dxtbx_hdf5_libs = ["libhdf5"]
  env_etc.dxtbx_includes.append(libtbx.env.under_base(os.path.join('HDF5-1.8.16', 'include')))
  env_etc.dxtbx_lib_paths = [
                              env_etc.libpath_python,
                              env_etc.libtbx_lib,
                              libtbx.env.under_base('libtiff'),
                            ]
  env_etc.dxtbx_hdf5_lib_paths = [
                                   libtbx.env.under_base(os.path.join('HDF5-1.8.16', 'lib'))
                                 ]

  env_etc.dxtbx_includes.append(os.path.join(env_etc.cctbx_include,"msvc9.0_include"))
  env_etc.dxtbx_includes.append(libtbx.env.under_base('libtiff'))

  if (libtbx.env.build_options.use_conda):
    env_etc.dxtbx_includes.extend(env_etc.conda_cpppath)
    env_etc.dxtbx_lib_paths.extend(env_etc.conda_libpath)
    # library changes
    # tiff.lib instead of libtiff.lib for newer libtiff conda packages
    env_etc.dxtbx_libs = ['tiff' if x == 'libtiff' else x
                          for x in env_etc.dxtbx_libs]
    # add zlib.lib for hdf5
    env_etc.dxtbx_hdf5_libs.append('zlib')

# for the hdf5.h file - look at where Python is coming from unless is OS X
# framework build... messy but appears to work on Linux and OS X
include_root = os.path.split(env_etc.python_include)[0]
if 'Python.framework' in include_root:
  include_root = os.path.join(
    include_root.split('Python.framework')[0], 'include')
if os.path.exists(os.path.join(include_root, 'hdf5.h')):
  env_etc.dxtbx_includes.append(include_root)
else:
  # check for PSDM installation. Example:
  # /reg/g/psdm/sw/external/hdf5/1.8.6/x86_64-rhel5-gcc41-opt/include
  psdm_hdf5_path = os.path.join(os.environ.get('SIT_ROOT',""),
                                'sw', 'external', 'hdf5', '1.8.6',
                                os.environ.get('SIT_ARCH',""), 'include')
  if os.path.exists(psdm_hdf5_path):
    env_etc.dxtbx_common_includes.append(psdm_hdf5_path)
  psdm_hdf5_path = os.path.join(os.environ.get('SIT_ROOT',""),
                                'sw', 'external', 'hdf5', '1.8.6',
                                os.environ.get('SIT_ARCH',""), 'lib')
  if os.path.exists(psdm_hdf5_path):
    env_etc.dxtbx_hdf5_lib_paths.append(psdm_hdf5_path)

if (not env_etc.no_boost_python and hasattr(env_etc, "boost_adaptbx_include")):
  Import("env_no_includes_boost_python_ext")
  env = env_no_includes_boost_python_ext.Clone()
  env_etc.enable_more_warnings(env=env)
  env_etc.include_registry.append(
    env=env,
    paths=env_etc.dxtbx_includes + env_etc.dxtbx_common_includes + [env_etc.python_include])

  env.Append(
    LIBS=env_etc.libm + [
    "cctbx",
    "scitbx_boost_python",
    ]+env_etc.dxtbx_libs, LIBPATH=env_etc.dxtbx_lib_paths)

  if env_etc.clang_version:
    wd = ["-Wno-unused-function"]
    env.Append(CCFLAGS=wd)

  env.SharedLibrary(
    target="#lib/dxtbx_ext",
    source=[
      "boost_python/to_ewald_sphere_helpers.cc",
      "boost_python/ext.cpp"],
      LIBS=env_etc.libs_python+env_etc.libm+env_etc.dxtbx_libs, LIBPATH=env_etc.dxtbx_lib_paths)

  nexus = env.SharedLibrary(
    target='#/lib/dxtbx_format_nexus_ext',
    source=[
      'format/boost_python/nexus_ext.cc'],
      LIBS=env_etc.libs_python+env_etc.libm+env_etc.dxtbx_libs+env_etc.dxtbx_hdf5_libs,
      LIBPATH=env_etc.dxtbx_lib_paths+env_etc.dxtbx_hdf5_lib_paths)

  imageset = env.SharedLibrary(
    target='#/lib/dxtbx_imageset_ext',
    source=[
      'boost_python/imageset_ext.cc'],
      LIBS=env_etc.libs_python+env_etc.libm+env_etc.dxtbx_libs+env_etc.dxtbx_hdf5_libs,
      LIBPATH=env_etc.dxtbx_lib_paths+env_etc.dxtbx_hdf5_lib_paths)

  image = env.SharedLibrary(
    target='#/lib/dxtbx_format_image_ext',
    source=[
      'format/boost_python/image_ext.cc'],
      LIBS=env_etc.libs_python+env_etc.libm+env_etc.dxtbx_libs+env_etc.dxtbx_hdf5_libs,
      LIBPATH=env_etc.dxtbx_lib_paths+env_etc.dxtbx_hdf5_lib_paths)

  model = env.SharedLibrary(
    target='#/lib/dxtbx_model_ext',
    source=[
      'model/boost_python/beam.cc',
      'model/boost_python/goniometer.cc',
      'model/boost_python/kappa_goniometer.cc',
      'model/boost_python/multi_axis_goniometer.cc',
      'model/boost_python/panel.cc',
      'model/boost_python/detector.cc',
      'model/boost_python/scan.cc',
      'model/boost_python/scan_helpers.cc',
      'model/boost_python/crystal.cc',
      'model/boost_python/parallax_correction.cc',
      'model/boost_python/pixel_to_millimeter.cc',
      'model/boost_python/experiment.cc',
      'model/boost_python/experiment_list.cc',
      'model/boost_python/model_ext.cc'],
      LIBS=env_etc.libs_python+env_etc.libm+env_etc.dxtbx_libs + env["LIBS"] , LIBPATH=env_etc.dxtbx_lib_paths)
