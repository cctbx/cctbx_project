import libtbx.config
import os

cache = libtbx.config.env()
dist_path = cache.dist_path
build_path = cache.LIBTBX_BUILD

def under_dist(package_name, path):
  return os.path.normpath(os.path.join(
    dist_path(package_name=package_name), path))

def under_build(path):
  return os.path.normpath(os.path.join(build_path, path))

bin_path = under_build("libtbx/bin")
