import libtbx.config
import os

env_cache = libtbx.config.env()
dist_root = env_cache.LIBTBX_DIST_ROOT
dist_path = env_cache.dist_path
build = env_cache.LIBTBX_BUILD

def under_dist_root(path):
  return os.path.normpath(os.path.join(dist_root, path))

def under_dist(package_name, path):
  return os.path.normpath(os.path.join(
    dist_path(package_name=package_name), path))

def under_build(path):
  return os.path.normpath(os.path.join(build, path))
