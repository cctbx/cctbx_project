"""
Script for copying or symbolically linking CCTBX builds into $PREFIX

The script assumes that the build directory has a complete build

By default, the installation will be placed into the sys.prefix
directory of the calling python.

Usage: libtbx.python install_modules.py
"""
from __future__ import absolute_import, division, print_function

import argparse
import os
import sys

from shutil import copy, copytree, ignore_patterns, rmtree

# =============================================================================
def copy_cmd(src, dst, link):
  """
  Copy or link

  Parameters
  ----------
    src: str
      The source path
    dst: str
      The destination path
    link: bool
      If True, linking is used instead of copying

  Returns
  -------
    Nothing
  """
  if link:
    os.symlink(src, dst)
  else:
    if os.path.isdir(src):
      copytree(src, dst, ignore=ignore_patterns('.git*', '.svn*'))
    else:
      copy(src, dst)

# =============================================================================
def remove_cmd(src):
  """
  Remove file, directory, or unlink

  Parameters
  ----------
    src: str
      The source path

  Returns
  -------
    Nothing
  """
  if os.path.islink(src):
    os.unlink(src)
  elif os.path.isdir(src):
    rmtree(src)
  else:
    os.remove(src)

# =============================================================================
def copy_build(env, prefix=None, ext_dir=None, link=False):
  """
  Copies the following items,
    1) Binaries from build/exe_dev to $PREFIX/bin
    2) Headers from build/include to $PREFIX/include
    3) Libraries from build/lib to $PREFIX/lib
    4) Python extensions from build/lib to $PREFIX/lib/$PYTHON/lib-dynload

  Parameters
  ----------
    env: libtbx.env_config.environment
      The libtbx environment
    prefix: str
      The destination $PREFIX directory
    ext_dir: str
      The destination directory for Python extensions
    link: bool
      If True, symbolic links are used instead of copying.

  Returns
  -------
    Nothing
  """

  # binaries and headers
  # ---------------------------------------------------------------------------
  for name, name_dir in [('binaries', 'exe_dev'), ('headers', 'include')]:
    src_path = os.path.join(abs(env.build_path), name_dir)
    if name_dir == 'exe_dev':
      name_dir = 'bin'
    dst_path = os.path.join(prefix, name_dir)
    cmd = 'Copying'
    if link:
      cmd = 'Linking'
    print(cmd + ' ' + name)
    print('-'*79)
    for src in os.listdir(src_path):
      dst = os.path.join(dst_path, src)
      src = os.path.join(src_path, src)
      print('  source:      ' + src)
      print('  destination: ' + dst)
      copy_cmd(src, dst, link)
    print('Done')
    print()

# =============================================================================
def copy_modules(env, sp_dir=None, link=False):
  """
  Copies configured modules to site-packages directory

  Parameters
  ----------
    env: libtbx.env_config.environment
      The libtbx environment
    sp_dir: str
      The destination site-packages directory
    link: bool
      If True, symbolic links are used instead of copying.

  Returns
  -------
    Nothing
  """

  # copy each module directory
  for module in env.module_dict:
    cmd = 'Copying'
    if link:
      cmd = 'Linking'
    print(cmd + ' ' + module)
    print('-'*79)
    try:
      for dist_path in env.module_dict[module].dist_paths:
        if dist_path is not None:
          src = abs(dist_path)
          dst = os.path.join(sp_dir, os.path.basename(src))
          if os.path.exists(dst):
            print('  {dst} already exists'.format(dst=dst))
            continue
          if os.path.isdir(src):
            print('  source:      ' + src)
            print('  destination: ' + dst)
            copy_cmd(src, dst, link)
          else:
            print('  Nothing done')
    except KeyError:
      print(dist_path)

    print('Done')
    print()

# =============================================================================
def remove_build(env, prefix=None):
  """
  Remove configured modules from site-packages directory

  Parameters
  ----------
    env: libtbx.env_config.environment
      The libtbx environment
    prefix: str
      The destination $PREFIX directory

  Returns
  -------
    Nothing
  """
  for name, name_dir in [('binaries', 'exe_dev'), ('headers', 'include')]:
    src_path = os.path.join(abs(env.build_path), name_dir)
    if name_dir == 'exe_dev':
      name_dir = 'bin'
    dst_path = os.path.join(prefix, name_dir)
    print('Removing ' + name)
    print('-'*79)
    for src in os.listdir(src_path):
      src = os.path.join(dst_path, src)
      if os.path.exists(src):
        print('  source: ' + src)
        remove_cmd(src)
      else:
        print('  {src} not found.'.format(src=src))
    print('Done')
    print()

# =============================================================================
def remove_modules(env, sp_dir=None):
  """
  Remove configured modules from site-packages directory

  Parameters
  ----------
    env: libtbx.env_config.environment
      The libtbx environment
    sp_dir: str
      The destination site-packages directory

  Returns
  -------
    Nothing
  """

  # load original libtbx_env
  if os.getenv('LIBTBX_BUILD') is None:
    raise RuntimeError('''\
The libtbx_env file must be in the original location specified \
by the LIBTBX_BUILD environment variable''')
  import libtbx.load_env
  env = libtbx.env

  for module in env.module_dict:
    src = os.path.join(sp_dir, module)
    cmd = 'Deleting'
    if os.path.islink(src):
      cmd = 'Unlinking'
    print(cmd + ' ' + module)
    print('-'*79)
    for name in [module, module + '_' + env.module_dict[module].mate_suffix]:
      src = os.path.join(sp_dir, name)
      if os.path.exists(src):
        print('  ' + src)
        remove_cmd(src)
      else:
        print('  {name} not found'.format(name=name))

    print('Done')
    print()

# =============================================================================
def run():
  parser = argparse.ArgumentParser(description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)

  default_sp_dir = None
  default_lib_dynload_dir = None
  for p in sys.path:
    if default_sp_dir is None and p.endswith('site-packages'):
      default_sp_dir = p
    if default_lib_dynload_dir is None and p.endswith('lib-dynload'):
      default_lib_dynload_dir = p

  parser.add_argument(
    '--prefix', default=sys.prefix, type=str,
    help="""The $PREFIX location, by default the directory is the sys.prefix
      location of the calling python.""")
  parser.add_argument(
    '--sp-dir', '--sp_dir', default=default_sp_dir, type=str,
    help="""The location where the modules will be installed, by default
      the directory is the site-packages location of the calling python.""")
  parser.add_argument(
    '--ext-dir', '--ext_dir', default=default_lib_dynload_dir, type=str,
    help="""The location where the Python extensions will be installed, by
      default the directory is the lib-dynload location of the calling python.""")
  parser.add_argument(
    '--link', action='store_true',
    help="""When set, instead of copying, symbolic links are created
      for directories and files.""")
  parser.add_argument(
    '--clean', action='store_true',
    help="""When set, directories and files are removed from the $PREFIX.""")

  namespace = parser.parse_args()

  # load original libtbx_env
  if os.getenv('LIBTBX_BUILD') is None:
    raise RuntimeError('''\
The libtbx_env file must be in the original location specified \
by the LIBTBX_BUILD environment variable''')
  import libtbx.load_env
  env = libtbx.env

  # copy or clean
  if namespace.clean:
    remove_build(env, prefix=namespace.prefix)
    remove_modules(env, sp_dir=namespace.sp_dir)
  else:
    copy_build(env, prefix=namespace.prefix, link=namespace.link)
    copy_modules(env, sp_dir=namespace.sp_dir, link=namespace.link)

  return 0

# =============================================================================
if __name__ == '__main__':
  sys.exit(run())
