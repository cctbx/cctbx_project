"""
Script for copying or symbolically linking CCTBX modules into
the site-packages directory.
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
def copy_headers(env, prefix=None, link=False):
  """
  Copy headers from build to $PREFIX/include

  Parameters
  ----------
    env: libtbx.env_config.environment
      The libtbx environment
    prefix: str
      The destination $PREFIX directory
    link: bool
      If True, symbolic links are used instead of copying.

  Returns
  -------
    0 if successful,
  """
  src_path = os.path.join(abs(env.build_path), 'include')
  dst_path = os.path.join(prefix, 'include')
  cmd = 'Copying'
  if link:
    cmd = 'Linking'
  print(cmd + ' headers')
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
    0 if successful,
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
  return 0

# =============================================================================
def remove_headers(env, prefix=None):
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
    0 if successful,
  """
  src_path = os.path.join(abs(env.build_path), 'include')
  dst_path = os.path.join(prefix, 'include')
  print('Removing headers')
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
    0 if successful,
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

  return 0

# =============================================================================
def run():
  parser = argparse.ArgumentParser(description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)

  default_sp_dir = None
  for p in sys.path:
    if p.endswith('site-packages'):
      default_sp_dir = p
      break

  parser.add_argument(
    '--prefix', default=sys.prefix, type=str, nargs='?',
    help="""The $PREFIX location, by default the directory is the sys.path
      location of the calling python.""")
  parser.add_argument(
    '--sp-dir', '--sp_dir', default=default_sp_dir, type=str, nargs='?',
    help="""The location where the modules will be installed, by default
      the directory is the site-packages location of the calling python.""")
  parser.add_argument(
    '--link', action='store_true',
    help="""When set, instead of copying, symbolic links are created
      for the repository directories.""")
  parser.add_argument(
    '--clean', action='store_true',
    help="""When set, repository directories are removed from the sp_dir.""")

  namespace = parser.parse_args()

  # show help if no arguments are provided
  if len(sys.argv) == 1 and namespace.sp_dir is None:
    parser.print_help()
    parser.exit()

  # load original libtbx_env
  if os.getenv('LIBTBX_BUILD') is None:
    raise RuntimeError('''\
The libtbx_env file must be in the original location specified \
by the LIBTBX_BUILD environment variable''')
  import libtbx.load_env
  env = libtbx.env

  # copy or clean
  if namespace.clean:
    remove_headers(env, prefix=namespace.prefix)
    remove_modules(env, sp_dir=namespace.sp_dir)
  else:
    copy_headers(env, prefix=namespace.prefix, link=namespace.link)
    copy_modules(env, sp_dir=namespace.sp_dir, link=namespace.link)

  return 0

# =============================================================================
if __name__ == '__main__':
  sys.exit(run())
