"""
Script for copying or symbolically linking CCTBX modules into
the site-packages directory.
"""
from __future__ import absolute_import, division, print_function

import argparse
import os
import sys

from shutil import copytree, ignore_patterns, rmtree

# =============================================================================
def copy_modules(sp_dir=None, link=False):
  """
  Copies configured modules to site-packages directory

  Parameters
  ----------
    sp_dir: str
      The destination site-packages directory
    link: bool
      If True, symbolic links are used instead of copying.

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

            if link:
              os.symlink(src, dst)
            else:
              copytree(src, dst, ignore=ignore_patterns('.git*', '.svn*'))
          else:
            print('  Nothing done')
    except KeyError:
      print(dist_path)

    print('Done')
    print()
  return 0

# =============================================================================
def remove_modules(sp_dir=None):
  """
  Remove configured modules from site-packages directory

  Parameters
  ----------
    sp_dir: str
      The destination site-packages directory
    link: bool
      If True, symbolic links are used instead of copying.

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
      if os.path.islink(src):
        os.unlink(src)
        print('  ' + src)
      elif os.path.isdir(src):
        rmtree(src)
        print('  ' + src)
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
    '--sp-dir', '--sp_dir', default=None, type=str, nargs='?',
    const=default_sp_dir,
    help="""The location where the modules will be installed, by default
      the directory is the site-packages location of the calling python.""")
  parser.add_argument(
    '--link', action='store_true',
    help="""When set, instead of copying, symbolic links are created
      for the repository directories.""")
  parser.add_argument(
    '--clean', action='store_true',
    help="""When set, repository directories are removed from the sp_dir.""")

  # show help if no arguments are provided
  if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

  namespace = parser.parse_args()

  if namespace.clean:
    remove_modules(sp_dir=namespace.sp_dir)
  else:
    copy_modules(sp_dir=namespace.sp_dir, link=namespace.link)

  return 0

# =============================================================================
if __name__ == '__main__':
  sys.exit(run())
