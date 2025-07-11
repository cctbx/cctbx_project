"""
Script for copying or symbolically linking CCTBX builds into $PREFIX

The script assumes that the build directory has a complete build

By default, the installation will be placed into the sys.prefix
directory of the calling python.

Usage: libtbx.python install_modules.py
"""
from __future__ import absolute_import, division, print_function

import argparse
import glob
import os
import sys

from shutil import copy, copytree, ignore_patterns, rmtree
from subprocess import check_output

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
  if not os.path.exists(src):
    print('  {src} does not exist, skipping'.format(src=src))
  elif os.path.exists(dst):
    print('  {dst} already exists'.format(dst=dst))
  else:
    print('  source:      ' + src)
    print('  destination: ' + dst)
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
def fix_rpath(src):
  """
  Fix relative paths for unix systems

  This is not necessary with conda-build.

  Parameters
  ----------
    src: str
      The source path

  Returns
  -------
    Nothing
  """

  platform = sys.platform
  if platform == 'win32':
    return
  if platform.startswith('linux'):
    libraries = []
    env = dict(os.environ)
    if env.get('LD_LIBRARY_PATH', None) is not None:
      del env['LD_LIBRARY_PATH']
    output = check_output(['ldd', src], env=env).decode('utf8').split('\n')
    for line in output:
      print(line)
      if 'not found' in line:
        libraries.append(line.split()[0])
    print(libraries)

# =============================================================================
def copy_build(env, prefix=None, ext_dir=None, sp_dir=None, copy_egg=False, link=False):
  """
  Copies the following items,
    1) Binaries from build/exe_dev and build/*/exe to $PREFIX/bin
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
    sp_dir: str
      The destination site-packages directory. This is only used if copy_egg is True
    copy_egg: bool
      If True, egg-info directories are copied to the site-packages directory.
      The sp_dir parameter needs to be set.
    link: bool
      If True, symbolic links are used instead of copying.

  Returns
  -------
    Nothing
  """

  # ---------------------------------------------------------------------------
  def loop_copy(src_path, dst_path, name, filenames):
    """
    Convenience function for looping over files to copy
    """
    cmd = 'Copying'
    if link:
      cmd = 'Linking'
    print(cmd + ' ' + name)
    print('-'*79)
    for src_file in filenames:
      src = os.path.join(src_path, src_file)
      if not os.path.exists(src):
        print('  {src} does not exist'.format(src=src))
        continue
      dst = os.path.join(dst_path, src_file)
      copy_cmd(src, dst, link)
    print('Done')
    print()
  # ---------------------------------------------------------------------------

  old_cwd = os.getcwd()

  # binaries and headers
  # may need to add rpath fixes for binaries when necessary
  # ---------------------------------------------------------------------------
  for name, name_dir in [('binaries', 'exe_dev'), ('headers', 'include')]:
    src_path = os.path.join(abs(env.build_path), name_dir)
    if name_dir == 'exe_dev':
      name_dir = 'bin'
    dst_path = os.path.join(prefix, name_dir)
    if os.path.isdir(src_path):
      filenames = os.listdir(src_path)
      loop_copy(src_path, dst_path, name, filenames)
  src_path = abs(env.build_path)
  dst_path = os.path.join(prefix, 'bin')
  os.chdir(src_path)
  all_names = glob.iglob('*/exe/*')
  module_names = set()
  for name in all_names:
    split_name = name.split(os.sep)
    module_names.add(split_name[0])
  for module_name in module_names:
    src_path = os.path.join(abs(env.build_path), module_name, 'exe')
    filenames = os.listdir(src_path)
    loop_copy(src_path, dst_path, '{} binaries'.format(module_name), filenames)

  # libraries
  # ---------------------------------------------------------------------------
  src_path = os.path.join(abs(env.build_path), 'lib')
  dst_path = os.path.join(prefix, 'lib')
  os.chdir(src_path)
  all_names = glob.glob('lib*') \
              + [f for f in glob.glob('*.lib') if not f.endswith('_ext.lib')]
  lib_names = []
  for name in all_names:
    if name.endswith('egg-info'):
      continue
    lib_names.append(name)
  loop_copy(src_path, dst_path, 'libraries', lib_names)

  # Python extensions
  # ---------------------------------------------------------------------------
  src_path = os.path.join(abs(env.build_path), 'lib')
  dst_path = ext_dir
  all_names = glob.iglob('*ext.*')
  ext_names = []
  for name in all_names:
    if name.endswith('egg-info'):
      continue
    ext_names.append(name)
  loop_copy(src_path, dst_path, 'Python extensions', ext_names)

  # .egg-info directories
  # ---------------------------------------------------------------------------
  if copy_egg:
    if sp_dir is None:
      raise RuntimeError('''\
The site-packages directory ("sp_dir") parameter must be set for copying egg directories.''')
    src_path = os.path.join(abs(env.build_path), 'lib')
    dst_path = sp_dir
    egg_names = glob.iglob('libtbx*egg-info')
    loop_copy(src_path, dst_path, 'egg-info directories', egg_names)

  # extra build stuff
  # ---------------------------------------------------------------------------
  share_dir = os.path.join(prefix, 'share', 'cctbx')
  if not os.path.exists(share_dir):
    os.mkdir(share_dir)
  directory_names = []
  for module in env.module_list:
    directory_names += module.names_active()
  loop_copy(abs(env.build_path), share_dir, 'extra build directories', directory_names)

  os.chdir(old_cwd)

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
          # skip boost
          if module == 'boost' and src.endswith('boost'):
            continue
          # copy subdirectories for some modules
          elif module in ['phenix', 'phaser_voyager']:
            if module == 'phenix':
              src_root = dist_path
              subdirs = ['phenix', 'wxGUI2']
            elif module == 'phaser_voyager':
              src_root = dist_path / 'src'
              subdirs = ['New_Voyager', 'Voyager']
            for subdir in subdirs:
              src = abs(src_root / subdir)
              dst = os.path.join(sp_dir, subdir)
              copy_cmd(src, dst, link)
              # phenix/wxGUI2 is also an expected path (fix upstream)
              if subdir == 'wxGUI2':
                dst = os.path.join(sp_dir, 'phenix', subdir)
                copy_cmd(src, dst, link)
            continue
          elif module in ['elbow', 'phaser', 'phasertng', 'PyQuante', 'reel']:
            if module == 'elbow':
              # yacc.py and resources
              for m in ['yacc.py', 'resources']:
                src = abs(dist_path / m)
                dst = os.path.join(sp_dir, m)
                copy_cmd(src, dst, link)
            src = abs(dist_path / module)
          dst = os.path.join(sp_dir, os.path.basename(src))
          copy_cmd(src, dst, link)
          if module == 'tntbx':
            for f in ['__init__.py', 'eigensystem.py']:
              src = abs(dist_path / module / f)
              dst = os.path.join(sp_dir, module, f)
              copy_cmd(src, dst, link)
    except KeyError:
      print(dist_path)

    print('Done')
    print()

# =============================================================================
def remove_build(env, prefix=None, ext_dir=None, sp_dir=None):
  """
  Remove configured modules from site-packages directory

  Parameters
  ----------
    env: libtbx.env_config.environment
      The libtbx environment
    prefix: str
      The destination $PREFIX directory
    ext_dir: str
      The destination directory for Python extensions
    sp_dir: str
      The destination site-pacakges directory

  Returns
  -------
    Nothing
  """

  # ---------------------------------------------------------------------------
  def loop_remove(dst_path, name, filenames):
    """
    Convenience function for looping over files to remove
    """
    print('Removing ' + name)
    print('-'*79)
    for src in filenames:
      src = os.path.join(dst_path, src)
      if os.path.exists(src):
        print('  source: ' + src)
        remove_cmd(src)
      else:
        print('  {src} not found.'.format(src=src))
    print('Done')
    print()
  # ---------------------------------------------------------------------------

  old_cwd = os.getcwd()

  # binaries and headers
  # ---------------------------------------------------------------------------
  for name, name_dir in [('binaries', 'exe_dev'), ('headers', 'include')]:
    src_path = os.path.join(abs(env.build_path), name_dir)
    if name_dir == 'exe_dev':
      name_dir = 'bin'
    dst_path = os.path.join(prefix, name_dir)
    filenames = os.listdir(src_path)
    loop_remove(dst_path, name, filenames)
  src_path = abs(env.build_path)
  dst_path = os.path.join(prefix, 'bin')
  os.chdir(src_path)
  all_names = glob.iglob('*/exe/*')
  module_names = set()
  for name in all_names:
    split_name = name.split(os.sep)
    module_names.add(split_name[0])
  for module_name in module_names:
    src_path = os.path.join(abs(env.build_path), module_name, 'exe')
    filenames = os.listdir(src_path)
    loop_remove(dst_path, '{} binaries'.format(module_name), filenames)

  # libraries
  # ---------------------------------------------------------------------------
  src_path = os.path.join(abs(env.build_path), 'lib')
  dst_path = os.path.join(prefix, 'lib')
  os.chdir(src_path)
  all_names = glob.glob('lib*') \
              + [f for f in glob.glob('*.lib') if not f.endswith('_ext.lib')]
  lib_names = []
  for name in all_names:
    lib_names.append(name)
  loop_remove(dst_path, 'libraries', lib_names)

  # extensions
  # ---------------------------------------------------------------------------
  dst_path = ext_dir
  all_names = glob.iglob('*ext.*')
  ext_names = []
  for name in all_names:
    ext_names.append(name)
  loop_remove(dst_path, 'Python extensions', ext_names)

  # .egg-info directories
  # ---------------------------------------------------------------------------
  dst_path = sp_dir
  egg_names = glob.iglob('libtbx*egg-info')
  loop_remove(dst_path, 'egg-info directories', egg_names)

  # extra build stuff
  # ---------------------------------------------------------------------------
  share_dir = os.path.join(prefix, 'share', 'cctbx')
  directory_names = []
  for module in env.module_list:
    directory_names += module.names_active()
  loop_remove(share_dir, 'extra build directories', directory_names)

  os.chdir(old_cwd)

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

    # extra stuff
    for module_name, extra_stuff in [
      ('elbow', 'yacc.py'),
      ('elbow', 'resources'),
      ('phenix', 'wxGUI2'),
      ('phaser_voyager', 'New_Voyager'),
      ('phaser_voyager', 'Voyager'),
      ]:
      if module == module_name:
        src = os.path.join(sp_dir, extra_stuff)
        if os.path.exists(src):
          print('  ' + src)
          remove_cmd(src)

    print('Done')
    print()

# =============================================================================
def run():
  parser = argparse.ArgumentParser(description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)

  default_sys_prefix = sys.prefix
  if sys.platform == 'darwin' and 'python.app' in default_sys_prefix:
    default_sys_prefix = default_sys_prefix.split('python.app')[0]
  default_sys_prefix = os.path.abspath(default_sys_prefix)
  default_sp_dir = None
  default_lib_dynload_dir = None
  for p in sys.path:
    if default_sp_dir is None \
      and p.startswith(default_sys_prefix) \
      and p.endswith('site-packages'):
      default_sp_dir = p
    if default_lib_dynload_dir is None \
      and p.startswith(default_sys_prefix) \
      and p.endswith('lib-dynload'):
      default_lib_dynload_dir = p

  parser.add_argument(
    '--prefix', default=default_sys_prefix, type=str,
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
  parser.add_argument('--preserve-egg-dir', '--preserve_egg_dir', action='store_true',
    help="""When set, the egg and egg-info directories are copied as well.""")
  parser.add_argument(
    '--fix-rpath', '--fix_rpath', action='store_true',
    help="""When set, the relative paths are fixed for library and extension
      files. NOT IMPLEMENTED YET""")
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
    remove_build(env, prefix=namespace.prefix, ext_dir=namespace.ext_dir,
                 sp_dir=namespace.sp_dir)
    remove_modules(env, sp_dir=namespace.sp_dir)
  else:
    copy_build(env, prefix=namespace.prefix, ext_dir=namespace.ext_dir,
               sp_dir=namespace.sp_dir, copy_egg=namespace.preserve_egg_dir,
               link=namespace.link)
    copy_modules(env, sp_dir=namespace.sp_dir, link=namespace.link)

  if namespace.fix_rpath:
    raise NotImplementedError

  return 0

# =============================================================================
if __name__ == '__main__':
  sys.exit(run())
