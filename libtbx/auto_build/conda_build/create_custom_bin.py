from __future__ import division, print_function
'''
Script for copying package-specific binaries to a custom directory to
avoid cluttering the user's path
'''

import argparse
import json
import os
import sys

from pathlib import Path
from shutil import copy

# =============================================================================
def relocate_bat_dispatcher(lines):
  '''
  Rewrite a copied Windows dispatcher so LIBTBX_PREFIX is still
  <prefix>\\Library.

  The conda-package dispatcher derives LIBTBX_PREFIX from its own location
  (%~dp0 is Library\\bin, so the prefix is its parent). A copy placed in
  <prefix>\\<custom_bin> would compute <prefix> instead, so the lines that
  derive it are replaced by one that resolves ..\\Library relative to the
  copy and normalizes the result.

  Parameters
  ----------
  lines : list of str
      The dispatcher's lines, with or without line endings

  Returns
  -------
  list of str
      The rewritten lines, without line endings
  '''
  out = []
  i = 0
  replaced = False
  derive = ('@set LIBTBX_PREFIX=', '@for %%F in ("%LIBTBX_PREFIX%")')
  while i < len(lines):
    line = lines[i].rstrip('\r\n')
    if not replaced and line.startswith('@set LIBTBX_PREFIX=%~dp0'):
      out.append('@for %%F in ("%~dp0..\\Library") do @set LIBTBX_PREFIX=%%~fF')
      i += 1
      while i < len(lines) and lines[i].rstrip('\r\n').startswith(derive):
        i += 1
      replaced = True
      continue
    out.append(line)
    i += 1
  return out

# =============================================================================
def copy_bin(prefix, custom_bin, packages=[]):
  '''
  Function to copy package-specific binaries to a custom directory to
  avoid cluttering the user's path

  A directory containing the binaries will be created in
  ${prefix}/${custom_bin}

  Parameters
  ----------
  prefix : str or Path
      The installation prefix
  custom_bin : str
      The name of the directory for the binaries you want copied
  packages : list of str
      The names of packages containing the binaries you want copied
  '''
  prefix = Path(os.path.abspath(prefix))
  bin_directories = ['bin', 'Library/bin', 'Scripts']
  meta = prefix/'conda-meta'
  if meta.exists():
    bin_files = []
    for package in packages:
      # <name>-<version>-<build>.json; the name is compared, not the file
      # name prefix, which would also match phenix_html or cctbx-base
      for json_file in os.listdir(meta):
        if json_file.startswith(package) and json_file.endswith('.json'):
          with open(meta/json_file, 'r') as f:
            info = json.load(f)
          if info.get('name') != package:
            continue
          for file in info.get('files', []):
            for bin_dir in bin_directories:
              if file.startswith(bin_dir):
                bin_files.append(Path(file))
    if len(bin_files) > 0:
      new_prefix = prefix/custom_bin
      os.makedirs(new_prefix, exist_ok=True)
      for bin_file in bin_files:
        bin_name = bin_file.name
        if (new_prefix/bin_name).exists():
          print(f'{new_prefix/bin_name} already exists, skipping')
        else:
          if sys.platform == 'win32':
            print(f'Copying {prefix/bin_file} to {new_prefix/bin_name}')
            copy(prefix/bin_file, new_prefix/bin_name)
            if (new_prefix/bin_name).suffix == '.bat':
              print(f'Fixing LIBTBX_PREFIX in {new_prefix/bin_name}')
              with open(new_prefix/bin_name, 'r') as f:
                lines = f.readlines()
              relocated = relocate_bat_dispatcher(lines)
              if relocated == [l.rstrip('\r\n') for l in lines]:
                print(f'Warning: no LIBTBX_PREFIX block found in {bin_name}; '
                      'the copy will compute the wrong prefix')
              with open(new_prefix/bin_name, 'w') as f:
                for line in relocated:
                  f.write(line)
                  f.write('\r\n')
          else:
            print(f'Linking {prefix/bin_file} to {new_prefix/bin_name}')
            os.symlink(prefix/bin_file, new_prefix/bin_name)
          if (new_prefix/bin_name).exists():
            print(f'  {new_prefix/bin_name} created')
          else:
            print(f'Warning: {new_prefix/bin_name} was not created')
    else:
      print('No binary files were found, skipping copy')
  else:
    print('The conda-meta directory does not exist. Cannot find package json.')

# =============================================================================
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('--prefix', type=str,
    help='''The prefix of the installation. The conda-meta directory must exist in this directory.''',
    required=True)
  parser.add_argument('--custom_bin', type=str,
    help='''The name of the custom binary directory. It will exist in the prefix directory.''',
    required=True)
  parser.add_argument('--packages', type=str, nargs='*',
    help='''The packages to search.''',
    required=True)

  namespace = parser.parse_args()

  copy_bin(namespace.prefix, namespace.custom_bin.lower(), namespace.packages)

# =============================================================================
# end
