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
    info = None
    bin_files = []
    for package in packages:
      for json_file in os.listdir(meta):
        if json_file.startswith(package) and json_file.endswith('.json'):
          with open(meta/json_file, 'r') as f:
            info = json.load(f)
      if info is not None:
        all_files = info.get('files', [])
        for file in all_files:
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
              with open(new_prefix/bin_name, 'w') as f:
                for line in lines:
                  line = line.strip()
                  if 'dp0' in line:  # make path relative to original location
                    line += '\\..\\Library\\bin'
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
