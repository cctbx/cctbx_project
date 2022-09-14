'''
Script for creating the files for setting up a user's environemnt after
running a conda-based installer.
'''
from __future__ import absolute_import, division, print_function

import argparse
import os
import sys

# =============================================================================
# unix files
bash_template = '''\
#!/usr/bin/env bash

export {program}_PREFIX={prefix}
export {program}_VERSION={version}
export PATH={bin_dir}:$PATH
'''

csh_template = '''\
#!/usr/bin/env csh

setenv {program}_PREFIX {prefix}
setenv {program}_VERSION {version}
setenv PATH {bin_dir}:$PATH
'''

# Windows file
bat_template = '''\
@echo off
set {program}_PREFIX=%~dp0
set {program}_VERSION={version}
set PATH={bin_dir};%PATH%
'''

# =============================================================================
def write_files(program=None, prefix=None, bin_dir=None, version=None,
  destination=None):
  '''
  Populate template with arguments

  Parameters
  ----------
  program : str
      see parser arguments in run()
  prefix : str
      see parser arguments in run()
  bin_dir : str
      see parser arguments in run()
  version : str
      see parser arguments in run()
  destination : str
      see parser arguments in run()
  '''
  if version is None:
    version = ''

  program = program.upper()
  filename = '{program}_env'.format(program=program.lower())
  filename = os.path.join(destination, filename)
  if sys.platform == 'win32':
    with open('.'.join([filename, 'bat']), 'w') as f:
      f.write(bat_template.format(
        program=program,
        prefix=prefix,
        bin_dir=bin_dir,
        version=version
      ))
  else:
    for template, ext in zip((bash_template, csh_template), ('sh', 'csh')):
      with open('.'.join([filename, ext]), 'w') as f:
        f.write(template.format(
          program=program,
          prefix=prefix,
          bin_dir=bin_dir,
          version=version
        ))

# -----------------------------------------------------------------------------
def run():
  parser = argparse.ArgumentParser(description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)

  default_prefix = sys.prefix
  if sys.platform == 'win32':
    default_bin = os.path.join(default_prefix, 'Library', 'bin')
  default_bin = os.path.join(default_prefix, 'bin')

  parser.add_argument(
    '--program', default='cctbx', type=str,
    help='''The program name for constructing the filenames. For example,
      the filename for bash will be <program>_env.sh. (cctbx)'''
  )
  parser.add_argument(
    '--prefix', default=default_prefix, type=str,
    help='''The $PREFIX location, by default the directory is the sys.prefix
      location of the calling python. ({default_prefix})'''.format(default_prefix=default_prefix)
  )
  parser.add_argument(
    '--bin-dir', default=default_bin, type=str,
    help='''The location to be added to $PATH, by default it is $PREFIX/bin.
      If the argument is a relative path, it will be appended to the
      --prefix argument. ({default_bin})'''.format(default_bin=default_bin)
  )
  parser.add_argument(
    '--version', default=None, type=str,
    help='''The version to be set (no default)'''
  )
  parser.add_argument(
    '--destination', default=os.getcwd(), type=str,
    help='''The directory to write the output ({destination})'''.format(destination=os.getcwd())
  )

  namespace = parser.parse_args()

  bin_dir = namespace.bin_dir
  if not os.path.isabs(bin_dir):
    bin_dir = os.path.join(namespace.prefix, namespace.bin_dir)
    bin_dir = os.path.abspath(bin_dir)

  prefix = namespace.prefix
  if not os.path.isabs(prefix):
    prefix = os.path.abspath(prefix)

  try:
    write_files(
      program=namespace.program,
      prefix=prefix,
      bin_dir=bin_dir,
      version=namespace.version,
      destination=namespace.destination
    )
  except IOError:
    print('''
There was an issue with writing the files. Please make sure you have
permissions to write to {destination} and that there is enough disk space.
'''.format(destination=namespace.destination))
    raise

# =============================================================================
if __name__ == '__main__':
  sys.exit(run())
