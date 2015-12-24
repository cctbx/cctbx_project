from os import path
import os
import subprocess
import shutil
import zipfile

import libtbx.load_env

licence_text = """\
This CCTBX build make use of OpenBLAS and of its dependency libgfortran.
We reproduce below the licence of both.

OpenBLAS
--------
%(openblas_licence)s

libgfortran
-----------
Libgfortran is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.

Libgfortran is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Under Section 7 of GPL version 3, you are granted additional
permissions described in the GCC Runtime Library Exception, version
3.1, as published by the Free Software Foundation.

You should have received a copy of the GNU General Public License and
a copy of the GCC Runtime Library Exception along with this program;
see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
<http://www.gnu.org/licenses/>.  */
"""

def run(build, stage, install, package, platform_info, bits):
  # Build a distro which can optimally run on any machine
  # with Intel or AMD processors
  if build:
    subprocess.check_call(['make',
                           'USE_THREAD=1',
                           'NUM_THREADS=16',
                           'DYNAMIC_ARCH=1'] +
                           (['BINARY=%i' % bits] if bits else []) +
                           ['NO_STATIC=1'])

  # Stage it one level up from the current build directory
  stage_dir = path.join(abs(libtbx.env.build_path.dirname()), 'openblas')
  if stage:
    subprocess.check_call(['make',
                           'PREFIX=%s' % stage_dir,
                           'install'])

  # Install the headers and the DLL's in the CCTBX build directory
  # Note that we need to install the runtime library for GNU Fortran and GCC
  # We also install a README
  if install:
    openblas_inc = abs(libtbx.env.include_path / 'openblas')
    if path.isdir(openblas_inc): shutil.rmtree(openblas_inc)
    shutil.copytree(path.join(stage_dir, 'include'), openblas_inc)
    if platform_info['platform'].startswith('mingw'):
      shutil.copy(path.join(stage_dir, 'bin', 'libopenblas.dll'),
                  abs(libtbx.env.lib_path))
      shutil.copy(path.join(stage_dir, 'lib', 'libopenblas.dll.a'),
                  path.join(abs(libtbx.env.lib_path), 'openblas.lib'))
      for dll in ('libgfortran-3.dll', 'libquadmath-0.dll',
                  'libgcc_s_dw2-1.dll'):
        shutil.copy(path.join('c:/mingw/bin', dll), abs(libtbx.env.lib_path))
    for f in ('COPYING3', 'COPYING.RUNTIME'):
      if platform_info['platform'].startswith('mingw'):
        fmt = {'filename':f}
        fmt.update(platform_info)
        shutil.copy(
          'c:/mingw/share/doc/gcc/%(gcc_version)s/%(filename)s' % fmt,
          abs(libtbx.env.build_path))
    licences = {
      'openblas_licence': open('LICENSE').read(),
    }
    with open(path.join(abs(libtbx.env.build_path),
                        'openblas_licence'), 'w') as license:
      license.write(licence_text % licences)

  # Package the files added to the CCTBX build directory
  if package:
    if platform_info['platform'].startswith('mingw'):
      archive = zipfile.ZipFile(
        abs(libtbx.env.build_path.dirname() / 'openblas.zip'), mode="w")
      try:
        openblas_inc = libtbx.env.include_path / 'openblas'
        for p in ([openblas_inc] +
                  [openblas_inc / f for f in os.listdir(abs(openblas_inc))] +
                  [libtbx.env.lib_path / lib
                   for lib in ('libopenblas.dll', 'openblas.lib',
                               'libgfortran-3.dll', 'libquadmath-0.dll',
                               'libgcc_s_dw2-1.dll')] +
                  [libtbx.env.build_path / f
                   for f in ('copying3', 'copying.runtime',
                             'openblas_licence')]):
          archive.write(
            abs(p),
            arcname=path.relpath(abs(p), abs(libtbx.env.build_path)))
      finally:
        archive.close()

if __name__ == '__main__':
  import argparse
  import sys

  # Gather platform information and check we support it
  platform_info = {
    'gcc_version': subprocess.check_output(['gcc', '-dumpversion']),
    'platform': subprocess.check_output(['gcc', '-dumpmachine'])
  }
  platform_info = dict((k,v.strip()) for k,v in platform_info.iteritems())
  supported_platforms = ('mingw32', 'mingw64')
  if platform_info['platform'] not in supported_platforms:
    print ("*** Only the following platforms are supported: " +
           ",".join(supported_platforms) + ' ***')
    sys.exit(1)

  # Parse arguments
  p = argparse.ArgumentParser(
    description=('Build OpenBLAS and prepare a package for distribution.\n'
                 'Please run this script from the top of an OpenBLAS working '
                 'directory, within MSYS shell on Windows. You will '
                 'need GNU C++ and Fortran compiler installed, using '
                 'MinGW to do so on Windows.')
  )
  features = ('build', 'stage', 'install', 'package')
  for arg in features:
    p.add_argument('--%s' % arg, dest=arg, action='store_true')
    p.add_argument('--no-%s' % arg, dest=arg, action='store_false')
  p.set_defaults(**dict((arg, False) for arg in features))
  p.add_argument('--bits', type=int, choices=(None, 32, 64), default=None,
                 help='Whether to build a 32- or 64-bit library '
                      '(None means that OpenBLAS build system shall decide)')
  args = vars(p.parse_args())

  # Run
  try:
    run(platform_info=platform_info, **args)
  except subprocess.CalledProcessError, e:
    print "\n*** Error %i ***\n" % e.returncode
    print "--- Reminder ---\n"
    print p.usage
