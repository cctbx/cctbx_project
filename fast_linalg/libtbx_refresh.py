from __future__ import division
import sys
import fileinput

if sys.platform == 'win32':
  # We need to run fast_blas.build_openblas from the MSYS shell
  # So we need to create a sh-dispatcher on Windows.
  # Unfortunately, method write_bin_sh_dispatcher produces one with a few
  # glitches that need fixing by hand.
  builder = self.env.bin_path / 'fast_linalg.build_openblas'
  self.env.write_bin_sh_dispatcher(
    source_file=
    self.env.under_dist('fast_linalg', 'command_line',
                          return_relocatable_path=True) / 'build_openblas.py',
    target_file=builder)
  for li in fileinput.FileInput([abs(builder)], inplace=True):
    if li.startswith('@REM'):
      li = '#' + li[4:]
    else:
      lpb = 'LIBTBX_PYEXE_BASENAME'
      li = li.replace('c:\\', '/c/').replace('\\', '/')\
        .replace(lpb.lower(), lpb)
    print li,
