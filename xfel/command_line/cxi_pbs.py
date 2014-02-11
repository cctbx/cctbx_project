# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.pbs
#
# $Id$

from __future__ import division

import sys


def run(argv=None):
  import libtbx.load_env

  from os import X_OK, access, path
  from libtbx import easy_run

  if argv is None:
    argv = sys.argv

  # Absolute path to the executable Bourne-shell script.
  pbs_sh = path.join(
    libtbx.env.find_in_repositories('cctbx_project/xfel/cxi'),
    'pbs.sh')
  assert access(pbs_sh, X_OK)

  return easy_run.call(' '.join([pbs_sh] + argv[1:]))


if __name__ == '__main__':
  sys.exit(run())
