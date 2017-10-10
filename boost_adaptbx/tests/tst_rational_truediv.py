from __future__ import division
from past.builtins import execfile
import libtbx.load_env
import os
execfile(libtbx.env.under_dist("boost_adaptbx",
                               os.path.join("tests", "tst_rational.py")))
