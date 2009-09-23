from __future__ import division
import libtbx.load_env
import os
execfile(libtbx.env.under_dist("boost_adaptbx",
                               os.path.join("tests", "tst_rational.py")))
