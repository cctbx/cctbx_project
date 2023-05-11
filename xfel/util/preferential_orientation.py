from __future__ import division
from dxtbx.model.experiment_list import ExperimentList
from xfel.util.drift import params_from_phil, path_lookup, \
    read_experiments, unique_elements
import numpy as np


message = """
This utility tool aims to determinate, characterise, and quantify the degree
of preferential orientation in crystals. To this aim, it first uses
the Iterative Reweighted Least Squares (IRLS) algorithm to estimate
the concentration matrix A of the Bingham distributio
from all experimental crystal matrices. Then it reports the type and degree
of anisotropy as expressed in the eigenvalues of said matrix.

This is work in progress.
""".strip()


phil_scope_str = """
  scrap {
    input {
      glob = None
        .type = str
        .multiple = True
        .help = glob which matches all expt files to be investigated.
      exclude = None
        .type = str
        .multiple = True
        .help = glob which matches all expt files to be excluded from input.
    }
"""

############################# ORIENTATION STORAGE #############################


############################ ORIENTATION SCRAPPING ############################


########################### ORIENTATION VISUALIZING ###########################


################################ ENTRY POINTS #################################


def run(params_):
    pass


params = []
if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
    print(message)
    exit()
  params = params_from_phil(sys.argv[1:])
  run(params)