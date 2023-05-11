from __future__ import division
import glob
import sys
from typing import List

# from dxtbx.model import ExperimentList
# from xfel.util.drift import params_from_phil, read_experiments

import numpy as np
import scipy.spatial.transform
from scipy.linalg import expm, logm


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


############################ ORIENTATION SCRAPPING ############################


# class OrientationScraper:
#   """Class responsible for scraping orientation matrices from expt files"""
#   def __init__(self, parameters) -> None:
#     self.parameters = parameters
#
#   @staticmethod
#   def assemble_orientation_stack(expts: ExperimentList) -> np.ndarray:
#     """Extract N orientation matrices from N expts into a Nx3x3 numpy array"""
#     ori_list = [np.array(expt.crystal.get_U()).reshape(3,3) for expt in expts]
#     return np.stack(arrays=ori_list, axis=0)
#
#   def locate_input_paths(self) -> List:
#     """Return a list of expt paths in scrap.input.glob, but not in exclude"""
#     input_paths, exclude_paths = [], []
#     for ig in self.parameters.scrap.input.glob:
#       input_paths.extend(glob.glob(ig))
#     for ie in self.parameters.scrap.input.exclude:
#       exclude_paths.extend(glob.glob(ie))
#     return [it for it in input_paths if it not in exclude_paths]
#
#   def scrap(self) -> np.ndarray:
#     """Read and return a Nx3x3 orientation matrix based on scrap parameters"""
#     expt_paths = self.locate_input_paths()
#     expts = read_experiments(expt_paths)
#     return self.assemble_orientation_stack(expts)


##################### PREFERENTIAL ORIENTATION CALCULATOR #####################

class PreferentialOrientationCalculator:
  """Calculator of the preferential orientation degree using IRLS algorithm"""
  EPSILON: float = 1e-6
  GUESS: np.ndarray = np.eye(3)
  MAX_ITER: int = 100

  def irls_orientations(self, orientations: np.ndarray) -> np.ndarray:
    """Estimates the average orientation matrix from a dataset of
    3D orientation matrices using IRLS algorithm."""

    # Initial guess for A (arithmetic mean)
    A_old = self.GUESS

    # Initialize weights to identity matrices
    W = np.tile(np.eye(3), (len(orientations), 1, 1))

    for i in range(self.MAX_ITER):
      # Compute residuals
      R = np.array([logm(A_old.T @ ori) for ori in orientations])
      r = np.array([R_i - R_i.T for R_i in R])
      b = np.array([r_i / np.linalg.norm(r_i) if np.linalg.norm(r_i) > 1e-10
                    else np.zeros(3) for r_i in r])
      B = np.array([np.outer(b_i, b_i) for b_i in b])

      # Update the weights matrix W using the Bingham distribution.
      B_inv = np.array([np.linalg.inv(B_i) for B_i in B])
      w = np.array(
        [np.exp(-b_i.T @ B_inv_i @ b_i) for b_i, B_inv_i in zip(b, B_inv)])
      W = np.diag(w)

      # Calculate the Jacobian J and check for convergence.
      J = np.zeros((3, 3))
      for d, w_i in zip(orientations, w):
        J += w_i * d.T @ A_old
      J /= np.sum(w)
      delta = np.linalg.norm(J - np.eye(3))
      if delta < self.EPSILON:
        break

      # Update the orientation matrix.
      A_new = expm((np.linalg.inv(J) @ np.sum(R * W[:, np.newaxis], axis=0)).T)
      A_old = A_new

    return A_new


########################### ORIENTATION VISUALIZING ###########################


################################ ENTRY POINTS #################################


# def run(params_):
#     os = OrientationScraper(parameters=params_)
#     poc = PreferentialOrientationCalculator()
#     ori = os.scrap()
#     avg_ori = poc.irls_orientations(orientations=ori)
#     print(avg_ori)


# params = []
# if __name__ == '__main__':
#   if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
#     print(message)
#     exit()
#   params = params_from_phil(sys.argv[1:])
#   run(params)


def main2():
  from scipy.spatial.transform import Rotation
  oris = Rotation.random(3).as_matrix()
  poc = PreferentialOrientationCalculator()
  avg_ori = poc.irls_orientations(orientations=oris)
  print(avg_ori)

if __name__ == '__main__':
  main2()