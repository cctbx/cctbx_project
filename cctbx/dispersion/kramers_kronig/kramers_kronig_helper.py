"""Penalty for f' and f" violating Kramers Kronig relations"""

from __future__ import division

import os
import pathlib
import numpy as np
import torch

import libtbx.load_env

from scipy.interpolate import interp1d

INTERP_FUNC = lambda x,y: interp1d(x,y,fill_value="extrapolate")

SAMPLE_DATA_PATH = USER = os.path.join(libtbx.env.find_in_repositories('ls49_big_data'), 'data_sherrell')

def parse_data(path, remove_first_line=False):
    """Load and parse input."""
    data_input=pathlib.Path(path).read_text().rstrip()
    lines = data_input.split('\n')
    if remove_first_line:
        start_ind=1
    else:
        start_ind=0
    sf = np.array([[float(p) for p in line.split()] for line in lines[start_ind:]])

    return(sf)

def interpolate(x_original, y_original, mode="scipy"):
    """
    Interpolate y_original onto a uniform grid, matching the smallest spacing of x_original
    mode can be "scipy" or "torch", torch allows for automatic differentiation to propagate by using PyTorch
    """

    if mode=="scipy":
        numerical_package = np
    elif mode=="torch":
        numerical_package = torch

    dx_original = x_original[1:]-x_original[:-1]
    dx = numerical_package.min(dx_original)

    x = numerical_package.arange(x_original[0],x_original[-1],dx)

    if mode == "scipy":
        interp = INTERP_FUNC(x_original, y_original)
        y = interp(x) # interpolated f_dp
    elif mode == "torch":
        y = interpolate_torch(x_original, y_original, x)
    return(x,y)

def interpolate_torch(x_original, y_original, x_new):
    """Linearly interpolate y_original(x_original), returning y = y_original(x_new), using PyTorch"""

    if torch.abs(x_original[-1]-x_new[-1])<1e-5:
      x_new = x_new[:-1]
      same_end = True
    else:
      same_end = False

    x_new_upper_position = torch.searchsorted(x_original, x_new)
    x_new_lower_position = x_new_upper_position - 1

    x_new_dist_lower_position = x_new - x_original[x_new_lower_position]
    x_new_dist_upper_position = x_original[x_new_upper_position] - x_new

    x_new_dist_lower_position_norm = x_new_dist_upper_position / (x_new_dist_lower_position + x_new_dist_upper_position)
    x_new_dist_upper_position_norm = x_new_dist_lower_position / (x_new_dist_lower_position + x_new_dist_upper_position)

    y_new = y_original[x_new_upper_position]*x_new_dist_upper_position_norm + y_original[x_new_lower_position]*x_new_dist_lower_position_norm

    if same_end:
      y_new = torch.concat((y_new,y_original[-1].expand(1)))

    return(y_new)
