"""Penalty for f' and f" violating Kramers Kronig relations"""

import os
import pathlib
import numpy as np

from scipy.interpolate import interp1d

INTERP_FUNC = lambda x,y: interp1d(x,y,fill_value="extrapolate")
SAMPLE_DATA_PATH = USER = os.getenv("MODULES") + "/ls49_big_data/data_sherrell"

def parse_data(path, remove_first_line=False):
    """Load and parse input."""
    data_input=pathlib.Path(path).read_text().rstrip()
    lines = data_input.split('\n')
    if remove_first_line:
        start_ind=1
    else:
        start_ind=0
    sf = np.array([[float(p) for p in line.split()] for line in lines[start_ind:]])
    
    """Check energy spacing is constant"""
    energy = sf[:,0]
    denergy = energy[1:]-energy[:-1]
    if np.any(denergy-denergy[0]):
        # interpolate
        sf_0, sf_1 = interpolate(sf[:,0], sf[:,1])
        sf_0, sf_2 = interpolate(sf[:,0], sf[:,2])
        sf = np.stack([sf_0,sf_1,sf_2],axis=1)
    return(sf)

def interpolate(x_original, y_original):
    dx_original = x_original[1:]-x_original[:-1]
    dx = np.min(dx_original)

    x = np.arange(x_original[0],x_original[-1],dx)
    interp = INTERP_FUNC(x_original, y_original)
    y = interp(x) # interpolated f_dp
    return(x,y)
