"""Penalty for f' and f" violating Kramers Kronig relations"""

import pathlib
import numpy as np

from scipy.signal import hilbert
from scipy.interpolate import CubicSpline,interp1d
from scipy.signal.windows import get_window

INTERP_FUNC = lambda x,y: interp1d(x,y,fill_value="extrapolate")


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


def hilbert_transform(f_dp):
    return(np.imag(hilbert(f_dp, N=None, axis=-1)))


def hilbert_transform_sherrell(f_dp):
    """Version of hilbert transform implemented in the Sherrell (2014)"""
    Hz = np.fft.fft(f_dp)
    H = np.zeros_like(Hz)
    H[0]=0
    H[1:(len(H)//2)]=-1
    H[(len(H)//2):]=1
    H = H*1j
    dF = np.fft.ifft(1j*H*Hz+Hz) #.reshape((1, npts))
    return(dF.imag)


def get_f_p(energy, f_dp, padn=5000,
            trim=0,
            Z = 26, # atomic number
            include_Z_term=False,
            hilbert_transform_func=hilbert_transform,
            window_type='cosine',
            known_response_energy=None,
            known_response_f_p=None,
            known_response_f_dp=None,
            ):
    """Derive f' from f" """
    
    denergy = energy[1:]-energy[:-1]
    dE = np.min(denergy)
    if np.any(denergy-denergy[0]): # nonuniform spacing
        uniform_mesh = np.arange(energy[0],energy[-1],dE)
        interp = INTERP_FUNC(energy, f_dp)
        f_in=interp(uniform_mesh) # interpolated f_dp
        energy_interp=uniform_mesh
    else:
        f_in = f_dp 
        energy_interp=energy

    if known_response_energy is not None:
        known_response_f_p_interp = INTERP_FUNC(known_response_energy, known_response_f_p)(energy_interp)
        known_response_f_dp_interp = INTERP_FUNC(known_response_energy, known_response_f_dp)(energy_interp)
    else:
        known_response_f_p_interp = np.zeros_like(f_in)
        known_response_f_dp_interp = np.zeros_like(f_in)
    
    f_in = f_in - known_response_f_dp_interp
        
    f_in = apply_window(f_in, padn, trim=trim, window_type=window_type)
    
    f_p_pred_padded = hilbert_transform_func(f_in)
    
    if padn != 0:
        f_p_pred = f_p_pred_padded[padn:-padn]
        start_energy = energy_interp[0]-padn*dE
        end_energy = energy_interp[-1]+(padn+1)*dE
        energy_interp_padded = np.arange(start_energy,end_energy,dE)
    else:
        f_p_pred = f_p_pred_padded
        energy_interp_padded = energy_interp
    
    if include_Z_term:
        Z_star = Z - (Z/82.5)**2.37
        f_p_pred = Z_star + f_p_pred
    
    f_p_pred = f_p_pred + known_response_f_p_interp
    return(energy_interp, f_p_pred, energy_interp_padded, f_p_pred_padded,f_in)


def get_f_dp(energy, f_p,
             padn=5000,
             trim=0,
             Z=26, # atomic number
             include_Z_term=False,
             hilbert_transform_func=hilbert_transform,
             window_type='cosine',
             known_response_energy=None,
             known_response_f_p=None,
             known_response_f_dp=None,
             ):
    
    """Derive f" from f' """
    
    if known_response_f_p is not None:
        known_response_f_p = -known_response_f_p
        
    energy_interp,f_dp_pred, energy_interp_padded, f_dp_pred_padded,_ = get_f_p(energy, -f_p, padn=padn,
                                                                                trim=trim,
                                                                                Z=Z, # atomic number
                                                                                include_Z_term=include_Z_term,
                                                                                hilbert_transform_func=hilbert_transform_func,
                                                                                window_type=window_type,
                                                                                known_response_energy=known_response_energy,
                                                                                known_response_f_p=known_response_f_dp,
                                                                                known_response_f_dp=known_response_f_p,
                                                                                )
    return(energy_interp,f_dp_pred, energy_interp_padded, f_dp_pred_padded)


def apply_window(f_in, padn,
                 trim=0,
                 window_type='cosine'):
    """
    Possible window_type:
    boxcar, triang, blackman, hamming, hann, bartlett, flattop, parzen, bohman,
    blackman, harris, nuttall, barthann, cosine, exponential, tukey, taylor,
    lanczos
    """

    window = get_window(window_type,(padn+trim)*2)
    window = window[0:padn+trim]
    window = np.expand_dims(window, 0)

    
    window_start = f_in[0]*np.ones(padn+trim)
    window_start[padn:]=f_in[0:trim]
    
    window_end = f_in[-1]*np.ones(padn+trim)
    window_end[0:trim]=f_in[len(f_in)-trim:]
    
    windowed_func = np.hstack((window_start * np.squeeze(window),
                               f_in[trim:len(f_in)-trim],
                               window_end * np.squeeze(np.fliplr(window))))

    return(windowed_func) 
    
    
    
def penalty(energy, f_p, f_dp, trim=0, padn=5000, Z=26, include_Z_term=False,
            hilbert_transform_func=hilbert_transform,window_type='cosine',
            known_response_energy=None,
            known_response_f_p=None,
            known_response_f_dp=None,
            ):
    """How close f' and f" are to obeying the Kramers Kronig relation?"""
    
    """Going from f_dp to f_p"""
    energy_interp,f_p_pred,_,f_p_pred_pad,_ = get_f_p(energy, f_dp, trim=trim, padn=padn,
                                                      Z = Z, # atomic number
                                                      include_Z_term=include_Z_term,
                                                      hilbert_transform_func=hilbert_transform_func,
                                                      window_type=window_type,
                                                      known_response_energy=known_response_energy,
                                                      known_response_f_p=known_response_f_p,
                                                      known_response_f_dp=known_response_f_dp,
                                                      )
    
    # add back DC term
    F_p_pred = np.fft.fft(f_p_pred_pad)
    F_p_pred[0] = np.fft.fft(f_p)[0]
    f_p_pred_pad = np.fft.ifft(F_p_pred).real
    
    f_p_pred_pad = f_p_pred_pad[padn:len(f_p_pred_pad)-padn]
    f_p_pred = f_p_pred_pad[trim:len(f_p_pred_pad)-trim]
    
    # interpolation in case of nonuniform energy
    f_p_interp = INTERP_FUNC(energy, f_p)(energy_interp)[trim:len(energy_interp)-trim]
    

    
    """Going from f_p to f_dp"""
    energy_interp,f_dp_pred,_,f_dp_pred_pad = get_f_dp(energy, f_p, trim=trim, padn=padn,
                                                       Z = Z, # atomic number
                                                       include_Z_term=include_Z_term,
                                                       hilbert_transform_func=hilbert_transform_func,
                                                       window_type=window_type,
                                                       known_response_energy=known_response_energy,
                                                       known_response_f_p=known_response_f_p,
                                                       known_response_f_dp=known_response_f_dp,)
    
    # add back DC term
    F_dp_pred = np.fft.fft(f_dp_pred_pad)
    F_dp_pred[0] = np.fft.fft(f_dp)[0]
    f_dp_pred_pad = np.fft.ifft(F_dp_pred).real
    
    f_dp_pred_pad = f_dp_pred_pad[padn:len(f_dp_pred_pad)-padn]
    f_dp_pred = f_dp_pred_pad[trim:len(f_dp_pred_pad)-trim]
    
    # interpolation in case of nonuniform energy
    f_dp_interp = INTERP_FUNC(energy, f_dp)(energy_interp)[trim:len(energy_interp)-trim]

    mse = np.mean((f_p_interp - f_p_pred)**2) + np.mean((f_dp_interp - f_dp_pred)**2)
    return(mse)


def filter_out_dc(f_dp):
    # filter out DC term
    F_dp = np.fft.fft(f_dp)
    F_dp[0] = 0
    f_dp = np.fft.ifft(F_dp).real
    return(f_dp)
