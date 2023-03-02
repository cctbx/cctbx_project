"""Functions to compute penalty for f' and f" violating Kramers Kronig relations, written in PyTorch"""

from __future__ import division

import numpy as np
import torch

from scipy.signal.windows import get_window
from . import kramers_kronig_helper

def get_hilbert_transform(x, axis=-1):
    """
    Perform the Hilbert transform on a real-valued input.
    This is a PyTorch implementation of the Hilbert transform in scipy.signal
    Reference: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.hilbert.html
    Like the scipy version, this version takes as input a real-valued function.
    Assuming the function is analytic, the full complex function is computed.
    Unlike the scipy version, this version only returns the imaginary part of the analytic function.
    scipy returns a complex vector where the real part is the input and the imaginary part is the same as
    the output of this function.
    Saying that the real and imaginary part of a function are related by Kramers-Kronig is the same as saying
    that the imaginary part of a function is the Hilbert transform of the real part.
    """

    N = x.shape[axis]
    if N <= 0:
        raise ValueError("N must be positive.")

    Xf = torch.fft.fft(x, N, axis=axis)
    h = torch.zeros(N, dtype=Xf.dtype, requires_grad=False)
    if N % 2 == 0:
        h[0] = h[N // 2] = 1
        h[1:N // 2] = 2
    else:
        h[0] = 1
        h[1:(N + 1) // 2] = 2

    if x.ndim > 1:
        ind = [np.newaxis] * x.ndim
        ind[axis] = slice(None)
        h = h[tuple(ind)]
    x = torch.fft.ifft(Xf * h, axis=axis)
    return(x.imag)

def get_f_p(energy,
            f_dp,
            padn=5000,
            trim=0,
            window_type='cosine',
            known_response_energy=None,
            known_response_f_p=None,
            known_response_f_dp=None,
            ):
    """
    This function calculates f' (f_p) from f" (f_dp).
    The input is f" as a function of energy. If the function is not on a uniform grid, f"
    is interpolated, using the smallest spacing of the given nonuniform grid.
    The Hilbert transform is linear. Thus, if a pair of f' and f" has already been calculated from the
    Hilbert transform, the known f" can be subtracted off of the input f". After f' is calculated from the input f"
    through the Hilbert transform, the known f' can be added back to the calculated f' for the final value.
    The Hilbert transform requires knowledge of the function at all energies. As we generally only have a truncated
    region, this subtraction and addition procedure allows the user to make the tails of the input function go to zero
    at the edges, getting a better approximation of the Hilbert transform.
    The subtracted input function is padded by zeros (of number padn) and windowed by a window of window_type.
    A number of start and endpoints of the input function (of number trim) are included in this window.
    """

    denergy = energy[1:]-energy[:-1]
    if torch.any(torch.abs(denergy-denergy[0])>1e-3): # Energy spacing is not constant.
        energy, f_dp = kramers_kronig_helper.interpolate(energy, f_dp, mode="torch")

    if known_response_energy is not None:
        known_response_f_p_interp = kramers_kronig_helper.INTERP_FUNC(known_response_energy, known_response_f_p)(energy)
        known_response_f_dp_interp = kramers_kronig_helper.INTERP_FUNC(known_response_energy, known_response_f_dp)(energy)
    else:
        known_response_f_p_interp = torch.zeros_like(f_dp)
        known_response_f_dp_interp = torch.zeros_like(f_dp)

    f_in = f_dp - torch.Tensor(known_response_f_dp_interp)
    f_in = apply_window(f_in, padn, trim=trim, window_type=window_type)

    f_p_pred_padded = get_hilbert_transform(f_in)
    f_p_pred_padded[padn:len(f_p_pred_padded)-padn] += torch.Tensor(known_response_f_p_interp)

    if padn != 0:
        f_p_pred = f_p_pred_padded[padn:-padn]
        dE = energy[1] - energy[0]
        start_energy = energy[0]-padn*dE
        if len(energy)%2: # odd
            end_energy = energy[-1]+(padn+1)*dE
        else:
            end_energy = energy[-1]+(padn+0)*dE
        energy_padded = torch.arange(start_energy,end_energy,dE)
    else:
        f_p_pred = f_p_pred_padded
        energy_padded = energy

    return(energy_padded,f_p_pred,f_p_pred_padded,f_in)

def get_f_dp(energy,
             f_p,
             padn=5000,
             trim=0,
             window_type='cosine',
             known_response_energy=None,
             known_response_f_p=None,
             known_response_f_dp=None,
             ):

    """Derive f" from f'. This is calculated from the negative of the Hilbert transform."""

    if known_response_f_p is not None:
        known_response_f_p = -known_response_f_p


    energy_padded,f_dp_pred,f_dp_pred_padded,f_in = get_f_p(energy, -f_p, padn=padn,
                                                            trim=trim,
                                                            window_type=window_type,
                                                            known_response_energy=known_response_energy,
                                                            known_response_f_p=known_response_f_dp,
                                                            known_response_f_dp=known_response_f_p,
                                                            )
    return(energy_padded,f_dp_pred,f_dp_pred_padded,f_in)

def apply_window(f_in, padn,
                 trim=0,
                 window_type='cosine'):
    """
    Apply windowing to f_in
    padn # of points are added to the beginning and end of f_in
    trim # of points are modified at the beginning and end of f_in

    Possible window_type:
    boxcar, triang, blackman, hamming, hann, bartlett, flattop, parzen, bohman,
    blackman, harris, nuttall, barthann, cosine, exponential, tukey, taylor,
    lanczos
    """

    window = get_window(window_type,(padn+trim)*2)
    window = window[0:padn+trim]
    window = torch.tensor(np.expand_dims(window, 0))

    window_start = f_in[0]*torch.ones(padn+trim)
    window_start[padn:]=f_in[0:trim]

    window_end = f_in[-1]*torch.ones(padn+trim)
    window_end[0:trim]=f_in[len(f_in)-trim:]

    windowed_func = torch.hstack((window_start * torch.squeeze(window),
                                  f_in[trim:len(f_in)-trim],
                                  window_end * torch.squeeze(torch.fliplr(window))))

    return(windowed_func)

def get_penalty(energy, f_p, f_dp, trim=0, padn=5000,window_type='cosine',
                known_response_energy=None,
                known_response_f_p=None,
                known_response_f_dp=None,
                ):
    """How close are f' and f" are to obeying the Kramers Kronig relation?"""

    #Going from f_dp to f_p

    energy_padded,f_p_pred,f_p_pred_padded,f_in = get_f_p(energy, f_dp, trim=trim, padn=padn,
                                                          window_type=window_type,
                                                          known_response_energy=known_response_energy,
                                                          known_response_f_p=known_response_f_p,
                                                          known_response_f_dp=known_response_f_dp,
                                                          )

    # The Hilbert transform filters out any DC term. We add back the DC term.
    F_p_pred = torch.fft.fft(f_p_pred_padded)
    F_p_pred[0] = torch.fft.fft(f_p)[0]
    f_p_pred_padded = torch.fft.ifft(F_p_pred).real

    f_p_pred_padded = f_p_pred_padded[padn:len(f_p_pred_padded)-padn]
    f_p_pred = f_p_pred_padded[trim:len(f_p_pred_padded)-trim]

    # Going from f_p to f_dp
    energy_padded,f_dp_pred,f_dp_pred_padded, f_in = get_f_dp(energy, f_p, trim=trim, padn=padn,
                                                              window_type=window_type,
                                                              known_response_energy=known_response_energy,
                                                              known_response_f_p=known_response_f_p,
                                                              known_response_f_dp=known_response_f_dp,
                                                              )

    # Add back DC term
    F_dp_pred = torch.fft.fft(f_dp_pred_padded)
    F_dp_pred[0] = torch.fft.fft(f_dp)[0]
    f_dp_pred_padded = torch.fft.ifft(F_dp_pred).real

    f_dp_pred_padded = f_dp_pred_padded[padn:len(f_dp_pred_padded)-padn]
    f_dp_pred = f_dp_pred_padded[trim:len(f_dp_pred_padded)-trim]

    # trim f_p and f_dp as the start and endpoints (of number trim) are
    # windowed during computation of the Hilbert transform
    f_p = f_p[trim:len(energy)-trim]
    f_dp = f_dp[trim:len(energy)-trim]
    mse = torch.mean((f_p - f_p_pred)**2) + torch.mean((f_dp - f_dp_pred)**2)
    return(mse)
