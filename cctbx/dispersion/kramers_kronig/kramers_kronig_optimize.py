"""Functions to create test optimizations with a kramers-kronig penalty."""

from __future__ import division

import numpy as np
import torch
import matplotlib.pyplot as plt

from . import kramers_kronig
from . import kramers_kronig_helper

def create_f(width=10,
             dE=.1,
             trim=0,
             slope = 1,
             padn=5000,
             uniform_energy=True
             ):
    """Create a simulated absorption edge"""

    if uniform_energy:
        energy = np.arange(-width,width,dE)
    else:
        energy_0 = np.arange(-width,0,dE)
        energy_1 = np.arange(0,width, 2*dE)
        energy = np.concatenate((energy_0, energy_1))

    # ramp/unit step function for f" to emulate a simple K-edge
    ramp = energy*slope
    ramp_start = np.argmin(np.abs(ramp))
    ramp_end = np.argmin(np.abs(ramp-1))
    ramp_size = ramp_end-ramp_start
    f_dp = np.heaviside(energy, 0.5)
    mid_ind = np.argmin(np.abs(energy))
    f_dp[mid_ind:mid_ind+ramp_size] = ramp[ramp_start:ramp_end]

    energy = torch.Tensor(energy)
    f_dp = torch.Tensor(f_dp)

    # get f' from the Hilbert transform
    energy_padded,_,f_p_padded,f_dp_padded = \
    kramers_kronig.get_f_p(energy, f_dp, padn=padn,
                           trim=trim,
                           )
    return(energy_padded,f_p_padded,f_dp_padded)

def sample(f_p,
           f_dp,
           loc=[0,0],
           scale=[1e-3,1e-3],
           ):
    """Add Gaussian noise to the f' and f" curves and sample"""
    f_p_dist = torch.distributions.normal.Normal(f_p + loc[0], scale[0])
    f_dp_dist = torch.distributions.normal.Normal(f_dp + loc[1], scale[1])
    return(f_p_dist.sample(),f_dp_dist.sample())

def subsample(energy,
              f_p,
              f_dp,
              spacing=2,
              ):
    """Subsample the f' and f" curves at an input spacing"""
    inds = np.arange(0,len(energy),spacing)
    energy = energy[inds]
    f_p = f_p[inds]
    f_dp = f_dp[inds]
    return(energy,f_p,f_dp,inds)

def get_loss(energy,
             f_p_opt,
             f_dp_opt,
             f_p_noisy_ss,
             f_dp_noisy_ss,
             inds,
             known_response_energy=None,
             known_response_f_p=None,
             known_response_f_dp=None,
             ):
    """Optimization loss is the MSE of the match to the noisy subsampled data summed with the penalty for
    violating the Kramers-Kronig relations"""
    data_loss = torch.mean((f_p_opt[inds]-f_p_noisy_ss)**2 + (f_dp_opt[inds]-f_dp_noisy_ss)**2)
    kk_loss = kramers_kronig.get_penalty(energy, f_p_opt, f_dp_opt, padn=0, trim=0,
                                         known_response_energy=known_response_energy,
                                         known_response_f_p=known_response_f_p,
                                         known_response_f_dp=known_response_f_dp,
                                         )

    return(data_loss + kk_loss)

def run_example_opt(width=5,
                    padn=100,
                    trim=30,
                    spacing=5,
                    noise_loc=[0,0],
                    noise_scale=[1e-1,1e-1],
                    learning_rate=1e-1,
                    num_iter=10000,
                    uniform_energy=True,
                    known_response_energy=None,
                    known_response_f_p=None,
                    known_response_f_dp=None,
                    ):
    """
    Run an example optimization with the following steps:
    1. Simulate f" based on a very simple model of the K-edge.
    2. Use the dispersion relations to calculate f′.
    3. Sample both of these curves with Gaussian noise to simulate experimental measurement of the two curves.
    4. Use restraint to optimize the parameters. Use automatic differentiation for first-derivatives.
    """

    energy,f_p,f_dp = create_f(width=width,
                               padn=padn,
                               trim=trim,
                               uniform_energy=uniform_energy)

    f_p_noisy,f_dp_noisy = sample(f_p,f_dp,
                                  loc=noise_loc,
                                  scale=noise_scale,
                                  )
    energy_ss,f_p_noisy_ss,f_dp_noisy_ss,inds = subsample(energy,f_p_noisy,f_dp_noisy,
                                                          spacing=spacing)

    # From energy_ss,f_p_noisy_ss, and f_dp_noisy_ss, determine f_p and f_dp, energy is given
    f_p_pred_0 = kramers_kronig_helper.INTERP_FUNC(energy_ss,f_p_noisy_ss)(energy)
    f_dp_pred_0 = kramers_kronig_helper.INTERP_FUNC(energy_ss,f_dp_noisy_ss)(energy)

    f_p_opt = torch.tensor(f_p_pred_0,requires_grad=True)
    f_dp_opt = torch.tensor(f_dp_pred_0, requires_grad=True)

    optimizer = torch.optim.SGD([f_p_opt,f_dp_opt],lr=learning_rate)

    loss_vec = []
    actual_loss_vec = []
    for i in range(num_iter):
        loss = get_loss(energy,
                        f_p_opt,
                        f_dp_opt,
                        f_p_noisy_ss,
                        f_dp_noisy_ss,
                        inds,
                        known_response_energy=known_response_energy,
                        known_response_f_p=known_response_f_p,
                        known_response_f_dp=known_response_f_dp,
                        )
        actual_loss = torch.mean((f_p_opt-f_p)**2 + (f_dp_opt-f_dp)**2)
        loss_vec.append(loss)
        actual_loss_vec.append(actual_loss)

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    return(energy,
           f_p,
           f_dp,
           energy_ss,
           f_p_noisy_ss,
           f_dp_noisy_ss,
           f_p_pred_0,
           f_dp_pred_0,
           f_p_opt,
           f_dp_opt,
           loss_vec,
           actual_loss_vec,
           )

def visualize(energy,
              f_p,
              f_dp,
              energy_ss,
              f_p_noisy_ss,
              f_dp_noisy_ss,
              f_p_pred_0,
              f_dp_pred_0,
              f_p_opt,
              f_dp_opt,
              loss_vec,
              actual_loss_vec,
              ):

    """Compare the optimized model to the initial ground truth
    and show result in matplotlib."""

    plt.figure()
    plt.title("Subsampled curves with noise")
    plt.plot(energy,f_dp,energy_ss,f_dp_noisy_ss,'.')
    plt.plot(energy,f_p,energy_ss,f_p_noisy_ss,'.')
    plt.xlim([energy[0],energy[-1]])
    plt.show()

    plt.figure()
    plt.title('Initial Guess for actual curves')
    plt.plot(energy,f_dp,energy,f_dp_pred_0,'.')
    plt.plot(energy,f_p,energy,f_p_pred_0,'.')
    plt.xlim([energy[0],energy[-1]])
    plt.show()

    plt.figure()
    plt.title('Final Guess for actual curves')
    plt.plot(energy,f_dp,energy,f_dp_opt.detach().numpy(),'.')
    plt.plot(energy,f_p,energy,f_p_opt.detach().numpy(),'.')
    plt.xlim([energy[0],energy[-1]])
    plt.show()

    plt.figure()
    plt.title('Objective Loss')
    plt.plot([loss.detach().numpy() for loss in loss_vec])
    plt.show()

    plt.figure()
    plt.title('Ground Truth Loss')
    plt.plot([actual_loss.detach().numpy() for actual_loss in actual_loss_vec])
    plt.show()
