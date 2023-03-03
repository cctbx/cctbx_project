"""
Example optimization using functions in kramers_kronig_optimze.py.
In this script, the following is completed:
1. f" is simulated based on a very simple model of the K-edge (a ramp function that saturates).
2. The Hilbert tranform is used to calculate f'.
3. f" and f' are contaminated with Gaussian noise.
4. f" and f' are subsampled.
5. The optimization attempts to determine f" and f' by both trying to match the subsampled values and
    minimize a penalty that restrains f" and f' to be the Hilbert transform and negative Hilbert transform
    of each other respectively.
6. The final results are plotted in matplotlib.
"""
from __future__ import division

import kramers_kronig.kramers_kronig_optimize as kramers_kronig_optimize

known_response_energy = None
known_response_f_p = None
known_response_f_dp = None

energy,\
f_p,\
f_dp,\
energy_ss,\
f_p_noisy_ss,\
f_dp_noisy_ss,\
f_p_pred_0,\
f_dp_pred_0,\
f_p_opt,\
f_dp_opt,\
loss_vec,\
actual_loss_vec  = kramers_kronig_optimize.run_example_opt(width=5,
                                                           padn=10,
                                                           trim=30,
                                                           spacing=10,
                                                           noise_loc=[0,0],
                                                           noise_scale=[1e-3,1e-3],
                                                           learning_rate=1e-1,
                                                           num_iter=1000,
                                                           uniform_energy=True,
                                                           known_response_energy=known_response_energy,
                                                           known_response_f_p=known_response_f_p,
                                                           known_response_f_dp=known_response_f_dp,
                                                           )

kramers_kronig_optimize.visualize(energy,
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
