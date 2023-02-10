"""
Example optimization using functions in kramers_kronig_optimze.
This script is meant to test importing and running the functions.
"""

import kramers_kronig.kramers_kronig_optimize as kramers_kronig_optimize

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
                                                           padn=100,
                                                           trim=30,
                                                           spacing=10,
                                                           noise_loc=[0,0],
                                                           noise_scale=[1e-3,1e-3],
                                                           learning_rate=1e-1,
                                                           num_iter=10000,
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