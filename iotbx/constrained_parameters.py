""" Tags to identify commonplace constraints (as those featured by ShelXL)
"""

(independent_parameter,
 constant_parameter,
 constant_times_independent_scalar_parameter_minus_1, # c*(x-1)
 constant_times_independent_scalar_parameter        , # c*x
                                                      # where c: constant
                                                      # and   x: parameter
 constant_times_u_eq) = xrange(5)
