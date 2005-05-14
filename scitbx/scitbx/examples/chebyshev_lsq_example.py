from cctbx.array_family import flex
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_lsq_fit
from iotbx import data_plots
import math



def example():
  x_obs = (flex.double(range(100))+1.0)/101.0
  y_ideal = flex.sin(x_obs*6.0*3.1415) + flex.exp(x_obs)
  y_obs = y_ideal + (flex.random_double(size=x_obs.size())-0.5)*0.5
  w_obs = flex.double(x_obs.size(),1)
  print "Trying to determine the best number of terms "
  print " via cross validation techniques"
  print
  n_terms = chebyshev_lsq_fit.cross_validate_to_determine_number_of_terms(
    x_obs,y_obs,w_obs,
    min_terms=5 ,max_terms=20,
    n_goes=20,n_free=20)
  print "Fitting with", n_terms, "terms"
  print
  fit = chebyshev_lsq_fit.chebyshev_lsq_fit(n_terms,x_obs,y_obs)
  print "Least Squares residual: %7.6f" %(fit.f)
  print "  R2-value            : %7.6f" %(fit.f/flex.sum(y_obs*y_obs))
  print
  fit_funct = chebyshev_polynome(
    n_terms, fit.low_limit, fit.high_limit, fit.coefs)

  y_fitted = fit_funct.f(x_obs)
  abs_deviation = flex.max(
    flex.abs( (y_ideal- y_fitted) ) )
  print "Maximum deviation between fitted and error free data:"
  print "    %4.3f" %(abs_deviation)
  abs_deviation = flex.mean(
    flex.abs( (y_ideal- y_fitted) ) )
  print "Mean deviation between fitted and error free data:"
  print "    %4.3f" %(abs_deviation)
  print
  abs_deviation = flex.max(
    flex.abs( (y_obs- y_fitted) ) )
  print "Maximum deviation between fitted and observed data:"
  print "    %4.3f" %(abs_deviation)
  abs_deviation = flex.mean(
    flex.abs( (y_obs- y_fitted) ) )
  print "Mean deviation between fitted and observed data:"
  print "    %4.3f" %(abs_deviation)
  print
  print "Showing 10 points"
  print "   x    y_obs y_ideal y_fit"
  for ii in range(10):
    print "%6.3f %6.3f %6.3f %6.3f" \
          %(x_obs[ii*9], y_obs[ii*9], y_ideal[ii*9], y_fitted[ii*9])

  print "Preparing output for loggraph in a file called"
  print "   chebyshev.loggraph"
  chebyshev_plot = data_plots.plot_data(plot_title='Chebyshev fitting',
                                        x_label = 'x values',
                                        y_label = 'y values',
                                        x_data = x_obs,
                                        y_data = y_obs,
                                        y_legend = 'Observed y values',
                                        comments = 'Chebyshev fit')
  chebyshev_plot.add_data(y_data=y_ideal,
                          y_legend='Error free y values')
  chebyshev_plot.add_data(y_data=y_fitted,
                          y_legend='Fitted chebyshev approximation')
  output_logfile=open('chebyshev.loggraph','w')
  data_plots.plot_data_loggraph(chebyshev_plot,output_logfile)





if (__name__ == "__main__"):
  example()
