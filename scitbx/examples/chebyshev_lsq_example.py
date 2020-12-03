from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_lsq_fit
from six.moves import cStringIO as StringIO
from six.moves import range
from six.moves import zip



def example():
  x_obs = (flex.double(range(100))+1.0)/101.0
  y_ideal = flex.sin(x_obs*6.0*3.1415) + flex.exp(x_obs)
  y_obs = y_ideal + (flex.random_double(size=x_obs.size())-0.5)*0.5
  w_obs = flex.double(x_obs.size(),1)
  print("Trying to determine the best number of terms ")
  print(" via cross validation techniques")
  print()
  n_terms = chebyshev_lsq_fit.cross_validate_to_determine_number_of_terms(
    x_obs,y_obs,w_obs,
    min_terms=5 ,max_terms=20,
    n_goes=20,n_free=20)
  print("Fitting with", n_terms, "terms")
  print()
  fit = chebyshev_lsq_fit.chebyshev_lsq_fit(n_terms,x_obs,y_obs)
  print("Least Squares residual: %7.6f" %(fit.f))
  print("  R2-value            : %7.6f" %(fit.f/flex.sum(y_obs*y_obs)))
  print()
  fit_funct = chebyshev_polynome(
    n_terms, fit.low_limit, fit.high_limit, fit.coefs)

  y_fitted = fit_funct.f(x_obs)
  abs_deviation = flex.max(
    flex.abs( (y_ideal- y_fitted) ) )
  print("Maximum deviation between fitted and error free data:")
  print("    %4.3f" %(abs_deviation))
  abs_deviation = flex.mean(
    flex.abs( (y_ideal- y_fitted) ) )
  print("Mean deviation between fitted and error free data:")
  print("    %4.3f" %(abs_deviation))
  print()
  abs_deviation = flex.max(
    flex.abs( (y_obs- y_fitted) ) )
  print("Maximum deviation between fitted and observed data:")
  print("    %4.3f" %(abs_deviation))
  abs_deviation = flex.mean(
    flex.abs( (y_obs- y_fitted) ) )
  print("Mean deviation between fitted and observed data:")
  print("    %4.3f" %(abs_deviation))
  print()
  print("Showing 10 points")
  print("   x    y_obs y_ideal y_fit")
  for ii in range(10):
    print("%6.3f %6.3f %6.3f %6.3f" \
          %(x_obs[ii*9], y_obs[ii*9], y_ideal[ii*9], y_fitted[ii*9]))

  try:
    from iotbx import data_plots
  except ImportError:
    pass
  else:
    print("Preparing output for loggraph in a file called")
    print("   chebyshev.loggraph")
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
    with open('chebyshev.loggraph', 'w') as output_logfile:
      f = StringIO()
      data_plots.plot_data_loggraph(chebyshev_plot, f)
      output_logfile.write(f.getvalue())


def another_example(np=41,nt=5):
  x = flex.double( range(np) )/(np-1)
  y = 0.99*flex.exp(-x*x*0.5)
  y = -flex.log(1.0/y-1)
  w = y*y/1.0
  d = (flex.random_double(np)-0.5)*w
  y_obs = y+d

  y = 1.0/( 1.0 + flex.exp(-y) )

  fit_w = chebyshev_lsq_fit.chebyshev_lsq_fit(nt,
                                              x,
                                              y_obs,
                                              w )
  fit_w_f = chebyshev_polynome(
    nt, fit_w.low_limit, fit_w.high_limit, fit_w.coefs)


  fit_nw = chebyshev_lsq_fit.chebyshev_lsq_fit(nt,
                                              x,
                                              y_obs)
  fit_nw_f = chebyshev_polynome(
    nt, fit_nw.low_limit, fit_nw.high_limit, fit_nw.coefs)
  print()
  print("Coefficients from weighted lsq")
  print(list( fit_w.coefs ))
  print("Coefficients from non-weighted lsq")
  print(list( fit_nw.coefs ))
  assert flex.max( flex.abs(fit_nw.coefs-fit_w.coefs) ) > 0

def runge_phenomenon(self,n=41,nt=35,print_it=False):
  x_e = 2.0*(flex.double( range(n) )/float(n-1)-0.5)
  y_e = 1/(1+x_e*x_e*25)
  fit_e = chebyshev_lsq_fit.chebyshev_lsq_fit(nt,
                                              x_e,
                                              y_e,
                                              )
  fit_e = chebyshev_polynome(
    nt, fit_e.low_limit, fit_e.high_limit, fit_e.coefs)


  x_c = chebyshev_lsq_fit.chebyshev_nodes(n, -1, 1, True)
  y_c = 1/(1+x_c*x_c*25)
  fit_c = chebyshev_lsq_fit.chebyshev_lsq_fit(nt,
                                              x_c,
                                              y_c,
                                              )
  fit_c = chebyshev_polynome(
    nt, fit_c.low_limit, fit_c.high_limit, fit_c.coefs)


  x_plot = 2.0*(flex.double( range(3*n) )/float(3*n-1)-0.5)
  y_plot_e = fit_e.f( x_plot )
  y_plot_c = fit_c.f( x_plot )
  y_id =  1/(1+x_plot*x_plot*25)
  if print_it:
    for x,y,yy,yyy in zip(x_plot,y_id,y_plot_e,y_plot_c):
      print(x,y,yy,yyy)







if (__name__ == "__main__"):
  example()
  another_example()
  runge_phenomenon(10)
