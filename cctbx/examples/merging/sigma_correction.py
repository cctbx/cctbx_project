from __future__ import absolute_import, division, print_function
from six.moves import range
from scitbx.array_family import flex
class ccp4_model(object):
  """Implement a sigma correction for semi-datasets, inspired by the classic
     treatment of SDFAC, SDB, SDADD as discussed by Phil Evans;
     Evans (2011) Acta Cryst D67, 282-292.
     Evans (2006) Acta Cryst D62, 72-82.
     http://ccp4wiki.org/~ccp4wiki/wiki/index.php?title=Symmetry%2C_Scale%2C_Merge#Analysis_of_Standard_Deviations
  """

  def __init__(self):
    #from libtbx import adopt_init_args
    #adopt_init_args(self, locals())
    return
    """An assumption of the class is that semi-datasets A and B are actually drawn from the
    same population.
    """

  @staticmethod
  def plots(a_data, b_data, a_sigmas, b_sigmas):

    # Diagnostic use of the (I - <I>) / sigma distribution, should have mean=0, std=1
    a_variance = a_sigmas * a_sigmas
    b_variance = b_sigmas * b_sigmas
    mean_num = (a_data/ (a_variance) ) + (b_data/ (b_variance) )
    mean_den = (1./ (a_variance) ) + (1./ (b_variance) )
    mean_values = mean_num / mean_den

    delta_I_a = a_data - mean_values
    normal_a = delta_I_a / (a_sigmas)
    stats_a = flex.mean_and_variance(normal_a)
    print("\nA mean %7.4f std %7.4f"%(stats_a.mean(),stats_a.unweighted_sample_standard_deviation()))
    order_a = flex.sort_permutation(normal_a)

    delta_I_b = b_data - mean_values
    normal_b = delta_I_b / (b_sigmas)
    stats_b = flex.mean_and_variance(normal_b)
    print("B mean %7.4f std %7.4f"%(stats_b.mean(),stats_b.unweighted_sample_standard_deviation()))
    order_b = flex.sort_permutation(normal_b)
    # plots for debugging
    from matplotlib import pyplot as plt
    cumnorm = plt.subplot(321)
    cumnorm.plot(range(len(order_a)),normal_a.select(order_a),"b.")
    cumnorm.plot(range(len(order_b)),normal_b.select(order_b),"r.")
    #plt.show()
    logger = plt.subplot(324)
    logger.loglog(a_data,b_data,"r.")
    delta = plt.subplot(322)
    delta.plot(a_data, delta_I_a, "g.")
    #plt.show()
    #nselection = (flex.abs(normal_a) < 2.).__and__(flex.abs(normal_b) < 2.)
    gam = plt.subplot(323)
    gam.plot(mean_values,normal_a,"b.")
    sigs = plt.subplot(326)
    sigs.plot(a_sigmas,b_sigmas,"g.")
    mean_order = flex.sort_permutation(mean_values)
    scatters = flex.double(50)
    scattersb = flex.double(50)
    for isubsection in range(50):
      subselect = mean_order[isubsection*len(mean_order)//50:(isubsection+1)*len(mean_order)//50]
      vals = normal_a.select(subselect)
      #scatters[isubsection] = flex.mean_and_variance(vals).unweighted_sample_standard_deviation()
      scatters[isubsection] = flex.mean_and_variance(vals).unweighted_sample_variance()

      valsb = normal_b.select(subselect)
      #scatters[isubsection] = flex.mean_and_variance(vals).unweighted_sample_standard_deviation()
      scattersb[isubsection] = flex.mean_and_variance(valsb).unweighted_sample_variance()
    aaronsplot = plt.subplot(325)
    aaronsplot.plot(range(50), 2. * scatters, "b.")
    plt.show()

  @staticmethod
  def apply_sd_error_params(vector, a_data, b_data, a_sigmas, b_sigmas):
    sdfac, sdb, sdadd = vector[0],0,vector[1]
    a_variance = a_sigmas * a_sigmas
    b_variance = b_sigmas * b_sigmas
    mean_num = (a_data/ (a_variance) ) + (b_data/ (b_variance) )
    mean_den = (1./ (a_variance) ) + (1./ (b_variance) )
    mean_values = mean_num / mean_den

    I_mean_dependent_part = sdb * mean_values + flex.pow(sdadd * mean_values, 2)
    a_new_variance = sdfac*sdfac * ( a_variance + I_mean_dependent_part )
    b_new_variance = sdfac*sdfac * ( b_variance + I_mean_dependent_part )
    return a_new_variance, b_new_variance

  @staticmethod
  def optimize(a_data, b_data, a_sigmas, b_sigmas):


    print("""Fit the parameters SDfac, SDB and SDAdd using Nelder-Mead simplex method.""")

    from scitbx.simplex import simplex_opt
    class simplex_minimizer(object):
      """Class for refining sdfac, sdb and sdadd"""
      def __init__(self):
        """
        """
        self.n = 2
        self.x = flex.double([0.5,0.0])
        self.starting_simplex = []
        for i in range(self.n+1):
          self.starting_simplex.append(flex.random_double(self.n))

        self.optimizer = simplex_opt( dimension = self.n,
                                      matrix    = self.starting_simplex,
                                      evaluator = self,
                                      tolerance = 1e-1)
        self.x = self.optimizer.get_solution()

      def target(self, vector):
        """ Compute the functional by first applying the current values for the sd parameters
        to the input data, then computing the complete set of normalized deviations and finally
        using those normalized deviations to compute the functional."""
        sdfac, sdb, sdadd = vector[0],0.0,vector[1]

        a_new_variance, b_new_variance = ccp4_model.apply_sd_error_params(
          vector, a_data, b_data, a_sigmas, b_sigmas)

        mean_num = (a_data/ (a_new_variance) ) + (b_data/ (b_new_variance) )
        mean_den = (1./ (a_new_variance) ) + (1./ (b_new_variance) )
        mean_values = mean_num / mean_den

        delta_I_a = a_data - mean_values
        normal_a = delta_I_a / flex.sqrt(a_new_variance)

        delta_I_b = b_data - mean_values
        normal_b = delta_I_b / flex.sqrt(b_new_variance)

        mean_order = flex.sort_permutation(mean_values)
        scatters = flex.double(50)
        scattersb = flex.double(50)
        for isubsection in range(50):
          subselect = mean_order[isubsection*len(mean_order)//50:(isubsection+1)*len(mean_order)//50]
          vals = normal_a.select(subselect)
          scatters[isubsection] = flex.mean_and_variance(vals).unweighted_sample_variance()

          valsb = normal_b.select(subselect)
          scattersb[isubsection] = flex.mean_and_variance(valsb).unweighted_sample_variance()

        f = flex.sum( flex.pow(1.-scatters, 2) )
        print("f: % 12.1f, sdfac: %8.5f, sdb: %8.5f, sdadd: %8.5f"%(f, sdfac, sdb, sdadd))
        return f
    optimizer = simplex_minimizer()
    print("new parameters:",list(optimizer.x))
    return ccp4_model.apply_sd_error_params(optimizer.x,a_data, b_data, a_sigmas, b_sigmas)
