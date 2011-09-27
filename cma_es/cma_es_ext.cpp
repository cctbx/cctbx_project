#include <string>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>

#include "cma/cmaes_interface.h"

/* ============================================================================
   Wrapper for the Covariance Matrix Adaptation Evolution Strategy (CMA-ES)

   The cma directory contains the original source code from

   http://www.lri.fr/~hansen/cmaesintro.html

   The only change is that cmaes.cpp is a copy of cmaes.c (for building the
   shared library in scons).

   This C++ class just keeps track of the cmaes_t struct and exposes the most
   basic functions for minimization to Python.  See example1.c, example2.c,
   and tst_cma_es.py for how to use the minimizer.

   Constructors:
     cma_es - dimensions (int), initial guess (double array), standard
              deviations for guesses (double array)
       This is the most basic constructor providing the bare minimum
       information required for the minimizer.  All other parameters are set
       to their defaults.
     cma_es - file name (string)
       This is the more advanced constructor where all parameters can be
       specified in the file.  See initials.par for the format.

   Member functions:
     sample_population - returns the points to sample
       This returns a 2-d array where the first dimension is the population
       size and the second dimension is the number of dimensions for the target
       function (size of guess)
     update_distribution - updated target function values (double array)
       Updates the minimizer with the new target function values
     converged - returns true is the minimization has converged
     get_result - returns the best guess ever encountered
   ----------------------------------------------------------------------------
*/

namespace cma_es {

  // wrapper class for calling CMA-ES minimizer based on
  // example1.c and example2.c
  class cma_es {

  public:
    cma_es(const int&, scitbx::af::ref<double>, scitbx::af::ref<double>);
    cma_es(std::string);
    ~cma_es();
    scitbx::af::versa<double,scitbx::af::c_grid<2> > sample_population();
    void update_distribution(const scitbx::af::const_ref<double>&);
    bool converged();
    scitbx::af::shared<double> get_result();

  private:
    int N, pop_size;
    cmaes_t evo;
    double *arFunvals, *const*pop, *xfinal;
  };

  // default parameters, only dimension and initial guess required
  cma_es::cma_es(const int& dimension, scitbx::af::ref<double> x,
                 scitbx::af::ref<double> std_dev) {
    N = dimension;
    arFunvals = cmaes_init(&evo,N,x.begin(),std_dev.begin(),0,0,"non");
    pop_size = cmaes_Get(&evo,"popsize");
  }

  // all parameters from file (see initials.par)
  cma_es::cma_es(std::string parameters) {
    arFunvals = cmaes_init(&evo,0,NULL,NULL,0,0,parameters.c_str());
    N = cmaes_Get(&evo,"dim");
    pop_size = cmaes_Get(&evo,"popsize");
  }

  cma_es::~cma_es() {
    cmaes_exit(&evo);
  }

  // returns new population for function evaluation
  scitbx::af::versa<double,scitbx::af::c_grid<2> >
  cma_es::cma_es::sample_population() {
    pop = cmaes_SamplePopulation(&evo);
    scitbx::af::versa<double,scitbx::af::c_grid<2> >
      p(scitbx::af::c_grid<2>(pop_size,N));
    for (int i=0; i<pop_size; i++) {
      for (int j=0; j<N; j++) {
        p(i,j) = pop[i][j];
      }
    }
    return p;
  }

  // updates minizmizer with new function values
  void cma_es::cma_es::update_distribution
  (const scitbx::af::const_ref<double>& new_function_values) {
    cmaes_UpdateDistribution(&evo,new_function_values.begin());
  }

  bool cma_es::cma_es::converged() {
    return cmaes_TestForTermination(&evo);
  }

  scitbx::af::shared<double> cma_es::cma_es::get_result() {
    xfinal = cmaes_GetNew(&evo,"xbestever");
    scitbx::af::shared<double> result(&xfinal[0],&xfinal[0] + N);
    return result;
  }

  // wrapper for Python
  namespace boost_python {
    struct cma_es_wrapper
    {
      static void
      wrap()
      {
        using namespace boost::python;
        class_<cma_es>("cma_es",
                       init<const int&,scitbx::af::ref<double>,
                            scitbx::af::ref<double> >() )
          .def(init<std::string>() )
          .def("sample_population",
               &cma_es::cma_es::sample_population)
          .def("update_distribution",
               &cma_es::cma_es::update_distribution)
          .def("converged",
               &cma_es::cma_es::converged)
          .def("get_result",
               &cma_es::cma_es::get_result)
          ;
      }
    };
  }
}

BOOST_PYTHON_MODULE(cma_es_ext)
{
  cma_es::boost_python::cma_es_wrapper::wrap();
}
