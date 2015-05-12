/*
 * ext.cc
 *
 *  Copyright (C) 2015 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/make_shared.hpp>
#include <scitbx/glmtbx/robust_glm.h>

namespace scitbx { namespace glmtbx { namespace boost_python {

  using namespace boost::python;

  template <typename Family>
  class_< robust_glm<Family> > wrap_robust_glm(const char *name) {

    typedef robust_glm<Family> glm_type;

    class_<glm_type> wrapper(name, no_init);

    wrapper
      .def("parameters",
          &glm_type::parameters)
      .def("niter",
          &glm_type::niter)
      .def("error",
          &glm_type::error)
      .def("converged",
          &glm_type::converged)
      .def("mu",
          &glm_type::mu, (
            arg("X")))
      ;

    return wrapper;
  }

  object robust_glm_selector(
      const af::const_ref< double, af::c_grid<2> > &X,
      const af::const_ref< double > &Y,
      const af::const_ref< double > &B,
      double c,
      double tolerance,
      std::size_t max_iter,
      const std::string &family) {
    object result;
    if (family == "poisson") {
      result = object(new robust_glm<poisson>(
            X,
            Y,
            B,
            c,
            tolerance,
            max_iter));
    } else {
      SCITBX_ERROR("Unknown distribution");
    }
    return result;
  }

  BOOST_PYTHON_MODULE(scitbx_glmtbx_ext)
  {

    wrap_robust_glm<poisson>("robust_glm_poisson");

    def("robust_glm",
        &robust_glm_selector, (
          arg("X"),
          arg("Y"),
          arg("B"),
          arg("c")=1.345,
          arg("tolerance")=1e-3,
          arg("max_iter")=10,
          arg("family")="poisson"));
  }

}}} // namespace = dials::model::boost_python
