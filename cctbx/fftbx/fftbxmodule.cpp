// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Nov 03: Created (R.W. Grosse-Kunstleve)
 */

#include <boost/python/class_builder.hpp>
#include <cctbx/fftbx/complex_to_complex_3d.h>
#include <cctbx/fftbx/real_to_complex_3d.h>
#include <cctbx/std_vector_bpl.h>
#include <cctbx/basic/boost_array_bpl.h>

using namespace cctbx;

namespace {

  typedef std::vector<double> vector_type;

  void rc_forward(fftbx::real_to_complex<vector_type>& t,
                  std::vector<double>& v) {
    t.forward(v);
  }
  void rc_backward(fftbx::real_to_complex<vector_type>& t,
                   std::vector<double>& v) {
    t.backward(v);
  }

  boost::array<std::size_t, 3>
  cc_3d_N(fftbx::complex_to_complex_3d<vector_type>& t) {
    return t.N();
  }
  void cc_3d_forward(fftbx::complex_to_complex_3d<vector_type>& t,
                     std::vector<double>& v) {
    t.forward(v);
  }
  void cc_3d_backward(fftbx::complex_to_complex_3d<vector_type>& t,
                      std::vector<double>& v) {
    t.backward(v);
  }

  void rc_3d_forward(fftbx::real_to_complex_3d<vector_type>& t,
                     std::vector<double>& v) {
    t.forward(v);
  }
  void rc_3d_backward(fftbx::real_to_complex_3d<vector_type>& t,
                      std::vector<double>& v) {
    t.backward(v);
  }

#   include <cctbx/basic/from_bpl_import.h>

  void init_module(python::module_builder& this_module)
  {
    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    (void) python::wrap_std_vector(this_module,
      "vector_of_std_size_t", std::size_t());
    (void) python::wrap_std_vector(this_module,
      "vector_of_double", double());

    class_builder<fftbx::factorization>
    py_factorization(this_module, "factorization");
    class_builder<fftbx::complex_to_complex<vector_type> >
    py_complex_to_complex(this_module, "complex_to_complex");
    class_builder<fftbx::real_to_complex<vector_type> >
    py_real_to_complex(this_module, "real_to_complex");
    class_builder<fftbx::complex_to_complex_3d<vector_type> >
    py_complex_to_complex_3d(this_module, "complex_to_complex_3d");
    class_builder<fftbx::real_to_complex_3d<vector_type> >
    py_real_to_complex_3d(this_module, "real_to_complex_3d");

    py_complex_to_complex.declare_base(
      py_factorization, boost::python::without_downcast);
    py_real_to_complex.declare_base(
      py_factorization, boost::python::without_downcast);

    py_factorization.def(constructor<>());
    py_factorization.def(constructor<std::size_t, bool>());
    py_factorization.def(&fftbx::factorization::N, "N");
    py_factorization.def(&fftbx::factorization::Factors, "Factors");

    py_complex_to_complex.def(constructor<>());
    py_complex_to_complex.def(constructor<std::size_t>());
    py_complex_to_complex.def(
      &fftbx::complex_to_complex<vector_type>::WA, "WA");
    py_complex_to_complex.def(
      &fftbx::complex_to_complex<vector_type>::forward, "forward");
    py_complex_to_complex.def(
      &fftbx::complex_to_complex<vector_type>::backward, "backward");

    py_real_to_complex.def(constructor<>());
    py_real_to_complex.def(constructor<std::size_t>());
    py_real_to_complex.def(&fftbx::real_to_complex<vector_type>::WA, "WA");
    py_real_to_complex.def(rc_forward, "forward");
    py_real_to_complex.def(rc_backward, "backward");

    py_complex_to_complex_3d.def(constructor<>());
    py_complex_to_complex_3d.def(
      constructor<std::size_t, std::size_t, std::size_t>());
    py_complex_to_complex_3d.def(
      constructor<const boost::array<std::size_t, 3>&>());
    py_complex_to_complex_3d.def(cc_3d_N, "N");
    py_complex_to_complex_3d.def(cc_3d_forward, "forward");
    py_complex_to_complex_3d.def(cc_3d_backward, "backward");

    py_real_to_complex_3d.def(constructor<>());
    py_real_to_complex_3d.def(
      constructor<std::size_t, std::size_t, std::size_t>());
    py_real_to_complex_3d.def(
      constructor<const boost::array<std::size_t, 3>&>());
    py_real_to_complex_3d.def(rc_3d_forward, "forward");
    py_real_to_complex_3d.def(rc_3d_backward, "backward");
  }

}

BOOST_PYTHON_MODULE_INIT(fftbx)
{
  try {
    boost::python::module_builder this_module("fftbx");
    init_module(this_module);
  }
  catch(...) {
    boost::python::handle_exception();
  }
}
