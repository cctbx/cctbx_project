// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Mar 2002: Using af::shared (rwgk)
     2001 Dec 21: Using iterator-based fftbx interface (rwgk)
     2001 Nov 03: Created (R.W. Grosse-Kunstleve)
 */

#include <boost/python/cross_module.hpp>
#include <cctbx/fftbx/gridding.h>
#include <cctbx/fftbx/complex_to_complex_3d.h>
#include <cctbx/fftbx/real_to_complex_3d.h>
#include <cctbx/array_family/tiny_bpl.h>

using namespace cctbx;

namespace {

  int adjust_gridding_2(int min_grid,
                                int max_prime) {
    return fftbx::adjust_gridding(min_grid, max_prime);
  }
  int adjust_gridding_3(int min_grid,
                                int max_prime,
                                int mandatory_factor) {
    return fftbx::adjust_gridding(min_grid, max_prime, mandatory_factor);
  }

  af::int3 adjust_gridding_triple_2(
    const af::int3& min_grid,
    int max_prime) {
    return fftbx::adjust_gridding_array(min_grid, max_prime);
  }

  af::int3 adjust_gridding_triple_3(
    const af::int3& min_grid,
    int max_prime,
    const af::int3& mandatory_factors) {
    return fftbx::adjust_gridding_array(min_grid, max_prime,
                                        mandatory_factors);
  }

  void throw_size_error() {
    PyErr_SetString(PyExc_RuntimeError, "Array is too small.");
    throw boost::python::error_already_set();
  }

  void throw_index_error() {
    PyErr_SetString(PyExc_IndexError, "Index is out of range.");
    throw boost::python::error_already_set();
  }

  typedef af::shared<double> shared_real_array;
  typedef af::ref<double, af::grid<3> > ref_3d_real_array;
  typedef af::shared<std::complex<double> > shared_complex_array;
  typedef af::ref<std::complex<double>, af::grid<3> > ref_3d_complex_array;

  void cc_forward_complex(fftbx::complex_to_complex<double>& fft,
                          shared_complex_array a) {
    if (a.size() < fft.N()) throw_size_error();
    fft.forward(a.begin());
  }
  void cc_backward_complex(fftbx::complex_to_complex<double>& fft,
                           shared_complex_array a) {
    if (a.size() < fft.N()) throw_size_error();
    fft.backward(a.begin());
  }
  void cc_forward_real(fftbx::complex_to_complex<double>& fft,
                       shared_real_array a) {
    if (a.size() < 2 * fft.N()) throw_size_error();
    fft.forward(a.begin());
  }
  void cc_backward_real(fftbx::complex_to_complex<double>& fft,
                        shared_real_array a) {
    if (a.size() < 2 * fft.N()) throw_size_error();
    fft.backward(a.begin());
  }

  void rc_forward_complex(fftbx::real_to_complex<double>& fft,
                          shared_complex_array a) {
    if (a.size() < fft.Ncomplex()) throw_size_error();
    fft.forward(a.begin());
  }
  void rc_backward_complex(fftbx::real_to_complex<double>& fft,
                           shared_complex_array a) {
    if (a.size() < fft.Ncomplex()) throw_size_error();
    fft.backward(a.begin());
  }
  void rc_forward_real(fftbx::real_to_complex<double>& fft,
                       shared_real_array a) {
    if (a.size() < fft.Mreal()) throw_size_error();
    fft.forward(a.begin());
  }
  void rc_backward_real(fftbx::real_to_complex<double>& fft,
                        shared_real_array a) {
    if (a.size() < fft.Mreal()) throw_size_error();
    fft.backward(a.begin());
  }

  void cc_3d_forward_complex(fftbx::complex_to_complex_3d<double>& fft,
                             shared_complex_array a) {
    ref_3d_complex_array map(a.begin(), af::grid<3>(fft.N()));
    if (a.size() < map.size()) throw_size_error();
    fft.forward(map);
  }
  void cc_3d_backward_complex(fftbx::complex_to_complex_3d<double>& fft,
                              shared_complex_array a) {
    ref_3d_complex_array map(a.begin(), af::grid<3>(fft.N()));
    if (a.size() < map.size()) throw_size_error();
    fft.backward(map);
  }
  void cc_3d_forward_real(fftbx::complex_to_complex_3d<double>& fft,
                          shared_real_array a) {
    ref_3d_real_array map(
      a.begin(), af::grid<3>(fftbx::Nreal_from_Ncomplex(fft.N())));
    if (a.size() < map.size()) throw_size_error();
    fft.forward(map);
  }
  void cc_3d_backward_real(fftbx::complex_to_complex_3d<double>& fft,
                           shared_real_array a) {
    ref_3d_real_array map(
      a.begin(), af::grid<3>(fftbx::Nreal_from_Ncomplex(fft.N())));
    if (a.size() < map.size()) throw_size_error();
    fft.backward(map);
  }

  void rc_3d_forward_complex(fftbx::real_to_complex_3d<double>& fft,
                             shared_complex_array a) {
    ref_3d_complex_array map(a.begin(), af::grid<3>(fft.Ncomplex()));
    if (a.size() < map.size()) throw_size_error();
    fft.forward(map);
  }
  void rc_3d_backward_complex(fftbx::real_to_complex_3d<double>& fft,
                              shared_complex_array a) {
    ref_3d_complex_array map(a.begin(), af::grid<3>(fft.Ncomplex()));
    if (a.size() < map.size()) throw_size_error();
    fft.backward(map);
  }
  void rc_3d_forward_real(fftbx::real_to_complex_3d<double>& fft,
                          shared_real_array a) {
    ref_3d_real_array map(a.begin(), af::grid<3>(fft.Mreal()));
    if (a.size() < map.size()) throw_size_error();
    fft.forward(map);
  }
  void rc_3d_backward_real(fftbx::real_to_complex_3d<double>& fft,
                           shared_real_array a) {
    ref_3d_real_array map(a.begin(), af::grid<3>(fft.Mreal()));
    if (a.size() < map.size()) throw_size_error();
    fft.backward(map);
  }

#   include <cctbx/basic/from_bpl_import.h>

  void init_module(python::module_builder& this_module)
  {
    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    python::import_converters<af::shared<int> >
    py_shared_int(
      "cctbx_boost.arraytbx.shared", "int");

    python::import_converters<af::shared<double> >
    py_shared_double(
      "cctbx_boost.arraytbx.shared", "double");

    python::import_converters<af::shared<std::complex<double> > >
    py_shared_complex_double(
      "cctbx_boost.arraytbx.shared", "complex_double");

    class_builder<fftbx::factorization>
    py_factorization(this_module, "factorization");
    class_builder<fftbx::complex_to_complex<double> >
    py_complex_to_complex(this_module, "complex_to_complex");
    class_builder<fftbx::real_to_complex<double> >
    py_real_to_complex(this_module, "real_to_complex");
    class_builder<fftbx::complex_to_complex_3d<double> >
    py_complex_to_complex_3d(this_module, "complex_to_complex_3d");
    class_builder<fftbx::real_to_complex_3d<double> >
    py_real_to_complex_3d(this_module, "real_to_complex_3d");

    py_complex_to_complex.declare_base(
      py_factorization, boost::python::without_downcast);
    py_real_to_complex.declare_base(
      py_factorization, boost::python::without_downcast);

    this_module.def(adjust_gridding_2, "adjust_gridding");
    this_module.def(adjust_gridding_3, "adjust_gridding");
    this_module.def(adjust_gridding_triple_2, "adjust_gridding_triple");
    this_module.def(adjust_gridding_triple_3, "adjust_gridding_triple");

    py_factorization.def(constructor<>());
    py_factorization.def(constructor<std::size_t, bool>());
    py_factorization.def(&fftbx::factorization::N, "N");
    py_factorization.def(&fftbx::factorization::Factors, "Factors");

    py_complex_to_complex.def(constructor<>());
    py_complex_to_complex.def(constructor<std::size_t>());
    py_complex_to_complex.def(
      &fftbx::complex_to_complex<double>::WA, "WA");
    py_complex_to_complex.def(cc_forward_complex, "forward");
    py_complex_to_complex.def(cc_backward_complex, "backward");
    py_complex_to_complex.def(cc_forward_real, "forward");
    py_complex_to_complex.def(cc_backward_real, "backward");

    py_real_to_complex.def(constructor<>());
    py_real_to_complex.def(constructor<std::size_t>());
    py_real_to_complex.def(&fftbx::real_to_complex<double>::Nreal, "Nreal");
    py_real_to_complex.def(&fftbx::real_to_complex<double>::Mreal, "Mreal");
    py_real_to_complex.def(
      &fftbx::real_to_complex<double>::Ncomplex, "Ncomplex");
    py_real_to_complex.def(&fftbx::real_to_complex<double>::WA, "WA");
    py_real_to_complex.def(rc_forward_complex, "forward");
    py_real_to_complex.def(rc_backward_complex, "backward");
    py_real_to_complex.def(rc_forward_real, "forward");
    py_real_to_complex.def(rc_backward_real, "backward");

    py_complex_to_complex_3d.def(constructor<>());
    py_complex_to_complex_3d.def(
      constructor<std::size_t, std::size_t, std::size_t>());
    py_complex_to_complex_3d.def(
      constructor<const af::int3&>());
    py_complex_to_complex_3d.def(
      &fftbx::complex_to_complex_3d<double>::N, "N");
    py_complex_to_complex_3d.def(cc_3d_forward_complex, "forward");
    py_complex_to_complex_3d.def(cc_3d_backward_complex, "backward");
    py_complex_to_complex_3d.def(cc_3d_forward_real, "forward");
    py_complex_to_complex_3d.def(cc_3d_backward_real, "backward");

    py_real_to_complex_3d.def(constructor<>());
    py_real_to_complex_3d.def(
      constructor<std::size_t, std::size_t, std::size_t>());
    py_real_to_complex_3d.def(
      constructor<const af::int3&>());
    py_real_to_complex_3d.def(
      &fftbx::real_to_complex_3d<double>::Nreal, "Nreal");
    py_real_to_complex_3d.def(
      &fftbx::real_to_complex_3d<double>::Mreal, "Mreal");
    py_real_to_complex_3d.def(
      &fftbx::real_to_complex_3d<double>::Ncomplex, "Ncomplex");
    py_real_to_complex_3d.def(rc_3d_forward_complex, "forward");
    py_real_to_complex_3d.def(rc_3d_backward_complex, "backward");
    py_real_to_complex_3d.def(rc_3d_forward_real, "forward");
    py_real_to_complex_3d.def(rc_3d_backward_real, "backward");
  }

}

BOOST_PYTHON_MODULE_INIT(fftbx)
{
  boost::python::module_builder this_module("fftbx");
  init_module(this_module);
}
